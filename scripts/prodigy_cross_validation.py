#!/usr/bin/env python3
"""
ROOTS -- PRODIGY-LIG cross-validation of Vina/QuickVina2 docking scores.

For every successfully docked compound, merge the receptor and first docked
ligand pose into one complex structure and score it with PRODIGY-LIG.

This is a cross-validation exercise, not a re-ranking exercise. The objective
is to compare the two scoring methods across the full score distribution.

Important implementation details:
- ProdigyLig receives chain identifiers as plain strings, e.g. ["A", "B:LIG"].
- The ligand is rewritten as HETATM records on chain B with residue name LIG.
- Only the first docked model is retained.
- The receptor is converted from PDBQT to PDB once and reused.
- Results are appended to a CSV so interrupted runs can resume.
"""

import argparse
import csv
import multiprocessing as mp
import random
import subprocess
import sys
import tarfile
import time
from pathlib import Path

from Bio.PDB import PDBParser
from prodigy_lig.prodigy_lig import ProdigyLig

LIGAND_CHAIN = "B"
LIGAND_RESNAME = "LIG"
CUTOFF = 10.5


def obabel_convert(pdbqt_path, pdb_path):
    """
    Convert PDBQT to PDB with OpenBabel.

    Returns (True, None) on success or (False, reason) on failure.
    """
    try:
        result = subprocess.run(
            ["obabel", str(pdbqt_path), "-O", str(pdb_path)],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (subprocess.SubprocessError, OSError) as exc:
        return False, str(exc)[:200]

    if (
        result.returncode != 0
        or not pdb_path.exists()
        or pdb_path.stat().st_size == 0
    ):
        return False, (
            result.stderr.strip() or "OpenBabel produced no output"
        )[:200]

    return True, None


def rewrite_ligand_as_chain(pdb_path, chain_id, resname):
    """
    Rewrite the ligand PDB as HETATM records on a specified chain and residue.

    Retain only the first model. OpenBabel may emit more than one MODEL header
    for a pose, so ENDMDL is used to identify the end of the first model.
    """
    lines = pdb_path.read_text().splitlines()
    atom_lines = []
    first_model_closed = False

    for line in lines:
        if first_model_closed:
            continue

        if line.startswith("ENDMDL"):
            first_model_closed = True
            continue

        if line.startswith(("ATOM  ", "HETATM")):
            # PDB fixed-width fields:
            # residue name = columns 18-20
            # chain ID = column 22
            line = (
                line[:17]
                + f"{resname:<3}"
                + line[20:21]
                + chain_id
                + line[22:]
            )
            line = "HETATM" + line[6:]
            atom_lines.append(line)

    if not atom_lines:
        raise ValueError(
            f"No atom records found in the first model of {pdb_path}. "
            "Check the OpenBabel output."
        )

    pdb_path.write_text("\n".join(atom_lines) + "\nEND\n")


def merge_complex(receptor_pdb_text, ligand_pdb_path, out_path):
    """
    Merge receptor records with the rewritten ligand pose.

    Remove receptor END, TER, CONECT and MASTER records before appending
    the ligand. In particular, an END record before the ligand can cause
    parsers to stop reading before chain B is registered.
    """
    footer_prefixes = ("END", "TER", "CONECT", "MASTER")

    receptor_lines = [
        line
        for line in receptor_pdb_text.splitlines()
        if not line.startswith(footer_prefixes)
    ]

    ligand_lines = [
        line
        for line in ligand_pdb_path.read_text().splitlines()
        if line.startswith("HETATM")
    ]

    if not receptor_lines:
        raise ValueError("The converted receptor PDB contains no usable records.")

    if not ligand_lines:
        raise ValueError("The rewritten ligand PDB contains no HETATM records.")

    out_path.write_text(
        "\n".join(receptor_lines)
        + "\n"
        + "\n".join(ligand_lines)
        + "\nEND\n"
    )


def score_one(task):
    """
    Score one docked pose.

    Returns a result dictionary rather than raising, so an individual bad
    pose does not terminate the entire multiprocessing job.
    """
    compound_id, ligand_pdbqt_path, receptor_pdb_text, receptor_chain, work_dir = task

    ligand_pdb = work_dir / f"{compound_id}_lig.pdb"
    complex_pdb = work_dir / f"{compound_id}_complex.pdb"

    try:
        ok, error = obabel_convert(ligand_pdbqt_path, ligand_pdb)

        if not ok:
            return {
                "compound_id": compound_id,
                "prodigy_dg": None,
                "status": f"obabel_failed: {error}",
            }

        rewrite_ligand_as_chain(
            ligand_pdb,
            LIGAND_CHAIN,
            LIGAND_RESNAME,
        )

        merge_complex(
            receptor_pdb_text,
            ligand_pdb,
            complex_pdb,
        )

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(compound_id, str(complex_pdb))

        chains = [
            receptor_chain,
            f"{LIGAND_CHAIN}:{LIGAND_RESNAME}",
        ]

        predictor = ProdigyLig(
            structure,
            chains,
            electrostatics=None,
            cutoff=CUTOFF,
        )
        predictor.predict()

        return {
            "compound_id": compound_id,
            "prodigy_dg": predictor.dg,
            "status": "ok",
        }

    except RuntimeWarning as exc:
        return {
            "compound_id": compound_id,
            "prodigy_dg": None,
            "status": f"no_contacts: {exc}"[:150],
        }

    except Exception as exc:
        return {
            "compound_id": compound_id,
            "prodigy_dg": None,
            "status": f"error: {type(exc).__name__}: {exc}"[:150],
        }

    finally:
        ligand_pdb.unlink(missing_ok=True)
        complex_pdb.unlink(missing_ok=True)


def extract_docked_poses(results_dir, extract_dir):
    """
    Extract docked poses from roots-docking-* result archives.

    Each archive is extracted once, using a marker file so the extraction
    can be resumed without repeatedly extracting previously processed archives.
    """
    results_dir = Path(results_dir)
    extract_dir = Path(extract_dir)
    extract_dir.mkdir(parents=True, exist_ok=True)

    for server_dir in sorted(results_dir.glob("roots-docking-*")):
        marker = extract_dir / f".{server_dir.name}.extracted"

        if marker.exists():
            continue

        archive = next(server_dir.glob("results_*.tar.gz"), None)

        if archive is None:
            print(
                f"  WARNING: no results archive found in {server_dir}; skipping."
            )
            continue

        print(f"  Extracting {archive.name}...")

        try:
            with tarfile.open(archive, "r:gz") as tar:
                # Reject absolute paths and parent-directory traversal.
                safe_members = []

                for member in tar.getmembers():
                    member_path = Path(member.name)

                    if member_path.is_absolute() or ".." in member_path.parts:
                        print(
                            f"  WARNING: skipping unsafe archive member: "
                            f"{member.name}"
                        )
                        continue

                    safe_members.append(member)

                tar.extractall(extract_dir, members=safe_members)

            marker.touch()

        except (tarfile.TarError, OSError) as exc:
            print(f"  WARNING: extraction failed for {archive}: {exc}")
            continue

    poses = sorted(
        (extract_dir / "docking_results").glob("*_docked.pdbqt")
    )

    return poses


def main():
    parser = argparse.ArgumentParser(
        description="Score docked ROOTS compounds with PRODIGY-LIG."
    )

    parser.add_argument(
        "--receptor-pdbqt",
        required=True,
        help="Path to the prepared receptor PDBQT.",
    )
    parser.add_argument(
        "--results-dir",
        required=True,
        help="Directory containing roots-docking-*/results_*.tar.gz archives.",
    )
    parser.add_argument(
        "--extract-dir",
        default="prodigy_work/poses",
        help="Directory in which docked poses are extracted.",
    )
    parser.add_argument(
        "--out",
        default="prodigy_scores.csv",
        help="CSV path for PRODIGY-LIG results.",
    )
    parser.add_argument(
        "--receptor-chain",
        default="A",
        help="Chain identifier assigned to the receptor in the converted PDB.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (mp.cpu_count() or 2) - 1),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Score a random sample of this size; 0 means score all remaining poses.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed used when --limit is set.",
    )

    args = parser.parse_args()

    if args.workers < 1:
        sys.exit("--workers must be at least 1.")

    if args.limit < 0:
        sys.exit("--limit cannot be negative.")

    receptor_pdbqt = Path(args.receptor_pdbqt).expanduser()

    if not receptor_pdbqt.is_file():
        sys.exit(f"Receptor not found: {receptor_pdbqt}")

    results_dir = Path(args.results_dir).expanduser()

    if not results_dir.is_dir():
        sys.exit(f"Results directory not found: {results_dir}")

    work_dir = Path("prodigy_work")
    work_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdb = work_dir / "receptor.pdb"

    if not receptor_pdb.exists() or receptor_pdb.stat().st_size == 0:
        print("Converting receptor PDBQT to PDB (once)...")

        ok, error = obabel_convert(receptor_pdbqt, receptor_pdb)

        if not ok:
            sys.exit(f"Receptor conversion failed: {error}")

    receptor_pdb_text = receptor_pdb.read_text()

    # Basic diagnostic only: receptor chain assignment can depend on the
    # source receptor and the way OpenBabel converted it.
    receptor_chain_present = any(
        line.startswith(("ATOM  ", "HETATM"))
        and len(line) > 21
        and line[21] == args.receptor_chain
        for line in receptor_pdb_text.splitlines()
    )

    if not receptor_chain_present:
        print(
            f"WARNING: receptor chain '{args.receptor_chain}' was not found "
            "in the converted receptor PDB. Check --receptor-chain if "
            "poses fail with chain errors."
        )

    print("Extracting docked poses from the pulled fleet results (resumable)...")

    poses = extract_docked_poses(
        results_dir,
        Path(args.extract_dir),
    )

    print(f"Total docked poses found: {len(poses):,}")

    out_path = Path(args.out).expanduser()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    done_ids = set()

    if out_path.exists() and out_path.stat().st_size > 0:
        with out_path.open(newline="") as file:
            reader = csv.DictReader(file)

            if reader.fieldnames and "compound_id" in reader.fieldnames:
                done_ids = {
                    row["compound_id"]
                    for row in reader
                    if row.get("compound_id")
                }
            else:
                sys.exit(
                    f"Existing output file {out_path} does not have the "
                    "expected compound_id CSV header. Back it up and inspect "
                    "it before resuming."
                )

        print(f"Resuming: {len(done_ids):,} already recorded; skipping those.")

    todo = [
        pose for pose in poses
        if pose.stem.removesuffix("_docked") not in done_ids
    ]

    if args.limit and len(todo) > args.limit:
        todo = random.Random(args.seed).sample(todo, args.limit)
        print(
            f"PILOT: random sample of {args.limit:,} poses "
            f"(seed {args.seed})."
        )

    print(f"Remaining to score: {len(todo):,} | workers: {args.workers}")

    if not todo:
        print("Nothing to do.")
        return

    tasks = [
        (
            pose.stem.removesuffix("_docked"),
            pose,
            receptor_pdb_text,
            args.receptor_chain,
            work_dir,
        )
        for pose in todo
    ]

    new_file = not out_path.exists() or out_path.stat().st_size == 0
    counts = {}
    start_time = time.time()

    context = mp.get_context("spawn")

    with out_path.open("a", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=["compound_id", "prodigy_dg", "status"],
        )

        if new_file:
            writer.writeheader()
            file.flush()

        with context.Pool(args.workers) as pool:
            for index, result in enumerate(
                pool.imap_unordered(score_one, tasks, chunksize=25),
                1,
            ):
                writer.writerow(result)

                status_group = result["status"].split(":")[0]
                counts[status_group] = counts.get(status_group, 0) + 1

                if index % 2000 == 0 or index == len(tasks):
                    file.flush()

                    elapsed = max(time.time() - start_time, 1e-9)
                    rate = index / elapsed
                    eta_minutes = (len(tasks) - index) / rate / 60

                    print(
                        f"  {index:,}/{len(tasks):,} "
                        f"| {rate:.1f}/s "
                        f"| ETA {eta_minutes:.0f} min "
                        f"| {counts}",
                        flush=True,
                    )

    elapsed_minutes = (time.time() - start_time) / 60

    print(f"\nDone in {elapsed_minutes:.1f} min.")
    print(f"Final tallies: {counts}")
    print(f"Output: {out_path}")


if __name__ == "__main__":
    main()