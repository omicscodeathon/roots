import subprocess
import sys
import os
import csv
import time
import json
import multiprocessing as mp
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError
from concurrent.futures.process import BrokenProcessPool

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")
MMFF_MAX_ITERS = 50

def process_chunk(tasks, output_dir, n_workers, per_task_timeout, chunk_wall_budget):
    pending = {cid: (cid, smi) for cid, smi in tasks}
    results = []
    start = time.time()

    while pending:
        if time.time() - start > chunk_wall_budget:
            for cid in pending:
                results.append((cid, False, "chunk_wall_time_exceeded"))
            break

        executor = ProcessPoolExecutor(max_workers=n_workers)
        futures = {executor.submit(prep_one, t, output_dir): cid for cid, t in pending.items()}
        broken = False

        for f in list(futures):
            cid = futures[f]
            try:
                _, ok, reason = f.result(timeout=per_task_timeout)
                results.append((cid, ok, reason))
                del pending[cid]
            except TimeoutError:
                results.append((cid, False, "task_timeout"))
                del pending[cid]
            except BrokenProcessPool:
                broken = True
                break
            except Exception as e:
                results.append((cid, False, f"crashed: {type(e).__name__}"))
                del pending[cid]

        executor.shutdown(wait=False, cancel_futures=True)
        if not broken:
            break

    return results


def prep_one(args, output_dir):
    compound_id, smi = args
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return compound_id, False, "smiles_parse_failed"
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) == -1:
            return compound_id, False, "embedding_failed"
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=MMFF_MAX_ITERS)
        except Exception:
            pass
        mk_prep = MoleculePreparation()
        molsetup_list = mk_prep(mol)
        if not molsetup_list:
            return compound_id, False, "meeko_prep_failed"
        write_result = PDBQTWriterLegacy.write_string(molsetup_list[0])
        if isinstance(write_result, tuple):
            pdbqt_string = write_result[0]
            is_ok = write_result[1] if len(write_result) > 1 else True
            if not is_ok or not pdbqt_string:
                return compound_id, False, "pdbqt_write_failed"
        else:
            pdbqt_string = write_result
        if not isinstance(pdbqt_string, str) or not pdbqt_string:
            return compound_id, False, "pdbqt_write_non_string"
        (Path(output_dir) / f"{compound_id}.pdbqt").write_text(pdbqt_string)
        return compound_id, True, None
    except Exception as e:
        return compound_id, False, f"unexpected: {type(e).__name__}: {e}"


def main():
    mp.set_start_method('spawn', force=True)

    BASE_DIR = Path("/workspace")
    STEP3_DIR = BASE_DIR / "step3"
    STEP4_DIR = BASE_DIR / "step4"
    STEP4_DIR.mkdir(parents=True, exist_ok=True)
    (STEP4_DIR / "logs").mkdir(exist_ok=True)

    LIBRARY_PATH = STEP3_DIR / "library_final.csv"
    if not LIBRARY_PATH.exists():
        raise FileNotFoundError(f"{LIBRARY_PATH} not found.")

    OUTPUT_DIR = STEP4_DIR / "ligands"
    OUTPUT_DIR.mkdir(exist_ok=True)
    PROGRESS_PATH = STEP4_DIR / "logs" / "progress.txt"
    FAIL_LOG_PATH = STEP4_DIR / "logs" / "failures.csv"

    # Tunable via env vars (set them in `docker run -e ...`) so box size can
    # change without editing this file. Defaults below match the values this
    # was last tuned for on an 11-worker (12 vCPU-ish) machine.
    CHUNK_SIZE = int(os.environ.get("CHUNK_SIZE", 200))
    N_WORKERS = int(os.environ.get("N_WORKERS", 11))
    CHUNK_TIMEOUT_S = int(os.environ.get("CHUNK_TIMEOUT_S", 400))
    PER_TASK_TIMEOUT_S = int(os.environ.get("PER_TASK_TIMEOUT_S", 45))  # one slow/hung molecule no longer takes the rest of the chunk with it

    print(f"Config: N_WORKERS={N_WORKERS} CHUNK_SIZE={CHUNK_SIZE} "
          f"CHUNK_TIMEOUT_S={CHUNK_TIMEOUT_S} PER_TASK_TIMEOUT_S={PER_TASK_TIMEOUT_S}")

    library_df = pd.read_csv(LIBRARY_PATH, low_memory=False)
    smiles_col = "canonical_smiles"
    id_col = "identifier" if "identifier" in library_df.columns else "id"

    existing_ids = set(f.stem for f in OUTPUT_DIR.glob("*.pdbqt"))
    print(f"Already done (on disk, permanently): {len(existing_ids)}")

    remaining_df = library_df[~library_df[id_col].astype(str).isin(existing_ids)].reset_index(drop=True)
    total = len(remaining_df)
    print(f"Remaining: {total}\n")

    n_success = 0
    n_fail = 0
    fail_log_exists = FAIL_LOG_PATH.exists()

    for start in range(0, total, CHUNK_SIZE):
        end = min(start + CHUNK_SIZE, total)
        chunk = remaining_df.iloc[start:end]
        tasks = list(zip(chunk[id_col].astype(str), chunk[smiles_col]))

        t0 = time.time()
        chunk_failures = []
        chunk_success = 0
        chunk_fail = 0
        timed_out = False

        chunk_results = process_chunk(
            tasks, str(OUTPUT_DIR),
            n_workers=N_WORKERS,
            per_task_timeout=PER_TASK_TIMEOUT_S,
            chunk_wall_budget=CHUNK_TIMEOUT_S,
        )
        for cid, ok, reason in chunk_results:
            if ok:
                chunk_success += 1
            else:
                chunk_fail += 1
                chunk_failures.append((cid, reason))
                if reason in ("task_timeout", "chunk_wall_time_exceeded"):
                    timed_out = True

        n_success += chunk_success
        n_fail += chunk_fail

        if chunk_failures:
            with open(FAIL_LOG_PATH, 'a', newline='') as fh:
                w = csv.writer(fh)
                if not fail_log_exists:
                    w.writerow(['compound_id', 'reason'])
                    fail_log_exists = True
                w.writerows(chunk_failures)

        PROGRESS_PATH.write_text(str(end))
        elapsed = time.time() - t0
        per_mol = elapsed / max(len(chunk), 1)
        total_done = len(existing_ids) + n_success
        eta_min = ((total - end) * per_mol) / 60
        status = "TIMEOUT" if timed_out else "ok"

        print(f"  {end}/{total} [{status}] | {total_done} total on disk | "
              f"{chunk_success} ok, {chunk_fail} failed | {per_mol:.3f}s/mol | ETA: {eta_min:.0f} min")

    print(f"\nDone. Session: {n_success} ok, {n_fail} failed")
    print(f"TOTAL PDBQT files on disk: {len(list(OUTPUT_DIR.glob('*.pdbqt')))}")


if __name__ == "__main__":
    main()
