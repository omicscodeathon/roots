#!/usr/bin/env python3
import fnmatch, importlib, os, sys, tarfile
from datetime import datetime
from pathlib import Path

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else ".").expanduser().resolve()
OUT = ROOT / "PreAdmetDiagnostics"

TARGET_PATTERNS = [
    ("*hit_reconciliation*.csv", "hit reconciliation (Vina/PRODIGY consensus)"),
    ("*vina_scores*.csv",        "extracted Vina scores"),
    ("*prodigy*scores*.csv",     "merged PRODIGY-LIG scores"),
    ("*manifest*.csv",           "library manifest (compound_id -> smiles)"),
    ("*library*.csv",            "library / SMILES table (fallback for manifest)"),
    ("*hit_smiles*.csv",         "hit compound_id -> SMILES extract"),
    ("*admet*predict*.csv",      "ADMET-AI raw predictions"),
    ("*tier2_shortlist*",        "tier-2 shortlist for pkCSM/SwissADME/ADMETlab"),
]
PACKAGES = ["admet_ai", "rdkit", "meeko", "Bio", "pandas", "gemmi"]


def human_size(n):
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}PB"


def scan(root):
    files = []
    for dirpath, _, filenames in os.walk(root):
        if Path(dirpath).name == "__pycache__":
            continue
        for fn in filenames:
            fp = Path(dirpath) / fn
            try:
                st = fp.stat()
            except OSError:
                continue
            files.append((fp, st.st_size, st.st_mtime))
    return files


def top_level_summary(root, files):
    buckets = {}
    for fp, size, mtime in files:
        try:
            rel = fp.relative_to(root)
        except ValueError:
            continue
        top = rel.parts[0] if len(rel.parts) > 1 else "(root)"
        b = buckets.setdefault(top, {"n": 0, "size": 0, "newest": 0})
        b["n"] += 1
        b["size"] += size
        b["newest"] = max(b["newest"], mtime)
    lines = []
    for name in sorted(buckets):
        b = buckets[name]
        newest = datetime.fromtimestamp(b["newest"]).strftime("%Y-%m-%d %H:%M")
        label = f"{name}/" if name != "(root)" else "(loose files in root)"
        lines.append(f"  {label:<30} {b['n']:>7} files, {human_size(b['size']):>10}, newest {newest}")
    return lines


def csv_rows_and_header(path, max_bytes=200_000_000):
    try:
        size = path.stat().st_size
        with open(path, "r", newline="", errors="replace") as f:
            header = f.readline().rstrip("\n")
            if size > max_bytes:
                return None, header
            n = sum(1 for _ in f)
        return n + 1, header
    except Exception as e:
        return None, f"<unreadable: {e}>"


def find_targets(root, files):
    lines = []
    for pattern, label in TARGET_PATTERNS:
        matches = [fp for fp, _, _ in files if fnmatch.fnmatch(fp.name, pattern)]
        if not matches:
            lines.append(f"  [MISSING] {label}  (pattern: {pattern})")
            continue
        for m in sorted(matches):
            rel = m.relative_to(root)
            if m.suffix.lower() == ".csv":
                n_rows, header = csv_rows_and_header(m)
                rows_str = f"{n_rows} rows incl. header" if n_rows is not None else "too large to count quickly"
                lines.append(f"  [FOUND]   {label}: {rel}  ({rows_str})")
                lines.append(f"            header: {header[:200]}")
            else:
                lines.append(f"  [FOUND]   {label}: {rel}  ({human_size(m.stat().st_size)})")
    return lines


def docking_archives(root):
    lines = []
    server_dirs = sorted(set(list(root.glob("roots-docking-*")) + list(root.glob("*/roots-docking-*"))))
    if not server_dirs:
        return ["  none found under this root"]
    for sd in server_dirs:
        archives = list(sd.glob("results_*.tar.gz"))
        if not archives:
            lines.append(f"  {sd.name}: no results_*.tar.gz present")
            continue
        for arc in archives:
            try:
                with tarfile.open(arc) as t:
                    n_docked = sum(1 for m in t.getmembers() if m.name.endswith("_docked.pdbqt"))
                lines.append(f"  {sd.name}/{arc.name}: {n_docked} docked poses, {human_size(arc.stat().st_size)}")
            except Exception as e:
                lines.append(f"  {sd.name}/{arc.name}: unreadable ({e})")
    return lines


def package_check():
    lines = [f"  python: {sys.version.split()[0]} ({sys.executable})"]
    for pkg in PACKAGES:
        try:
            mod = importlib.import_module(pkg)
            ver = getattr(mod, "__version__", "version n/a")
            lines.append(f"  {pkg}: installed ({ver})")
        except ImportError:
            lines.append(f"  {pkg}: NOT INSTALLED")
    return lines


def main():
    files = scan(ROOT)
    report = []
    report.append("=== ROOTS pipeline -- pre-ADMET diagnostic ===")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"Root: {ROOT}")
    report.append(f"Total files under root: {len(files)}, total size: {human_size(sum(s for _, s, _ in files))}")
    report.append("")
    report.append("--- Top-level directory summary ---")
    report.extend(top_level_summary(ROOT, files))
    report.append("")
    report.append("--- Key artifacts for the ADMET-AI step ---")
    report.extend(find_targets(ROOT, files))
    report.append("")
    report.append("--- Docking-fleet archives (roots-docking-* dirs) ---")
    report.extend(docking_archives(ROOT))
    report.append("")
    report.append("--- Python environment (packages needed downstream) ---")
    report.extend(package_check())
    report.append("")
    report.append("=== End of report ===")

    text = "\n".join(report)
    print(text)
    OUT.write_text(text + "\n")
    print(f"\nWritten to: {OUT}")


if __name__ == "__main__":
    main()
