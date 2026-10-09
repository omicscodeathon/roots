#!/usr/bin/env python3
"""ROOTS analysis via gmx with integer group indices (reliable across versions)."""
import subprocess, numpy as np, pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid'); sns.set_context('talk')

BASE = Path('/root/roots')
OUT  = BASE / 'results' / 'analysis'
OUT.mkdir(parents=True, exist_ok=True); (OUT / 'comparative').mkdir(exist_ok=True)
SYSTEMS = ['Sys1', 'Sys2', 'Sys3', 'Sys4']
GRP_PROTEIN='1'; GRP_CALPHA='3'; GRP_BACKBONE='4'; GRP_UNL='13'

def read_xvg(p):
    if not Path(p).exists(): return None, None
    rows = []
    for line in open(p):
        if line.startswith(('@', '#')): continue
        parts = line.split()
        if parts:
            try: rows.append([float(x) for x in parts])
            except: continue
    if not rows: return None, None
    a = np.array(rows); return a[:, 0], a[:, 1:]

def gmx(args, stdin=None, cwd=None, timeout=900):
    try:
        r = subprocess.run(['gmx'] + args, input=stdin, text=True,
                           capture_output=True, cwd=str(cwd), timeout=timeout)
        if r.returncode != 0:
            err = r.stderr or r.stdout
            for line in err.split('\n'):
                if any(k in line.lower() for k in ['error','fatal','cannot','invalid']):
                    print(f"      {line.strip()}", flush=True)
        return r.returncode == 0
    except Exception as e:
        print(f"    gmx {args[0]} exception: {e}", flush=True); return False

all_stats = []
for s in SYSTEMS:
    print(f"\n=== {s} ===", flush=True)
    sdir = OUT / s; sdir.mkdir(exist_ok=True)
    hdir = BASE / 'harvest' / s
    tpr  = next(hdir.glob('*.tpr'))
    xtc  = sorted(hdir.glob('*.xtc'))[-1]
    ndx  = BASE / 'results' / 'mmpbsa' / s / 'complex.ndx'
    st = {'system': s}

    print("  RMSD backbone...", flush=True)
    if gmx(['rms','-s',str(tpr),'-f',str(xtc),'-n',str(ndx),'-o','rmsd_bb.xvg','-fit','rot+trans'],
           f'{GRP_BACKBONE}\n{GRP_BACKBONE}\n', sdir):
        t,v = read_xvg(sdir/'rmsd_bb.xvg')
        if v is not None:
            st['bb_rmsd_mean']=float(v[:,0].mean()); st['bb_rmsd_std']=float(v[:,0].std())
            print(f"    {st['bb_rmsd_mean']:.3f} ± {st['bb_rmsd_std']:.3f} nm", flush=True)

    print("  RMSD ligand (fit on backbone)...", flush=True)
    if gmx(['rms','-s',str(tpr),'-f',str(xtc),'-n',str(ndx),'-o','rmsd_lig.xvg','-fit','rot+trans'],
           f'{GRP_BACKBONE}\n{GRP_UNL}\n', sdir):
        t,v = read_xvg(sdir/'rmsd_lig.xvg')
        if v is not None:
            st['lig_rmsd_mean']=float(v[:,0].mean()); st['lig_rmsd_std']=float(v[:,0].std())
            print(f"    {st['lig_rmsd_mean']:.3f} ± {st['lig_rmsd_std']:.3f} nm", flush=True)

    print("  RMSF per-residue...", flush=True)
    gmx(['rmsf','-s',str(tpr),'-f',str(xtc),'-n',str(ndx),'-o','rmsf.xvg','-res'],
        f'{GRP_CALPHA}\n', sdir)

    print("  Rg...", flush=True)
    if gmx(['gyrate','-s',str(tpr),'-f',str(xtc),'-n',str(ndx),'-o','rg.xvg'],
           f'{GRP_PROTEIN}\n', sdir):
        t,v = read_xvg(sdir/'rg.xvg')
        if v is not None: st['rg_mean']=float(v[:,0].mean())

    print("  SASA...", flush=True)
    gmx(['sasa','-s',str(tpr),'-f',str(xtc),'-n',str(ndx),'-o','sasa.xvg'],
        f'{GRP_PROTEIN}\n', sdir)

    print("  H-bonds P->UNL (non-fatal)...", flush=True)
    gmx(['hbond','-s',str(tpr),'-f',str(xtc),'-n',str(ndx),'-num','hbnum.xvg'],
        f'{GRP_PROTEIN}\n{GRP_UNL}\n', sdir)
    t,v = read_xvg(sdir/'hbnum.xvg')
    st['hbonds_mean'] = float(v[:,0].mean()) if v is not None else np.nan
    if v is not None: print(f"    {st['hbonds_mean']:.2f} H-bonds", flush=True)

    print("  PCA...", flush=True)
    if gmx(['covar','-s',str(tpr),'-f',str(xtc),'-n',str(ndx),'-o','eigenval.xvg','-v','eigenvec.trr'],
           f'{GRP_CALPHA}\n{GRP_CALPHA}\n', sdir):
        gmx(['anaeig','-s',str(tpr),'-f',str(xtc),'-n',str(ndx),'-v','eigenvec.trr',
             '-2d','pca_2d.xvg','-first','1','-last','2'],
            f'{GRP_CALPHA}\n{GRP_CALPHA}\n', sdir)
    all_stats.append(st)

def plot_ts(xvg,out,title,ylabel):
    t,v = read_xvg(xvg)
    if t is None: return
    fig,ax=plt.subplots(figsize=(9,5))
    ax.plot(t/1000.0, v[:,0], lw=1.1, color='steelblue')
    ax.set_xlabel('Time (ns)'); ax.set_ylabel(ylabel); ax.set_title(title)
    plt.tight_layout(); plt.savefig(out, dpi=300); plt.close()

print("\n=== Plots ===", flush=True)
for s in SYSTEMS:
    sd = OUT/s
    plot_ts(sd/'rmsd_bb.xvg',  sd/'rmsd_backbone.png', f'{s} Backbone RMSD', 'RMSD (nm)')
    plot_ts(sd/'rmsd_lig.xvg', sd/'rmsd_ligand.png',   f'{s} Ligand RMSD',   'RMSD (nm)')
    plot_ts(sd/'rg.xvg',       sd/'rg.png',            f'{s} Radius of Gyration', 'Rg (nm)')
    plot_ts(sd/'sasa.xvg',     sd/'sasa.png',          f'{s} SASA',           'SASA (nm²)')
    plot_ts(sd/'hbnum.xvg',    sd/'hbonds.png',        f'{s} Protein-Ligand H-bonds', 'Count')
    t,v = read_xvg(sd/'rmsf.xvg')
    if t is not None:
        fig,ax=plt.subplots(figsize=(12,5))
        ax.plot(t, v[:,0], lw=1.1, color='steelblue')
        ax.fill_between(t, 0, v[:,0], alpha=0.3, color='steelblue')
        ax.set_xlabel('Residue'); ax.set_ylabel('RMSF (nm)'); ax.set_title(f'{s} RMSF')
        plt.tight_layout(); plt.savefig(sd/'rmsf.png', dpi=300); plt.close()
    t,v = read_xvg(sd/'pca_2d.xvg')
    if t is not None and v is not None and v.shape[1] >= 1:
        fig,ax=plt.subplots(figsize=(7,6))
        sc=ax.scatter(t, v[:,0], c=np.arange(len(t)), cmap='viridis', s=15, alpha=0.75)
        plt.colorbar(sc, ax=ax, label='Frame')
        ax.set_xlabel('PC1 (nm)'); ax.set_ylabel('PC2 (nm)'); ax.set_title(f'{s} PCA')
        plt.tight_layout(); plt.savefig(sd/'pca.png', dpi=300); plt.close()

fig,ax=plt.subplots(figsize=(10,6))
for s in SYSTEMS:
    t,v = read_xvg(OUT/s/'rmsd_bb.xvg')
    if t is not None: ax.plot(t/1000.0, v[:,0], label=s, lw=1.3, alpha=0.85)
ax.set_xlabel('Time (ns)'); ax.set_ylabel('Backbone RMSD (nm)'); ax.set_title('Backbone RMSD'); ax.legend()
plt.tight_layout(); plt.savefig(OUT/'comparative'/'comparative_rmsd.png', dpi=300); plt.close()

fig,ax=plt.subplots(figsize=(14,6))
for s in SYSTEMS:
    t,v = read_xvg(OUT/s/'rmsf.xvg')
    if t is not None: ax.plot(t, v[:,0], label=s, lw=1.2, alpha=0.85)
ax.set_xlabel('Residue'); ax.set_ylabel('RMSF (nm)'); ax.set_title('Per-residue RMSF'); ax.legend()
plt.tight_layout(); plt.savefig(OUT/'comparative'/'comparative_rmsf.png', dpi=300); plt.close()

summary = pd.DataFrame(all_stats)
summary.to_csv(OUT/'comparative'/'summary_statistics.csv', index=False)
print(summary.to_string(index=False), flush=True)
print(f"\n=== DONE === output: {OUT}", flush=True)
