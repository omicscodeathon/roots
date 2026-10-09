#!/usr/bin/env python3
"""ROOTS trajectory analysis: RMSD, RMSF, Rg, H-bonds, SASA, PCA + comparative plots."""
import os, sys, glob, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

warnings.filterwarnings('ignore')
sns.set_context('talk')
sns.set_style('whitegrid')

BASE = Path('/root/roots')
OUT  = BASE / 'results' / 'analysis'
OUT.mkdir(parents=True, exist_ok=True)
(OUT / 'comparative').mkdir(exist_ok=True)
SYSTEMS = ['Sys1', 'Sys2', 'Sys3', 'Sys4']
LIG = 'UNL'
START_NS = 34  # equilibration cutoff for part0003 trajectories

all_stats = []

for sys_name in SYSTEMS:
    print(f"\n{'='*60}\n{sys_name}\n{'='*60}", flush=True)
    sdir = OUT / sys_name
    sdir.mkdir(exist_ok=True)
    tpr = BASE / 'harvest' / sys_name / 'first_frame.pdb'
    xtc = sorted((BASE / 'harvest' / sys_name).glob('*.xtc'))[-1]
    print(f"TPR: {tpr}\nXTC: {xtc}", flush=True)

    u = mda.Universe(str(tpr), str(xtc))
    print(f"Loaded: {len(u.trajectory)} frames, {u.atoms.n_atoms} atoms", flush=True)

    protein = u.select_atoms('protein')
    ligand  = u.select_atoms(f'resname {LIG}')
    backbone = u.select_atoms('protein and backbone')
    calpha  = u.select_atoms('protein and name CA')
    print(f"Protein: {protein.n_atoms}, Ligand: {ligand.n_atoms}", flush=True)

    # --- Align to first frame ---
    ref = mda.Universe(str(tpr), str(xtc))
    align.AlignTraj(u, ref, select='protein and name CA', in_memory=True).run()

    # --- RMSD ---
    R = rms.RMSD(u, ref,
                 select='protein and name CA',
                 groupselections=['protein and backbone', f'resname {LIG}']).run()
    rmsd = pd.DataFrame(R.results.rmsd,
                        columns=['frame', 'time_ps', 'calpha', 'backbone', 'ligand'])
    rmsd['time_ns'] = rmsd['time_ps'] / 1000.0
    rmsd.to_csv(sdir / 'rmsd.csv', index=False)
    fig, ax = plt.subplots(figsize=(9, 5))
    for col, lab in [('calpha', 'C-alpha'), ('backbone', 'Backbone'), ('ligand', 'Ligand')]:
        ax.plot(rmsd['time_ns'], rmsd[col], label=lab, lw=1.2, alpha=0.85)
    ax.axvline(START_NS, ls='--', color='gray', alpha=0.5, label=f'Equil. ({START_NS}ns)')
    ax.set_xlabel('Time (ns)'); ax.set_ylabel('RMSD (Å)'); ax.set_title(f'{sys_name} RMSD')
    ax.legend(); plt.tight_layout()
    plt.savefig(sdir / 'rmsd.png', dpi=300); plt.close()

    equil = rmsd[rmsd['time_ns'] >= START_NS]
    rmsd_bb_mean = equil['backbone'].mean(); rmsd_bb_std = equil['backbone'].std()
    rmsd_lig_mean = equil['ligand'].mean();  rmsd_lig_std = equil['ligand'].std()
    print(f"  Backbone RMSD: {rmsd_bb_mean:.2f} +/- {rmsd_bb_std:.2f} A", flush=True)
    print(f"  Ligand RMSD:   {rmsd_lig_mean:.2f} +/- {rmsd_lig_std:.2f} A", flush=True)

    # --- RMSF (per-residue, equilibrated portion) ---
    ca_group = u.select_atoms('protein and name CA')
    rmsf_vals = []
    # select frames after START_NS
    equil_frames = [ts.frame for ts in u.trajectory if ts.time / 1000.0 >= START_NS]
    if equil_frames:
        aln = align.AlignTraj(u, ref, select='protein and name CA', in_memory=True)
        aln.run(start=equil_frames[0])
        R2 = rms.RMSF(ca_group).run(start=equil_frames[0])
        rmsf_df = pd.DataFrame({'resid': ca_group.resids, 'resname': ca_group.resnames,
                                'rmsf': R2.results.rmsf})
        rmsf_df.to_csv(sdir / 'rmsf.csv', index=False)
        fig, ax = plt.subplots(figsize=(12, 5))
        ax.plot(rmsf_df['resid'], rmsf_df['rmsf'], lw=1.2, color='steelblue')
        ax.fill_between(rmsf_df['resid'], 0, rmsf_df['rmsf'], alpha=0.3, color='steelblue')
        ax.set_xlabel('Residue'); ax.set_ylabel('RMSF (Å)')
        ax.set_title(f'{sys_name} per-residue RMSF (C-alpha, from {START_NS}ns)')
        plt.tight_layout(); plt.savefig(sdir / 'rmsf.png', dpi=300); plt.close()
        print(f"  Mean RMSF: {rmsf_df['rmsf'].mean():.2f} A", flush=True)

    # --- Radius of gyration ---
    rgs, times = [], []
    for ts in u.trajectory:
        rgs.append(protein.radius_of_gyration())
        times.append(ts.time / 1000.0)
    rg = pd.DataFrame({'time_ns': times, 'rg': rgs})
    rg.to_csv(sdir / 'rg.csv', index=False)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(rg['time_ns'], rg['rg'], lw=1.0, color='teal')
    ax.axvline(START_NS, ls='--', color='gray', alpha=0.5)
    ax.set_xlabel('Time (ns)'); ax.set_ylabel('Rg (Å)'); ax.set_title(f'{sys_name} Radius of Gyration')
    plt.tight_layout(); plt.savefig(sdir / 'rg.png', dpi=300); plt.close()
    rg_mean = rg[rg['time_ns'] >= START_NS]['rg'].mean()
    print(f"  Rg: {rg_mean:.2f} A", flush=True)

    # --- H-bonds (protein-ligand) ---
    try:
        hba = HBA(universe=u,
                  donors_sel=f'protein',
                  acceptors_sel=f'resname {LIG}',
                  d_h_cutoff=1.2, d_a_cutoff=3.5, d_h_a_angle_cutoff=150)
        hba.run()
        hbonds = hba.results.hbonds
        hb_count_per_frame = []
        for ts in u.trajectory:
            if ts.time / 1000.0 >= START_NS:
                mask = hbonds[:, 0] == ts.frame
                hb_count_per_frame.append(mask.sum())
        hb_mean = np.mean(hb_count_per_frame) if hb_count_per_frame else 0
        print(f"  Mean H-bonds (protein->ligand): {hb_mean:.2f}", flush=True)

        # reverse direction
        hba2 = HBA(universe=u, donors_sel=f'resname {LIG}',
                   acceptors_sel='protein',
                   d_h_cutoff=1.2, d_a_cutoff=3.5, d_h_a_angle_cutoff=150)
        hba2.run()
        hbonds2 = hba2.results.hbonds
        hb2_count = []
        for ts in u.trajectory:
            if ts.time / 1000.0 >= START_NS:
                hb2_count.append((hbonds2[:, 0] == ts.frame).sum())
        hb2_mean = np.mean(hb2_count) if hb2_count else 0
        print(f"  Mean H-bonds (ligand->protein): {hb2_mean:.2f}", flush=True)
        hb_total = hb_mean + hb2_mean
    except Exception as e:
        print(f"  H-bond analysis failed: {e}", flush=True)
        hb_total = np.nan

    # --- PCA (first 2 PCs of C-alpha) ---
    try:
        from MDAnalysis.analysis import pca
        p = pca.PCA(u, select='protein and name CA', align=False).run()
        print(f"  PC1 variance: {p.results.variance[0]:.2f}%, PC2: {p.results.variance[1]:.2f}%", flush=True)
        transformed = p.transform(ca_group, n_components=2)
        fig, ax = plt.subplots(figsize=(7, 6))
        sc = ax.scatter(transformed[:, 0], transformed[:, 1],
                        c=np.arange(len(transformed)), cmap='viridis', s=15, alpha=0.7)
        plt.colorbar(sc, ax=ax, label='Frame')
        ax.set_xlabel(f'PC1 ({p.results.variance[0]:.1f}%)')
        ax.set_ylabel(f'PC2 ({p.results.variance[1]:.1f}%)')
        ax.set_title(f'{sys_name} PCA (C-alpha)')
        plt.tight_layout(); plt.savefig(sdir / 'pca.png', dpi=300); plt.close()
    except Exception as e:
        print(f"  PCA failed: {e}", flush=True)

    all_stats.append({
        'system': sys_name,
        'backbone_rmsd_mean': rmsd_bb_mean,
        'backbone_rmsd_std': rmsd_bb_std,
        'ligand_rmsd_mean': rmsd_lig_mean,
        'ligand_rmsd_std': rmsd_lig_std,
        'rg_mean': rg_mean,
        'hbonds_mean': hb_total,
        'n_frames': len(u.trajectory),
    })

# --- Comparative plots ---
print("\n=== Comparative plots ===", flush=True)
summary = pd.DataFrame(all_stats)
summary.to_csv(OUT / 'comparative' / 'summary_statistics.csv', index=False)
print(summary.to_string(index=False))

# Comparative RMSD overlay
fig, ax = plt.subplots(figsize=(10, 6))
for s in SYSTEMS:
    df = pd.read_csv(OUT / s / 'rmsd.csv')
    ax.plot(df['time_ns'], df['backbone'], label=s, lw=1.3, alpha=0.85)
ax.axvline(START_NS, ls='--', color='gray', alpha=0.5)
ax.set_xlabel('Time (ns)'); ax.set_ylabel('Backbone RMSD (Å)')
ax.set_title('Backbone RMSD across systems'); ax.legend()
plt.tight_layout(); plt.savefig(OUT / 'comparative' / 'comparative_rmsd.png', dpi=300); plt.close()

# Comparative RMSF
fig, ax = plt.subplots(figsize=(14, 6))
for s in SYSTEMS:
    f = OUT / s / 'rmsf.csv'
    if f.exists():
        df = pd.read_csv(f)
        ax.plot(df['resid'], df['rmsf'], label=s, lw=1.2, alpha=0.85)
ax.set_xlabel('Residue'); ax.set_ylabel('RMSF (Å)')
ax.set_title('Per-residue RMSF across systems'); ax.legend()
plt.tight_layout(); plt.savefig(OUT / 'comparative' / 'comparative_rmsf.png', dpi=300); plt.close()

# Barplot of key metrics
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
summary.plot.bar(x='system', y='backbone_rmsd_mean', yerr='backbone_rmsd_std',
                 ax=axes[0], legend=False, color='steelblue')
axes[0].set_ylabel('Backbone RMSD (Å)'); axes[0].set_title('Protein stability')
summary.plot.bar(x='system', y='ligand_rmsd_mean', yerr='ligand_rmsd_std',
                 ax=axes[1], legend=False, color='coral')
axes[1].set_ylabel('Ligand RMSD (Å)'); axes[1].set_title('Ligand positioning')
summary.plot.bar(x='system', y='hbonds_mean',
                 ax=axes[2], legend=False, color='seagreen')
axes[2].set_ylabel('Avg H-bonds'); axes[2].set_title('Protein-ligand H-bonds')
for ax in axes: ax.set_xlabel('')
plt.tight_layout(); plt.savefig(OUT / 'comparative' / 'metrics_barplot.png', dpi=300); plt.close()

print(f"\n=== ANALYSIS COMPLETE ===\nOutput: {OUT}", flush=True)
