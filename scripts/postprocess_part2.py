#!/usr/bin/env python3
"""Continuation: trajectory analyses via mdtraj (handles CHARMM-style prmtop)."""
import os, re, warnings
from pathlib import Path
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid'); sns.set_context('talk')

RES = Path.home() / 'roots_pipeline1' / 'results'
OUT = RES / 'analysis' / 'local'
SYSTEMS = ['Sys1','Sys2','Sys3','Sys4']
LABELS_SHORT = {'Sys1':'CNP0220660','Sys2':'CNP0358193','Sys3':'CNP0384182','Sys4':'CNP0088713'}

df_dg = pd.read_csv(OUT/'binding_energies.csv')
df_comp = pd.read_csv(OUT/'energy_components.csv')
decomp_data = {s: pd.read_csv(OUT/f'{s}_decomp.csv') for s in SYSTEMS if (OUT/f'{s}_decomp.csv').exists()}

import mdtraj as md

# ============ 5. RMSD (backbone + Cα + ligand) via mdtraj ============
print("\n[ 5/14] RMSD via mdtraj...")
rmsd_stats=[]
trajs={}
for s in SYSTEMS:
    pt = RES/'mmpbsa'/s/'COM.prmtop'; xt = RES/'mmpbsa'/s/'complex_traj.xtc'
    if not (pt.exists() and xt.exists()): continue
    try:
        traj = md.load(str(xt), top=str(pt))
        trajs[s] = traj
        print(f"  {s}: {traj.n_frames} frames, {traj.n_atoms} atoms")
        bb_idx = traj.topology.select('backbone')
        ca_idx = traj.topology.select('name CA')
        lig_idx = traj.topology.select('resname UNL')
        traj.superpose(traj[0], atom_indices=bb_idx)
        bb_rmsd = md.rmsd(traj, traj[0], atom_indices=bb_idx) * 10  # nm to Å
        ca_rmsd = md.rmsd(traj, traj[0], atom_indices=ca_idx) * 10
        lig_rmsd = md.rmsd(traj, traj[0], atom_indices=lig_idx) * 10
        frames = np.arange(traj.n_frames)
        pd.DataFrame({'frame':frames,'time_ps':traj.time,
                      'backbone':bb_rmsd,'calpha':ca_rmsd,'ligand':lig_rmsd}
                    ).to_csv(OUT/f'{s}_rmsd.csv', index=False)
        fig,ax = plt.subplots(figsize=(10,6))
        ax.plot(frames, bb_rmsd, label='Backbone', lw=1.3, color='steelblue')
        ax.plot(frames, ca_rmsd, label='C-α', lw=1.3, color='teal')
        ax.plot(frames, lig_rmsd, label='Ligand (UNL)', lw=1.3, color='coral')
        ax.set_xlabel('Frame'); ax.set_ylabel('RMSD (Å)')
        ax.set_title(f'{LABELS_SHORT[s]} RMSD'); ax.legend()
        plt.tight_layout(); plt.savefig(OUT/f'{s}_rmsd.png', dpi=300); plt.close()
        rmsd_stats.append({'System':s,'Compound':LABELS_SHORT[s],
            'backbone_A':bb_rmsd.mean(),'backbone_SD':bb_rmsd.std(),
            'calpha_A':ca_rmsd.mean(),'calpha_SD':ca_rmsd.std(),
            'ligand_A':lig_rmsd.mean(),'ligand_SD':lig_rmsd.std()})
        print(f"    bb={bb_rmsd.mean():.2f}±{bb_rmsd.std():.2f}Å, lig={lig_rmsd.mean():.2f}±{lig_rmsd.std():.2f}Å")
    except Exception as e:
        print(f"  {s}: FAILED {e}")
df_rmsd = pd.DataFrame(rmsd_stats); df_rmsd.to_csv(OUT/'rmsd_summary.csv', index=False)
print(df_rmsd.round(3).to_string(index=False))

# ============ 6. RMSD CONVERGENCE ============
print("\n[ 6/14] RMSD convergence check...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
for ax, s in zip(axes.flat, SYSTEMS):
    f = OUT/f'{s}_rmsd.csv'
    if not f.exists(): continue
    df = pd.read_csv(f); win = max(10, len(df)//10)
    rm = df.backbone.rolling(win, min_periods=1).mean()
    rs = df.backbone.rolling(win, min_periods=1).std()
    ax.plot(df.frame, df.backbone, lw=0.8, alpha=0.4, color='steelblue', label='Instantaneous')
    ax.plot(df.frame, rm, lw=2, color='darkblue', label=f'Rolling mean (w={win})')
    ax.fill_between(df.frame, rm-rs, rm+rs, alpha=0.2, color='steelblue')
    ax.set_xlabel('Frame'); ax.set_ylabel('Backbone RMSD (Å)')
    ax.set_title(f'{LABELS_SHORT[s]} convergence'); ax.legend(fontsize=10)
plt.suptitle('RMSD convergence (rolling window)', y=1.00, fontsize=14)
plt.tight_layout(); plt.savefig(OUT/'Fig3_rmsd_convergence.png', dpi=300); plt.close()

# Comparative ligand RMSD
fig,ax = plt.subplots(figsize=(10,6))
for s in SYSTEMS:
    f = OUT/f'{s}_rmsd.csv'
    if f.exists():
        df = pd.read_csv(f); ax.plot(df.frame, df.ligand, label=LABELS_SHORT[s], lw=1.3, alpha=0.85)
ax.set_xlabel('Frame'); ax.set_ylabel('Ligand RMSD (Å)')
ax.set_title('Ligand RMSD across compounds'); ax.legend()
plt.tight_layout(); plt.savefig(OUT/'Fig_comparative_ligand_rmsd.png', dpi=300); plt.close()

# ============ 7. H-BONDS via mdtraj baker_hubbard ============
print("\n[ 7/14] H-bond persistence via mdtraj...")
hb_rows=[]
for s in SYSTEMS:
    if s not in trajs: continue
    traj = trajs[s]
    try:
        # baker_hubbard returns (donor, H, acceptor) atom indices
        hbonds = md.baker_hubbard(traj, periodic=False, freq=0.0, 
                                    distance_cutoff=0.35, angle_cutoff=150)
        prot_idx = set(traj.topology.select('protein'))
        lig_idx = set(traj.topology.select('resname UNL'))
        # Count protein-ligand H-bonds (either direction)
        pl_hb = 0
        for d, h, a in hbonds:
            if (d in prot_idx and a in lig_idx) or (d in lig_idx and a in prot_idx):
                pl_hb += 1
        # Need per-frame count - rerun with freq=0 trick doesn't give per-frame
        # Compute per-frame H-bonds via geometric criterion
        from scipy.spatial.distance import cdist
        prot_no = traj.topology.select('protein and (symbol N or symbol O)')
        lig_no = traj.topology.select('resname UNL and (symbol N or symbol O)')
        per_frame = []
        for i in range(traj.n_frames):
            prot_pos = traj.xyz[i, prot_no] * 10  # nm to Å
            lig_pos = traj.xyz[i, lig_no] * 10
            d = cdist(prot_pos, lig_pos)
            per_frame.append((d < 3.5).sum())
        hb_rows.append({'System':s,'Compound':LABELS_SHORT[s],
                        'n_frames':traj.n_frames,
                        'baker_hubbard_total':pl_hb,
                        'HB_per_frame':np.mean(per_frame)})
        print(f"  {s}: {np.mean(per_frame):.2f} N/O-N/O pairs <3.5Å per frame")
    except Exception as e:
        print(f"  {s}: FAILED {e}")
df_hb = pd.DataFrame(hb_rows); df_hb.to_csv(OUT/'hbonds.csv', index=False)
fig,ax = plt.subplots(figsize=(10,6))
ax.bar(df_hb.Compound, df_hb.HB_per_frame, color='seagreen', edgecolor='k')
for i,v in enumerate(df_hb.HB_per_frame): ax.text(i, v+0.05, f'{v:.2f}', ha='center')
ax.set_ylabel('Avg H-bond-capable N/O pairs < 3.5Å'); ax.set_title('Protein-ligand H-bond persistence')
plt.xticks(rotation=15); plt.tight_layout(); plt.savefig(OUT/'Fig4_hbonds.png', dpi=300); plt.close()

# ============ 8. DISTANCE TRACKING (top-5 residues) ============
print("\n[ 8/14] Distance tracking...")
from scipy.spatial.distance import cdist
for s in SYSTEMS:
    if s not in decomp_data or s not in trajs: continue
    traj = trajs[s]
    top5 = decomp_data[s].nsmallest(5,'delta_total')
    lig_idx = traj.topology.select('resname UNL')
    d = {'frame': np.arange(traj.n_frames)}
    for _,r in top5.iterrows():
        ca_sel = traj.topology.select(f'resid {r.resid-1} and name CA')  # mdtraj uses 0-indexed
        if len(ca_sel) == 0:
            ca_sel = traj.topology.select(f'residue {r.resid} and name CA')
        label = f'{r.resname}{r.resid}'
        if len(ca_sel) == 0:
            d[label] = [np.nan]*traj.n_frames
            continue
        dists = []
        for i in range(traj.n_frames):
            ca_pos = traj.xyz[i, ca_sel] * 10
            lig_pos = traj.xyz[i, lig_idx] * 10
            dists.append(cdist(ca_pos, lig_pos).min())
        d[label] = dists
    df_d = pd.DataFrame(d); df_d.to_csv(OUT/f'{s}_distances.csv', index=False)
    fig,ax = plt.subplots(figsize=(10,6))
    for col in df_d.columns[1:]: ax.plot(df_d.frame, df_d[col], label=col, lw=1.2, alpha=0.85)
    ax.set_xlabel('Frame'); ax.set_ylabel('Min CA-ligand distance (Å)')
    ax.set_title(f'{LABELS_SHORT[s]} — top-5 binding residue distances'); ax.legend(ncol=2)
    plt.tight_layout(); plt.savefig(OUT/f'{s}_distances.png', dpi=300); plt.close()
    print(f"  {s}: tracked {len(top5)} residues")

# ============ 9. REPRESENTATIVE FRAME + PLIP ============
print("\n[ 9/14] Representative frame extraction + PLIP...")
for s in SYSTEMS:
    if s not in trajs: continue
    traj = trajs[s]; mid = traj.n_frames//2
    traj[mid].save_pdb(str(OUT/f'{s}_representative_frame.pdb'))
    print(f"  {s}: wrote frame {mid}")

try:
    from plip.structure.preparation import PDBComplex
    plip_rows=[]
    for s in SYSTEMS:
        for pdb_label, pdb_path in [('first', RES/'mmpbsa'/s/'_GMXMMPBSA_COM.pdb'),
                                     ('rep',   OUT/f'{s}_representative_frame.pdb')]:
            if not pdb_path.exists(): continue
            try:
                mol = PDBComplex(); mol.load_pdb(str(pdb_path)); mol.analyze()
                for lk, inter in mol.interaction_sets.items():
                    if 'UNL' in lk:
                        plip_rows.append({'System':s,'Compound':LABELS_SHORT[s],'Frame':pdb_label,
                            'HB': len(inter.hbonds_pdon)+len(inter.hbonds_ldon),
                            'Hydrophobic': len(inter.hydrophobic_contacts),
                            'SaltBridge': len(inter.saltbridge_lneg)+len(inter.saltbridge_pneg),
                            'Pi_stack': len(inter.pistacking),
                            'Pi_cation': len(inter.pication_laro)+len(inter.pication_paro),
                            'Halogen': len(inter.halogen_bonds),
                            'WaterBridge': len(inter.water_bridges)})
                        if pdb_label == 'rep':
                            with open(OUT/f'{s}_plip_representative.txt','w') as g:
                                g.write(f"PLIP — {LABELS_SHORT[s]} (middle frame)\n{'='*60}\n\n")
                                g.write("H-BONDS (protein donor):\n")
                                for hb in inter.hbonds_pdon:
                                    g.write(f"  {hb.restype}{hb.resnr} dist={hb.distance_ad:.2f}Å\n")
                                g.write("\nH-BONDS (ligand donor):\n")
                                for hb in inter.hbonds_ldon:
                                    g.write(f"  {hb.restype}{hb.resnr} dist={hb.distance_ad:.2f}Å\n")
                                g.write("\nHYDROPHOBIC:\n")
                                for hc in inter.hydrophobic_contacts:
                                    g.write(f"  {hc.restype}{hc.resnr} dist={hc.distance:.2f}Å\n")
                                g.write("\nSALT BRIDGES:\n")
                                for sb in inter.saltbridge_lneg + inter.saltbridge_pneg:
                                    g.write(f"  {sb.restype}{sb.resnr} dist={sb.distance:.2f}Å\n")
                                g.write("\nπ-STACKING:\n")
                                for ps in inter.pistacking:
                                    g.write(f"  {ps.restype}{ps.resnr} dist={ps.distance:.2f}Å\n")
                        break
            except Exception as ee: print(f"    {s} ({pdb_label}): {ee}")
    df_plip = pd.DataFrame(plip_rows)
    if not df_plip.empty:
        df_plip.to_csv(OUT/'plip_summary.csv', index=False)
        print(df_plip.to_string(index=False))
        rep_df = df_plip[df_plip.Frame=='rep'].reset_index(drop=True)
        if not rep_df.empty:
            itypes = ['HB','Hydrophobic','SaltBridge','Pi_stack','Pi_cation','Halogen','WaterBridge']
            cs = sns.color_palette('Set2', len(itypes))
            fig,ax = plt.subplots(figsize=(11,6)); bot = np.zeros(len(rep_df))
            for i,t in enumerate(itypes):
                if t in rep_df:
                    ax.bar(rep_df.Compound, rep_df[t], bottom=bot, label=t,
                           color=cs[i], edgecolor='k', lw=0.5)
                    bot += rep_df[t].values
            ax.set_ylabel('# interactions'); ax.set_title('PLIP interactions (representative frame)')
            ax.legend(loc='upper left', bbox_to_anchor=(1,1))
            plt.xticks(rotation=15); plt.tight_layout()
            plt.savefig(OUT/'Fig5_plip_representative.png', dpi=300, bbox_inches='tight'); plt.close()
except Exception as e: print(f"    PLIP: {e}")

# ============ 10. Py3Dmol VIEWERS ============
print("\n[10/14] Py3Dmol HTML viewers...")
try:
    import py3Dmol
    for s in SYSTEMS:
        pdb = OUT/f'{s}_representative_frame.pdb'
        if not pdb.exists(): pdb = RES/'mmpbsa'/s/'_GMXMMPBSA_COM.pdb'
        if not pdb.exists(): continue
        view = py3Dmol.view(width=900, height=650)
        view.addModel(open(pdb).read(), 'pdb')
        view.setStyle({}, {'cartoon':{'color':'lightgray'}})
        view.setStyle({'resn':'UNL'}, {'stick':{'colorscheme':'yellowCarbon','radius':0.3}})
        view.addStyle({'byres':True,'within':{'distance':5,'sel':{'resn':'UNL'}}},
                      {'stick':{'colorscheme':'cyanCarbon','radius':0.15}})
        view.zoomTo({'resn':'UNL'}); view.spin(True)
        with open(OUT/f'{s}_binding_site.html','w') as g: g.write(view._make_html())
        print(f"  → {s}_binding_site.html")
except Exception as e: print(f"    FAILED: {e}")

# ============ 11. CROSS-TARGET COMPARATIVE FIGURE ============
print("\n[11/14] Cross-target comparative figure...")
try:
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    compounds = [LABELS_SHORT[s] for s in SYSTEMS]
    xp = np.arange(len(SYSTEMS))
    # A: ΔG
    ax = axes[0,0]; wp=0.35
    ax.bar(xp-wp/2, df_dg.GB_mean, wp, yerr=df_dg.GB_SEM, label='GB', 
           color='steelblue', edgecolor='k', capsize=4)
    ax.bar(xp+wp/2, df_dg.PB_mean, wp, yerr=df_dg.PB_SEM, label='PB', 
           color='coral', edgecolor='k', capsize=4)
    ax.set_xticks(xp); ax.set_xticklabels(compounds, rotation=15, fontsize=10)
    ax.set_ylabel('ΔG (kcal/mol)'); ax.set_title('A. Binding free energy')
    ax.axhline(0, color='k', lw=0.5); ax.legend(fontsize=10)
    # B: BB RMSD
    ax = axes[0,1]
    ax.bar(compounds, df_rmsd.backbone_A, yerr=df_rmsd.backbone_SD,
           color='teal', edgecolor='k', capsize=5)
    ax.set_ylabel('Backbone RMSD (Å)'); ax.set_title('B. Protein stability')
    ax.tick_params(axis='x', rotation=15, labelsize=10)
    # C: Lig RMSD
    ax = axes[0,2]
    ax.bar(compounds, df_rmsd.ligand_A, yerr=df_rmsd.ligand_SD,
           color='coral', edgecolor='k', capsize=5)
    ax.set_ylabel('Ligand RMSD (Å)'); ax.set_title('C. Ligand positioning')
    ax.tick_params(axis='x', rotation=15, labelsize=10)
    # D: H-bonds
    ax = axes[1,0]
    ax.bar(compounds, df_hb.HB_per_frame, color='seagreen', edgecolor='k')
    ax.set_ylabel('H-bond pairs / frame'); ax.set_title('D. H-bond persistence')
    ax.tick_params(axis='x', rotation=15, labelsize=10)
    # E: Energy drivers
    ax = axes[1,1]; wp=0.25
    ax.bar(xp-wp, df_comp['ΔVDW'], wp, label='ΔVDW', color='#4a7fb0', edgecolor='k')
    ax.bar(xp, df_comp['ΔEEL'], wp, label='ΔEEL', color='#d16b6b', edgecolor='k')
    ax.bar(xp+wp, df_comp['ΔGSOLV_GB'], wp, label='ΔG_solv(GB)', color='#5fa06a', edgecolor='k')
    ax.set_xticks(xp); ax.set_xticklabels(compounds, rotation=15, fontsize=10)
    ax.set_ylabel('kcal/mol'); ax.set_title('E. Energy drivers')
    ax.axhline(0, color='k', lw=0.5); ax.legend(fontsize=9)
    # F: Correlation
    ax = axes[1,2]
    ax.scatter(df_rmsd.ligand_A, df_dg.GB_mean, s=180, c='steelblue', 
               edgecolor='k', label='GB', zorder=3)
    ax.scatter(df_rmsd.ligand_A, df_dg.PB_mean, s=180, c='coral', 
               edgecolor='k', label='PB', zorder=3)
    for i, sysc in enumerate(compounds):
        ax.annotate(sysc, (df_rmsd.ligand_A.iloc[i], df_dg.GB_mean.iloc[i]),
                    xytext=(8,5), textcoords='offset points', fontsize=9)
    try:
        from scipy.stats import pearsonr
        r,p = pearsonr(df_rmsd.ligand_A, df_dg.GB_mean)
        ax.text(0.05, 0.95, f'Pearson r = {r:.2f}\np = {p:.3f}', transform=ax.transAxes,
                verticalalignment='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    except: pass
    ax.set_xlabel('Ligand RMSD (Å)'); ax.set_ylabel('ΔG binding (kcal/mol)')
    ax.set_title('F. Binding–stability correlation'); ax.legend(fontsize=10)
    plt.suptitle('ROOTS: ALDH3B2·NAD — cross-compound binding and dynamics summary',
                 fontsize=15, y=1.01)
    plt.tight_layout()
    plt.savefig(OUT/'FIG_MAIN_cross_target_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  → FIG_MAIN_cross_target_comparison.png")
except Exception as e: print(f"    FAILED: {e}")

# ============ 12. INTEGRATED TABLE ============
print("\n[12/14] Integrated publication Table...")
try:
    integ = pd.DataFrame({
        'Compound': [LABELS_SHORT[s] for s in SYSTEMS],
        'ΔG_GB': [f"{m:.2f} ± {s:.2f}" for m,s in zip(df_dg.GB_mean,df_dg.GB_SD)],
        'ΔG_PB': [f"{m:.2f} ± {s:.2f}" for m,s in zip(df_dg.PB_mean,df_dg.PB_SD)],
        'Backbone RMSD (Å)': [f"{m:.2f} ± {s:.2f}" for m,s in zip(df_rmsd.backbone_A, df_rmsd.backbone_SD)],
        'Ligand RMSD (Å)': [f"{m:.2f} ± {s:.2f}" for m,s in zip(df_rmsd.ligand_A, df_rmsd.ligand_SD)],
        'H-bond pairs/frame': [f"{v:.2f}" for v in df_hb.HB_per_frame],
    })
    integ.to_csv(OUT/'Table1_INTEGRATED.csv', index=False)
    print(integ.to_string(index=False))
except Exception as e: print(f"    FAILED: {e}")

# ============ 13. README ============
print("\n[13/14] Writing README...")
readme = """ROOTS Local Post-Processing — ALDH3B2·NAD study
==================================================
Compounds:
  Sys1 = CNP0220660.1 (consensus)
  Sys2 = CNP0358193.1 (consensus)
  Sys3 = CNP0384182.2 (selective/flagged)
  Sys4 = CNP0088713.1 (selective/borderline)

Key Publication Files
---------------------
Table1_INTEGRATED.csv              - ΔG + stability + H-bonds (main table)
Table1_binding_energies.csv        - GB + PB side by side
FIG_MAIN_cross_target_comparison.png - 6-panel publication centerpiece
Fig1_binding_energies.png          - GB vs PB barplot
Fig2_decomp_heatmap.png            - Top-15 residue contributions
Fig3_rmsd_convergence.png          - Rolling-window equilibration check
Fig4_hbonds.png                    - H-bond persistence
Fig5_plip_representative.png       - PLIP interaction stacked bar
Fig_energy_components.png          - ΔVDW/ΔEEL/ΔEGB-EPB/ΔESURF-ENPOLAR
Fig_comparative_ligand_rmsd.png    - Ligand RMSD overlay

Per-system
----------
<Sys>_rmsd.{csv,png}               - Backbone+Cα+Ligand RMSD
<Sys>_distances.{csv,png}          - Top-5 binding residue distances
<Sys>_decomp.csv                   - Per-residue ΔG contributions
<Sys>_representative_frame.pdb     - Middle-frame structure
<Sys>_plip_representative.txt      - PLIP details
<Sys>_binding_site.html            - Interactive 3D viewer

Methods
-------
MM-PBSA: gmx_MMPBSA v1.5.0.3 (Valdes-Tresanco et al. JCTC 2021)
  GROMACS 2026.3, CHARMM36, GB (igb=8) & PB (linearized)
  157 frames per system, 310K, saltcon=0.150 M
MD analysis: mdtraj 1.x (prmtop CHARMM-compat), PLIP 2.x, Py3Dmol
H-bond proxy: geometric N/O-N/O pairs < 3.5Å (per frame avg)
"""
with open(OUT/'README.md','w') as f: f.write(readme)

# ============ 14. FINAL LISTING ============
print("\n[14/14] Final listing...")
files = sorted(OUT.glob('*'))
print(f"\n{len(files)} files in {OUT}:")
for f in files: print(f"  {f.name}  ({f.stat().st_size:>10,} bytes)")
print("\n" + "="*72); print("COMPLETE"); print("="*72)
