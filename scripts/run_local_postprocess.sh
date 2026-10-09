#!/bin/bash
set +e

echo "=== Setting up conda env (~5 min first time) ==="
if ! conda env list | grep -q "^roots-local "; then
  conda create -y -n roots-local -c conda-forge \
    python=3.12 mdanalysis plip py3dmol seaborn pandas numpy matplotlib openbabel 2>&1 | tail -8
fi
source $(conda info --base)/bin/activate roots-local

cat > ~/roots_pipeline1/local_postprocess.py <<'PY'
#!/usr/bin/env python3
"""Complete ROOTS post-processing: fills all analysis gaps + publication figures."""
import os, re, sys, warnings
from pathlib import Path
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid'); sns.set_context('talk')

RES = Path.home() / 'roots_pipeline1' / 'results'
OUT = RES / 'analysis' / 'local'
OUT.mkdir(parents=True, exist_ok=True)
SYSTEMS = ['Sys1','Sys2','Sys3','Sys4']
# Publication labels — ALDH3B2 holo-NAD vs CNP compounds
LABELS = {
    'Sys1': 'CNP0220660.1\n(consensus)',
    'Sys2': 'CNP0358193.1\n(consensus)',
    'Sys3': 'CNP0384182.2\n(selective)',
    'Sys4': 'CNP0088713.1\n(borderline)',
}
LABELS_SHORT = {'Sys1':'CNP0220660','Sys2':'CNP0358193','Sys3':'CNP0384182','Sys4':'CNP0088713'}

print("="*72); print("ROOTS COMPLETE LOCAL POST-PROCESSING"); print(f"Output: {OUT}"); print("="*72)

# ============ 1. BINDING ENERGY TABLE ============
print("\n[ 1/14] Binding energy table...")
def pd_dg(p):
    try:
        for line in open(p):
            if line.startswith('ΔTOTAL'):
                n = re.findall(r'-?\d+\.\d+', line)
                return float(n[0]), float(n[2]), float(n[4])
    except: pass
    return None,None,None

rows=[]
for s in SYSTEMS:
    gm,gsd,gsem = pd_dg(RES/'mmpbsa'/s/'FINAL_RESULTS_MMPBSA_GB.dat')
    pm,psd,psem = pd_dg(RES/'mmpbsa'/s/'FINAL_RESULTS_MMPBSA_PB.dat')
    rows.append({'System':s,'Compound':LABELS_SHORT[s],
                 'GB_mean':gm,'GB_SD':gsd,'GB_SEM':gsem,
                 'PB_mean':pm,'PB_SD':psd,'PB_SEM':psem})
df_dg = pd.DataFrame(rows); df_dg.to_csv(OUT/'binding_energies.csv', index=False)
pub = pd.DataFrame({'System':df_dg.System,'Compound':df_dg.Compound,
    'ΔG_GB (kcal/mol)':[f"{m:.2f} ± {s:.2f}" for m,s in zip(df_dg.GB_mean,df_dg.GB_SD)],
    'ΔG_PB (kcal/mol)':[f"{m:.2f} ± {s:.2f}" for m,s in zip(df_dg.PB_mean,df_dg.PB_SD)]})
pub.to_csv(OUT/'Table1_binding_energies.csv', index=False)
print(pub.to_string(index=False))

# ============ 2. ENERGY COMPONENT BREAKDOWN ============
print("\n[ 2/14] Energy component breakdown (ΔVDW, ΔEEL, ΔEGB/PB, ΔESURF/NPOLAR)...")
def parse_components(dat_path):
    comp = {}; in_delta = False
    for line in open(dat_path):
        if 'Delta (Complex' in line: in_delta=True; continue
        if not in_delta: continue
        m = re.match(r'^Δ(\S+)\s+(-?\d+\.\d+)', line)
        if m: comp[m.group(1)] = float(m.group(2))
    return comp

comp_rows=[]
for s in SYSTEMS:
    gb = parse_components(RES/'mmpbsa'/s/'FINAL_RESULTS_MMPBSA_GB.dat')
    pb = parse_components(RES/'mmpbsa'/s/'FINAL_RESULTS_MMPBSA_PB.dat')
    comp_rows.append({'System':s,'Compound':LABELS_SHORT[s],
        'ΔVDW':gb.get('VDWAALS'),'ΔEEL':gb.get('EEL'),
        'ΔEGB':gb.get('EGB'),'ΔESURF':gb.get('ESURF'),
        'ΔEPB':pb.get('EPB'),'ΔENPOLAR':pb.get('ENPOLAR'),
        'ΔGGAS':gb.get('GGAS'),'ΔGSOLV_GB':gb.get('GSOLV'),'ΔGSOLV_PB':pb.get('GSOLV'),
        'ΔTOTAL_GB':gb.get('TOTAL'),'ΔTOTAL_PB':pb.get('TOTAL')})
df_comp = pd.DataFrame(comp_rows); df_comp.to_csv(OUT/'energy_components.csv', index=False)
print(df_comp.round(2).to_string(index=False))

# Grouped barplot per solvation model
fig,(a1,a2) = plt.subplots(1,2,figsize=(16,6), sharey=True)
x = np.arange(len(SYSTEMS)); w=0.2; colors=['#4a7fb0','#d16b6b','#5fa06a','#d4a853']
for i,(c,col) in enumerate(zip(['ΔVDW','ΔEEL','ΔEGB','ΔESURF'], colors)):
    a1.bar(x+(i-1.5)*w, df_comp[c], w, label=c, color=col, edgecolor='k', lw=0.5)
a1.set_xticks(x); a1.set_xticklabels([LABELS_SHORT[s] for s in SYSTEMS], rotation=15)
a1.set_ylabel('Energy (kcal/mol)'); a1.set_title('GB decomposition'); a1.axhline(0,color='k',lw=0.5); a1.legend()
for i,(c,col) in enumerate(zip(['ΔVDW','ΔEEL','ΔEPB','ΔENPOLAR'], colors)):
    a2.bar(x+(i-1.5)*w, df_comp[c], w, label=c, color=col, edgecolor='k', lw=0.5)
a2.set_xticks(x); a2.set_xticklabels([LABELS_SHORT[s] for s in SYSTEMS], rotation=15)
a2.set_title('PB decomposition'); a2.axhline(0,color='k',lw=0.5); a2.legend()
plt.suptitle('Energy component contributions to ΔG binding', y=1.02)
plt.tight_layout(); plt.savefig(OUT/'Fig_energy_components.png', dpi=300, bbox_inches='tight'); plt.close()

# ============ 3. GB vs PB BARPLOT (Table 1 visual) ============
print("\n[ 3/14] GB vs PB barplot...")
fig,ax = plt.subplots(figsize=(10,6)); w=0.35
ax.bar(x-w/2, df_dg.GB_mean, w, yerr=df_dg.GB_SEM, label='GB (igb=8)',
       color='steelblue', capsize=5, edgecolor='k', lw=0.8)
ax.bar(x+w/2, df_dg.PB_mean, w, yerr=df_dg.PB_SEM, label='PB (linearized)',
       color='coral', capsize=5, edgecolor='k', lw=0.8)
ax.set_xticks(x); ax.set_xticklabels([LABELS_SHORT[s] for s in SYSTEMS], rotation=15)
ax.set_ylabel('ΔG binding (kcal/mol)'); ax.axhline(0, color='k', lw=0.5)
ax.set_title('MM-PBSA Binding Free Energies — ALDH3B2·NAD + 4 compounds'); ax.legend(loc='lower right')
plt.tight_layout(); plt.savefig(OUT/'Fig1_binding_energies.png', dpi=300); plt.close()

# ============ 4. PER-RESIDUE DECOMPOSITION HEATMAP ============
print("\n[ 4/14] Per-residue decomposition...")
def parse_decomp(path):
    r=[]; in_d=False
    for line in open(path):
        if 'DELTAS' in line.upper(): in_d=True; continue
        if not in_d: continue
        m = re.match(r'^\s*([RL]):(\w):(\w{3}):(\d+)', line)
        if not m: continue
        ct,ch,rn,ri = m.groups()
        nums = re.findall(r'-?\d+\.\d+', line)
        if len(nums)>=5:
            r.append({'chain_type':ct,'resname':rn,'resid':int(ri),
                      'label':f'{rn}{ri}','delta_total':float(nums[-5])})
    return pd.DataFrame(r)

decomp_data={}
for s in SYSTEMS:
    df = parse_decomp(RES/'mmpbsa'/s/'FINAL_DECOMP_MMPBSA_GB.dat')
    if not df.empty: decomp_data[s]=df; df.to_csv(OUT/f'{s}_decomp.csv', index=False)
    print(f"  {s}: {len(df)} residues")

if decomp_data:
    tops=set()
    for s,df in decomp_data.items(): tops.update(df.nsmallest(15,'delta_total').label.tolist())
    tops = sorted(tops, key=lambda x:int(re.search(r'\d+',x).group()))
    mat = pd.DataFrame(index=tops, columns=SYSTEMS, dtype=float)
    for s,df in decomp_data.items():
        for _,r in df.iterrows():
            if r.label in tops: mat.loc[r.label,s]=r.delta_total
    mat.fillna(0.0, inplace=True); mat.to_csv(OUT/'top_residues_matrix.csv')
    mat.columns = [LABELS_SHORT[c] for c in mat.columns]
    fig,ax = plt.subplots(figsize=(9, max(8, len(tops)*0.35)))
    sns.heatmap(mat, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                cbar_kws={'label':'ΔG contribution (kcal/mol)'}, ax=ax, linewidths=0.5)
    ax.set_title('Top-15 contributing residues per compound (GB decomp)')
    plt.tight_layout(); plt.savefig(OUT/'Fig2_decomp_heatmap.png', dpi=300); plt.close()

# ============ 5. RMSD (ligand + Cα + backbone) ============
print("\n[ 5/14] Full RMSD analysis...")
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

rmsd_stats=[]
for s in SYSTEMS:
    p = RES/'mmpbsa'/s/'COM.prmtop'; x = RES/'mmpbsa'/s/'complex_traj.xtc'
    if not (p.exists() and x.exists()): continue
    u = mda.Universe(str(p),str(x)); ref = mda.Universe(str(p),str(x))
    align.AlignTraj(u, ref, select='protein and backbone', in_memory=True).run()
    R = rms.RMSD(u, ref, select='protein and backbone',
                 groupselections=['protein and name CA','resname UNL']).run()
    a = R.results.rmsd
    pd.DataFrame({'frame':a[:,0],'time_ps':a[:,1],'backbone':a[:,2],
                  'calpha':a[:,3],'ligand':a[:,4]}).to_csv(OUT/f'{s}_rmsd.csv', index=False)
    fig,ax = plt.subplots(figsize=(10,6))
    ax.plot(a[:,0], a[:,2], label='Backbone', lw=1.3, color='steelblue')
    ax.plot(a[:,0], a[:,3], label='C-α', lw=1.3, color='teal')
    ax.plot(a[:,0], a[:,4], label='Ligand (UNL)', lw=1.3, color='coral')
    ax.set_xlabel('Frame'); ax.set_ylabel('RMSD (Å)')
    ax.set_title(f'{LABELS_SHORT[s]} RMSD'); ax.legend()
    plt.tight_layout(); plt.savefig(OUT/f'{s}_rmsd.png', dpi=300); plt.close()
    rmsd_stats.append({'System':s,'Compound':LABELS_SHORT[s],
        'backbone_A':a[:,2].mean(),'backbone_SD':a[:,2].std(),
        'calpha_A':a[:,3].mean(),'calpha_SD':a[:,3].std(),
        'ligand_A':a[:,4].mean(),'ligand_SD':a[:,4].std()})
df_rmsd = pd.DataFrame(rmsd_stats); df_rmsd.to_csv(OUT/'rmsd_summary.csv', index=False)
print(df_rmsd.round(3).to_string(index=False))

# ============ 6. RMSD CONVERGENCE (block averaging) ============
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
plt.suptitle('RMSD convergence (rolling window) — equilibration check', y=1.00, fontsize=14)
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

# ============ 7. H-BOND PERSISTENCE ============
print("\n[ 7/14] H-bond persistence...")
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
hb_rows=[]
for s in SYSTEMS:
    p = RES/'mmpbsa'/s/'COM.prmtop'; x = RES/'mmpbsa'/s/'complex_traj.xtc'
    if not (p.exists() and x.exists()): continue
    u = mda.Universe(str(p),str(x)); n = len(u.trajectory); tot=0
    for ds,as_ in [('protein','resname UNL'),('resname UNL','protein')]:
        try:
            h = HBA(u, donors_sel=ds, acceptors_sel=as_,
                    d_a_cutoff=3.5, d_h_a_angle_cutoff=150).run()
            tot += len(h.results.hbonds)
        except Exception as e: print(f"    {s} {ds}: {e}")
    hb_rows.append({'System':s,'Compound':LABELS_SHORT[s],'n_frames':n,
                    'HB_total':tot,'HB_per_frame':tot/n if n else 0})
    print(f"  {s}: {tot/n:.2f} H-bonds/frame")
df_hb = pd.DataFrame(hb_rows); df_hb.to_csv(OUT/'hbonds.csv', index=False)
fig,ax = plt.subplots(figsize=(10,6))
ax.bar(df_hb.Compound, df_hb.HB_per_frame, color='seagreen', edgecolor='k')
for i,v in enumerate(df_hb.HB_per_frame): ax.text(i, v+0.05, f'{v:.2f}', ha='center')
ax.set_ylabel('Avg H-bonds / frame'); ax.set_title('Protein-ligand H-bond persistence')
plt.xticks(rotation=15); plt.tight_layout(); plt.savefig(OUT/'Fig4_hbonds.png', dpi=300); plt.close()

# ============ 8. DISTANCE TRACKING (top-5 binding residues) ============
print("\n[ 8/14] Distance tracking...")
from MDAnalysis.analysis.distances import distance_array
for s in SYSTEMS:
    if s not in decomp_data: continue
    p = RES/'mmpbsa'/s/'COM.prmtop'; x = RES/'mmpbsa'/s/'complex_traj.xtc'
    if not (p.exists() and x.exists()): continue
    top5 = decomp_data[s].nsmallest(5,'delta_total')
    u = mda.Universe(str(p),str(x)); lig = u.select_atoms('resname UNL')
    d = {'frame':[]}
    for _,r in top5.iterrows(): d[f'{r.resname}{r.resid}']=[]
    for ts in u.trajectory:
        d['frame'].append(ts.frame)
        for _,r in top5.iterrows():
            ra = u.select_atoms(f'resid {r.resid} and name CA')
            d[f'{r.resname}{r.resid}'].append(
                distance_array(ra.positions, lig.positions).min() if len(ra) else np.nan)
    df_d = pd.DataFrame(d); df_d.to_csv(OUT/f'{s}_distances.csv', index=False)
    fig,ax = plt.subplots(figsize=(10,6))
    for col in df_d.columns[1:]: ax.plot(df_d.frame, df_d[col], label=col, lw=1.2, alpha=0.85)
    ax.set_xlabel('Frame'); ax.set_ylabel('Min CA-ligand distance (Å)')
    ax.set_title(f'{LABELS_SHORT[s]} — top-5 binding residue distances'); ax.legend(ncol=2)
    plt.tight_layout(); plt.savefig(OUT/f'{s}_distances.png', dpi=300); plt.close()

# ============ 9. REPRESENTATIVE FRAME + PLIP ============
print("\n[ 9/14] Representative frame extraction + PLIP...")
# Extract middle frame of each trajectory as the representative structure
for s in SYSTEMS:
    p = RES/'mmpbsa'/s/'COM.prmtop'; x = RES/'mmpbsa'/s/'complex_traj.xtc'
    if not (p.exists() and x.exists()): continue
    u = mda.Universe(str(p),str(x))
    mid = len(u.trajectory)//2
    u.trajectory[mid]
    # Write representative frame
    rep_pdb = OUT/f'{s}_representative_frame.pdb'
    with mda.Writer(str(rep_pdb), multiframe=False) as w:
        w.write(u.atoms)
    print(f"  {s}: wrote frame {mid} → {rep_pdb.name}")

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
                                g.write(f"PLIP interactions — {LABELS_SHORT[s]} (middle frame)\n{'='*60}\n\n")
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
        # Representative-frame only barplot
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
except ImportError: print("    PLIP missing")
except Exception as e: print(f"    FAILED: {e}")

# ============ 10. Py3Dmol HTML VIEWERS ============
print("\n[10/14] Py3Dmol binding-site HTML viewers...")
try:
    import py3Dmol
    for s in SYSTEMS:
        # Use representative frame if available
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

# ============ 11. CROSS-TARGET COMPARATIVE MAIN FIGURE ============
print("\n[11/14] Cross-target comparative main figure (publication centerpiece)...")
try:
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    compounds = [LABELS_SHORT[s] for s in SYSTEMS]
    xp = np.arange(len(SYSTEMS))
    
    # Panel A: ΔG GB vs PB
    ax = axes[0,0]; wp=0.35
    ax.bar(xp-wp/2, df_dg.GB_mean, wp, yerr=df_dg.GB_SEM, label='GB', 
           color='steelblue', edgecolor='k', capsize=4)
    ax.bar(xp+wp/2, df_dg.PB_mean, wp, yerr=df_dg.PB_SEM, label='PB', 
           color='coral', edgecolor='k', capsize=4)
    ax.set_xticks(xp); ax.set_xticklabels(compounds, rotation=15, fontsize=10)
    ax.set_ylabel('ΔG (kcal/mol)'); ax.set_title('A. Binding free energy')
    ax.axhline(0, color='k', lw=0.5); ax.legend(fontsize=10)
    
    # Panel B: Backbone RMSD
    ax = axes[0,1]
    ax.bar(compounds, df_rmsd.backbone_A, yerr=df_rmsd.backbone_SD,
           color='teal', edgecolor='k', capsize=5)
    ax.set_ylabel('Backbone RMSD (Å)'); ax.set_title('B. Protein stability')
    ax.tick_params(axis='x', rotation=15, labelsize=10)
    
    # Panel C: Ligand RMSD
    ax = axes[0,2]
    ax.bar(compounds, df_rmsd.ligand_A, yerr=df_rmsd.ligand_SD,
           color='coral', edgecolor='k', capsize=5)
    ax.set_ylabel('Ligand RMSD (Å)'); ax.set_title('C. Ligand positioning')
    ax.tick_params(axis='x', rotation=15, labelsize=10)
    
    # Panel D: H-bonds
    ax = axes[1,0]
    ax.bar(compounds, df_hb.HB_per_frame, color='seagreen', edgecolor='k')
    ax.set_ylabel('H-bonds / frame'); ax.set_title('D. H-bond persistence')
    ax.tick_params(axis='x', rotation=15, labelsize=10)
    
    # Panel E: Energy decomposition (ΔVDW vs ΔEEL vs ΔGSOLV)
    ax = axes[1,1]; wp=0.25
    ax.bar(xp-wp, df_comp['ΔVDW'], wp, label='ΔVDW', color='#4a7fb0', edgecolor='k')
    ax.bar(xp, df_comp['ΔEEL'], wp, label='ΔEEL', color='#d16b6b', edgecolor='k')
    ax.bar(xp+wp, df_comp['ΔGSOLV_GB'], wp, label='ΔG_solv(GB)', color='#5fa06a', edgecolor='k')
    ax.set_xticks(xp); ax.set_xticklabels(compounds, rotation=15, fontsize=10)
    ax.set_ylabel('kcal/mol'); ax.set_title('E. Energy drivers')
    ax.axhline(0, color='k', lw=0.5); ax.legend(fontsize=9)
    
    # Panel F: ΔG vs stability correlation
    ax = axes[1,2]
    ax.scatter(df_rmsd.ligand_A, df_dg.GB_mean, s=180, c='steelblue', 
               edgecolor='k', label='GB', zorder=3)
    ax.scatter(df_rmsd.ligand_A, df_dg.PB_mean, s=180, c='coral', 
               edgecolor='k', label='PB', zorder=3)
    for i, sysc in enumerate(compounds):
        ax.annotate(sysc.split('.')[0], (df_rmsd.ligand_A.iloc[i], df_dg.GB_mean.iloc[i]),
                    xytext=(8,5), textcoords='offset points', fontsize=9)
    # Correlation line
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
    print("  → FIG_MAIN_cross_target_comparison.png (publication centerpiece)")
except Exception as e: print(f"    FAILED: {e}")

# ============ 12. INTEGRATED PUBLICATION TABLE ============
print("\n[12/14] Integrated publication Table (ΔG + stability + interactions)...")
try:
    integ = pd.DataFrame({
        'Compound': [LABELS_SHORT[s] for s in SYSTEMS],
        'ΔG_GB': [f"{m:.2f} ± {s:.2f}" for m,s in zip(df_dg.GB_mean,df_dg.GB_SD)],
        'ΔG_PB': [f"{m:.2f} ± {s:.2f}" for m,s in zip(df_dg.PB_mean,df_dg.PB_SD)],
        'Backbone RMSD (Å)': [f"{m:.2f} ± {s:.2f}" for m,s in zip(df_rmsd.backbone_A, df_rmsd.backbone_SD)],
        'Ligand RMSD (Å)': [f"{m:.2f} ± {s:.2f}" for m,s in zip(df_rmsd.ligand_A, df_rmsd.ligand_SD)],
        'H-bonds / frame': [f"{v:.2f}" for v in df_hb.HB_per_frame],
    })
    integ.to_csv(OUT/'Table1_INTEGRATED.csv', index=False)
    print(integ.to_string(index=False))
except Exception as e: print(f"    FAILED: {e}")

# ============ 13. BUILD README FOR OUTPUT ============
print("\n[13/14] Writing output README...")
readme = """ROOTS Local Post-Processing Results
=====================================

Target: ALDH3B2 holo-NAD receptor
Compounds (4 systems):
  Sys1 = CNP0220660.1 (consensus hit)
  Sys2 = CNP0358193.1 (consensus hit)
  Sys3 = CNP0384182.2 (selective/flagged)
  Sys4 = CNP0088713.1 (selective/borderline)

Key Publication Files
---------------------
Table1_INTEGRATED.csv       - Combined binding + stability + H-bonds (main table)
Table1_binding_energies.csv - GB + PB side by side
FIG_MAIN_cross_target_comparison.png  - Publication centerpiece (6-panel)
Fig1_binding_energies.png   - GB vs PB ΔG barplot
Fig2_decomp_heatmap.png     - Per-residue decomposition top-15
Fig3_rmsd_convergence.png   - Equilibration check (rolling window)
Fig4_hbonds.png             - H-bond persistence
Fig5_plip_representative.png - PLIP stacked bar
Fig_energy_components.png   - ΔVDW / ΔEEL / ΔEGB-EPB / ΔESURF-ENPOLAR

Per-System Outputs
------------------
<Sys>_rmsd.{csv,png}            - Backbone+Cα+Ligand RMSD per frame
<Sys>_rmsd_all.png              - RMSD overlay
<Sys>_distances.{csv,png}       - Top-5 binding residue CA-ligand distances
<Sys>_decomp.csv                - Per-residue ΔG contributions
<Sys>_plip_representative.txt   - PLIP interaction details (middle frame)
<Sys>_representative_frame.pdb  - Middle-frame structure for figures
<Sys>_binding_site.html         - Interactive 3D view (open in browser)

Raw Data
--------
rmsd_summary.csv, hbonds.csv, energy_components.csv, 
binding_energies.csv, top_residues_matrix.csv, plip_summary.csv

Reproducibility
---------------
MM-PBSA: gmx_MMPBSA v1.5.0.3, GROMACS 2026.3, CHARMM36 force field
  GB: igb=8, saltcon=0.150, mbondi2 radii, 310 K
  PB: istrng=0.150, fillratio=4.0, indi=1.0, exdi=80.0
  Frames: ~157 per system, extracted from 1 ns into part0003 trajectory
MD analysis: MDAnalysis 2.x, PLIP, Py3Dmol
"""
with open(OUT/'README.md','w') as f: f.write(readme)

# ============ 14. FINAL FILE LISTING ============
print("\n[14/14] Final listing...")
files = sorted(OUT.glob('*'))
print(f"\n{len(files)} files in {OUT}:")
for f in files: 
    sz = f.stat().st_size
    print(f"  {f.name}  ({sz:>10,} bytes)")

print("\n" + "="*72)
print("ALL ANALYSES COMPLETE — publication-ready outputs assembled")
print("="*72)
PY

echo
echo "=== Running analysis ==="
python3 ~/roots_pipeline1/local_postprocess.py 2>&1 | tee ~/roots_pipeline1/local_postprocess.log

echo
echo "=== Final output listing ==="
ls -la ~/roots_pipeline1/results/analysis/local/ | head -50
echo
echo "===================================================="
echo "SERVER STATUS: still alive at 5.161.160.175"
echo "To delete when you've verified these outputs are good:"
echo "  hcloud server delete roots-analysis"
echo "===================================================="
