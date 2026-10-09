ROOTS Local Post-Processing — ALDH3B2·NAD study
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
