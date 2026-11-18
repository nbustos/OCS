# resources

Public reference atlases and derived summaries used to choose seeds. No
participant data lives here.

## atlases/CIT168_Pauli2018

CIT168 reinforcement-learning subcortical atlas (16 nuclei), as distributed via
OSF and fetched with `nilearn.datasets.fetch_atlas_pauli_2017`.

| File | Contents | Space / grid |
|---|---|---|
| `CIT168_prob_MNI152NLin2009cAsym_res-1mm.nii.gz` | 4D probabilistic maps, 16 volumes | MNI152NLin2009cAsym, 1 mm (193×229×193) |
| `CIT168_det_MNI152NLin2009cAsym_res-0.7mm.nii.gz` | 3D deterministic labels | MNI152NLin2009cAsym, 0.7 mm |
| `labels.txt` | volume index → nucleus | |

Label order (prob volume index 0–15; det label value = index + 1, 0 = background):
Pu, Ca, NAC, EXA, GPe, GPi, **SNc**, RN, SNr, **PBP**, **VTA**, VeP, HN, HTH, MN, STH.

XCP-D's 4S atlases merge SNc + PBP + VTA into a single `LH/RH-SNc_PBP_VTA` parcel.

- Source: https://osf.io/r2hvk/ (OSF downloads: maps `w8zq2`, `5mqfx`; labels `6qrcb`)
- License: CC-BY 4.0 (per the OSF project; confirm there before redistribution)
- Cite: Pauli WM, Nili AN, Tyszka JM (2018). A high-resolution probabilistic in vivo
  atlas of human subcortical brain nuclei. *Scientific Data* 5:180063.
  https://doi.org/10.1038/sdata.2018.63

## atlases/Murty2014

Probabilistic VTA and SN maps (hand-drawn ROIs averaged across 50 participants; values step by 0.02).

| File | NeuroVault image |
|---|---|
| `Murty2014_VTA_bilateral_prob_MNI152_res-2mm.nii.gz` | 46838 |
| `Murty2014_VTA_L_prob_MNI152_res-2mm.nii.gz` | 46836 |
| `Murty2014_VTA_R_prob_MNI152_res-2mm.nii.gz` | 46837 |
| `Murty2014_SN_bilateral_prob_MNI152_res-2mm.nii.gz` | 46835 |

Space: FSL MNI152 2 mm (91×109×91, stored with a negative x step). This is the
same grid as MNI152NLin6Asym 2 mm, which is the CIFTI 91k subcortical voxel grid.
`*_neurovault-*.json` and `neurovault_collection_2485.json` are the original
NeuroVault metadata.

- Source: https://neurovault.org/collections/2485/
- License: CC0 (NeuroVault public collections)
- Cite: Murty VP, Shermohammed M, Smith DV, Carter RM, Huettel SA, Adcock RA (2014).
  Resting state networks distinguish human ventral tegmental area from substantia
  nigra. *NeuroImage* 100:580–589. https://doi.org/10.1016/j.neuroimage.2014.06.047

## figures/ and vta_atlas_metrics.csv

Output of the VTA atlas comparison notebook (`OCS_VTA_atlas_compare.ipynb`):

- `figures/vta_atlas_side_by_side.png`: CIT168 VTA prob | CIT168 SNc+PBP+VTA composite | Murty VTA prob
- `figures/vta_atlas_overlay.png`: outlines overlaid, 1 mm (2009cAsym) and 2 mm CIFTI grid (NLin6Asym)
- `vta_atlas_metrics.csv`: per definition × hemisphere × grid: voxel count, volume,
  prob-weighted centroid (MNI), Dice vs Murty VTA (p≥.5), fraction inside Murty SN (p≥.5)

The composite here is rebuilt from the CIT168 deterministic labels, not read from
XCP-D's own 4S dlabel, so its edges may differ by about one voxel.
