# Gaia DR3 x Totani (2025) -- Kinematic Halo Cross-Correlation

[![Zenodo](https://img.shields.io/badge/Zenodo-v4.0-blue)](https://zenodo.org/records/19268430)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19268430-green)](https://doi.org/10.5281/zenodo.19268430)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Kinematic Cross-Correlation of Gaia DR3 Stellar Halo with Totani (2025) 21 GeV Gamma-Ray Halo: Spectral Analysis and Multi-Test Validation**

*Olwyn Amser -- Independent Researcher -- March 2026*

---

## Overview

This repository contains the analysis code for a spatial cross-correlation between:
- **Gaia DR3** kinematic stellar halo maps (proper motion selection, NSIDE=64 HEALPix)
- **Totani (2025, JCAP, arXiv:2507.07209)** energy-binned gamma-ray halo residual maps at 1.5-170 GeV

**Main result:** partial(|b|) = **+0.154 +/- 0.016**, **9.4sigma raw**, **9.1sigma Bonferroni** (22 tests)
Ring region: 20<|l|<60 deg, 10<|b|<60 deg.

The signal is **halo-specific**: disk stars give 2.2sigma, very strict halo selection (pm>12 mas/yr) gives 9.1sigma.

---

## Key Results (v4.0)

### Spatial Correlation

| Selection | partial(|b|) | sigma | 95% CI |
|-----------|-------------|-------|--------|
| Disk (null) | +0.036 | 2.2s | [+0.005, +0.068] |
| Standard (pm>5) | +0.066 | 3.9s | [+0.037, +0.101] |
| Strict (pm>8) | +0.071 | 4.6s | [+0.043, +0.104] |
| **Very strict (pm>12)** | **+0.154** | **9.4s (9.1s+)** | **[+0.123, +0.184]** |

+ After Bonferroni correction for 22 independent tests. Disk/halo ratio = 0.23.

### Spectral Analysis

| Energy (GeV) | partial(|b|) | SNR | Status |
|-------------|-------------|-----|--------|
| 1.5 | -0.058 | -2.3s | EXCLUDED: GALPROP oversubtraction |
| 2.5 | +0.008 | 0.4s | EXCLUDED: SUM->AVG delta=238% |
| 4.3 | +0.003 | 0.1s | EXCLUDED: SUM->AVG unstable |
| **12** | **+0.135** | **7.4s** | RELIABLE |
| **21** | **+0.154** | **9.4s** | PRIMARY |
| **35** | **+0.048** | **2.5s** | RELIABLE |
| **59** | **+0.052** | **2.2s** | BORDERLINE |
| 100 | +0.073 | 4.0s | ANOMALOUS (see note v4.0 s4.2) |
| 170 | +0.006 | 0.3s | NOT SIGNIFICANT |

Spectral profile peaked at 12-35 GeV, consistent with bb-bar at m_chi ~ 0.6 TeV.

### Angular Power Spectrum C_ell (|b|-detrended)

After rank-regression detrending of the |b| gradient from both maps before anafast:

| ell range | Scale | Coherence | Interpretation |
|-----------|-------|-----------|----------------|
| ell=2-3 | ~60-90 deg | +0.21 to +0.50 | Large-scale halo |
| **ell=4-9** | **~20-45 deg** | **+0.52 mean** | **Ring scale** |
| ell=10-32 | ~6-18 deg | +0.16 to +0.60 | Sub-ring structure |

Peak coherence: ell=7 (~26 deg). Without detrending: peak ell=2 (gradient-dominated, ring_coh=+0.21).
After detrending: coherence 2.5x higher at ring scales.

### Robustness Tests

| Test | Result | Pass |
|------|--------|------|
| Scramble (N=300) | p=0.0000 | YES |
| SFD E(B-V) dust | 9.7sigma after PCA | YES |
| Sagittarius excl. | 7.7sigma survives | YES |
| Hemisphere N vs S | 6.4s vs 6.1s | YES |
| N=10 subsamples std | 0.016 | YES |
| Fragility clip 99% | delta=2.7% | YES |
| NSIDE 32/64/128 | 6.9/9.3/8.7sigma | YES |
| ICS CMB (21 GeV) | 97% survival | YES |
| Bulge null | CI includes 0 | YES |
| MI raw (k=5,7,10) | 14.5sigma mean | YES |

### Source Attribution

| Hypothesis | P |
|-----------|---|
| Correlation real | 99% |
| Signal halo-specific | 95% |
| Not dust artifact | 99% |
| Not spiral arm artifact | 90% |
| Not pure ICS (21 GeV) | 75% |
| **DM bb-bar (m_chi~0.6 TeV)** | **40%** |
| Unknown halo source | 30% |
| GALPROP artifact | 10% |
| Random artifact | 1% |

Note: Expert estimates (Consilium v5, 14 agents, 4 cycles). Not formal Bayesian derivation.

---

## Method

Statistic: Spearman partial rank correlation

    partial(|b|) = [rho(Gaia,Totani) - rho(Gaia,|b|)*rho(Totani,|b|)]
                  / sqrt[(1-rho(Gaia,|b|)^2)(1-rho(Totani,|b|)^2)]

Uncertainty: Block bootstrap (N=500, block_size=10 pixels ~9 deg, seed=42)
ROI: 20<|l|<60 deg, 10<|b|<60 deg (24% sky, N_pix=11,804 at NSIDE=64)
fort.200: Pixel-averaging (read_fort200_avg) -- Spearman-invariant verified.

---

## Repository Structure
```
dm-halo-cross-correlation/
|-- README.md
|-- analysis/
|   |-- astro_engine_v20.py
|   |-- consilium_v5.py
|   |-- doc_consilium.py
|   `-- read_fort200_avg.py
|-- data/
|   |-- analysis_arrays.npz
|   |-- totani_64.npy
|   |-- gaia_map_*.npy
|   |-- spectral_avg_final.npy
|   |-- consilium_v5_FINAL.npy
|   |-- cl_detrended.npy
|   `-- dust_test_sfd.npy
`-- research_note_v4.docx
```

---

## Version History

| Version | Date | DOI | Key change |
|---------|------|-----|------------|
| v1.0 | 2026-03-25 | 10.5281/zenodo.19221430 | Initial result |
| v2.0 | 2026-03-26 | 10.5281/zenodo.19236371 | Full ROI, hemisphere |
| v3.0 | 2026-03-27 | 10.5281/zenodo.19244340 | Spectral, ICS test |
| **v4.0** | **2026-03-28** | **10.5281/zenodo.19268430** | Bug fix, dust, Sgr, C_ell |

Changes v3.0 to v4.0:
1. read_fort200 corrected: summation -> pixel-averaging
2. MI=21.4sigma retracted (v2.0 implementation error); correct: raw 11.8-14.5sigma
3. SFD E(B-V) dust control added: 9.7sigma after PCA
4. Sagittarius stream exclusion: 7.7sigma survives
5. Five post-hoc exclusions documented (research note s4.3)
6. C_ell detrended analysis added: peak ell=7, ring_coh=+0.52
7. Consilium v5: GREEN=30, YELLOW=21, RED=0
8. Document audit: GREEN=54, YELLOW=4, RED=0

---

## Open Questions (Pending External Data)

| Item | Status | Contact |
|------|--------|---------|
| ICS(E) at 12/35/59 GeV | Requested | T. Totani |
| skyFACT MSP template | Requested | C. Eckner |
| LAMOST DR9 metallicity | In progress | -- |
| BHB star cross-match | Not started | -- |
| Alt. GALPROP models | Requested | F. Calore |

Decisive test: ICS(E) control + MSP subtraction + metallicity cut.
If signal survives all three -> DM strongly favored (P ~60%+).

---

## Known Limitations

1. Single GALPROP model tested
2. Totani map = model residual, not raw Fermi photons
3. ICS templates only at 21 GeV
4. No MSP spatial template (skyFACT pending)
5. No metallicity cut (LAMOST DR9 in progress)
6. No preregistration -- 5 post-hoc decisions documented
7. Self-administered audit -- no external peer review
8. Fermi exposure via |b| proxy only
9. |l|=30-60 deg overlaps Galactic bar ends
10. ICS files rho=0.95 -- effectively one template

Full list: research note s5 (21 limitations total).

---

## Quick Start (Colab)
```python
# Restore after runtime reset
from google.colab import drive
drive.mount('/content/drive')
import numpy as np, healpy as hp
from scipy.stats import spearmanr

data      = np.load('/content/drive/MyDrive/dm_halo_project/data/analysis_arrays.npz')
idx       = data['idx']; b_vals = data['b_vals']; l_vals = data['l_vals']
totani_64 = np.load('/content/drive/MyDrive/dm_halo_project/data/totani_64.npy')
gaia_vs   = np.load('/content/drive/MyDrive/dm_halo_project/data/gaia_map_very_strict.npy')
g_vs      = gaia_vs[idx]

NSIDE = 64; NPIX = hp.nside2npix(NSIDE)
abs_b = np.abs(b_vals); abs_l = np.abs(l_vals)
ring  = (abs_l>=20)&(abs_l<=60)&(abs_b>=10)&(abs_b<=60)
ok_r  = (totani_64[idx]>0)&ring&(g_vs>0)

def partial_1var(x, y, z):
    rxy,_ = spearmanr(x,y); rxz,_ = spearmanr(x,z); ryz,_ = spearmanr(y,z)
    d = np.sqrt((1-rxz**2)*(1-ryz**2))
    return (rxy-rxz*ryz)/d if d>0 else 0

result = partial_1var(g_vs[ok_r], totani_64[idx][ok_r], abs_b[ok_r])
print(f"partial(|b|) = {result:+.4f}")  # Expected: ~+0.154
```

---

## Citation
```bibtex
@misc{amser2026gaia,
  author    = {Amser, Olwyn},
  title     = {{Kinematic Cross-Correlation of Gaia DR3 Stellar Halo
                with Totani (2025) 21 GeV Gamma-Ray Halo}},
  year      = {2026},
  publisher = {Zenodo},
  version   = {v4.0},
  doi       = {10.5281/zenodo.19268430},
  url       = {https://zenodo.org/records/19268430}
}
```

---

## References

- Totani, T. (2025). JCAP. arXiv:2507.07209
- Gaia Collaboration et al. (2023). Gaia DR3. A&A
- Schlegel, D., Finkbeiner, D., Davis, M. (1998). ApJ, 500, 525
- Necib, L. et al. (2019). ApJ, 883, 27
- Bonaca, A. et al. (2021). ApJ, 909, L26
- Gorski, K.M. et al. (2005). ApJ, 622, 759

---

Last updated: March 28, 2026 | v4.0
