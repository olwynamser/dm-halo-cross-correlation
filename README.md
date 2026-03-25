# Dark Matter Cross-Correlation Pipeline

## Goal

Verify the Totani (2025) dark matter signal by cross-correlating the Fermi LAT
gamma-ray residual at the Galactic Center with an independent dark matter halo
map derived from Gaia DR3 stellar kinematics.

**If the gamma-ray excess has the same spatial shape as the gravitational dark
matter halo → two independent methods see the same thing → strong evidence for
dark matter annihilation.**

---

## Quick Start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Download Fermi data (follow manual instructions in output)
python scripts/01_download_fermi.py

# 3. Build gamma-ray residual map
python scripts/02_fermi_residual.py

# 4. Build dark matter halo map from Gaia
python scripts/03_gaia_halo.py

# 5. Cross-correlate and measure significance
python scripts/04_cross_correlation.py
```

All results will be in the `results/` folder.

---

## What You Need

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| RAM | 8 GB | 32 GB |
| Disk | 20 GB | 100 GB |
| Internet | Yes (for data download) | Fast |
| Python | 3.9+ | 3.11 |
| GPU | Not needed | Not needed |
| Time | 1-2 days | 3-5 days (with full fermipy) |

Works on: Linux, Mac, Windows (WSL2).
Google Colab Pro ($10/mo) is sufficient.

---

## Pipeline Overview

### Step A: Fermi Gamma-Ray Residual (Scripts 01, 02)

1. Download Fermi LAT photon data (10-50 GeV, 40° around GC)
2. Subtract known sources (4FGL catalog) and diffuse Galactic model
3. What remains = potential dark matter signal (the "Galactic Center Excess")
4. Output: HEALPix map of residual gamma-ray emission

The script has two modes:
- **Full mode**: Uses `fermipy` for proper likelihood analysis (recommended)
- **Simplified mode**: Direct photon binning with smoothing subtraction (faster)

### Step B: Gaia Dark Halo Map (Script 03)

1. Query Gaia DR3 for halo stars with 3D velocities (5-50 kpc from GC)
2. Compute radial velocity dispersion σ_r in each angular bin
3. σ² is a proxy for the depth of the gravitational potential → dark matter density
4. Output: HEALPix map of projected dark matter density

If Gaia query fails, the script generates synthetic NFW halo data for pipeline testing.

### Step C: Cross-Correlation (Script 04)

1. **Pixel correlation**: Pearson r and Spearman ρ between maps
2. **Angular cross-spectrum**: C_ℓ decomposition → which angular scales correlate?
3. **Shuffle test**: Randomize one map 10,000× to estimate significance in σ
4. **Radial profile comparison**: Do both maps fall off the same way from GC?

Output: significance in σ, p-value, and diagnostic plots.

---

## Interpreting Results

| Significance | Meaning | Action |
|-------------|---------|--------|
| < 2σ | No correlation | Signal is likely not dark matter |
| 2-3σ | Hint | Promising, needs more data / better method |
| 3-5σ | Evidence | Strong indication; write paper |
| > 5σ | Discovery | Contact arXiv immediately |

---

## Key References

- Totani (2025), "20 GeV halo-like excess", JCAP 2025(11):080
- Ackermann et al. (2017), "The Fermi Galactic Center GeV Excess"
- McMillan (2017), "The mass distribution of the Milky Way"
- Fermipy: Wood et al. (2017), PoS ICRC2017, 824

---

## Output Files

| File | Script | Description |
|------|--------|-------------|
| `results/fermi_residual_healpix.fits` | 02 | Gamma-ray residual HEALPix map |
| `results/fermi_residual_plot.png` | 02 | Visualization of residual |
| `results/gaia_dm_halo_healpix.fits` | 03 | DM halo density HEALPix map |
| `results/gaia_halo_plot.png` | 03 | Visualization of halo |
| `results/cross_correlation_results.txt` | 04 | Significance summary |
| `results/cross_correlation_plot.png` | 04 | Diagnostic plots |

---

## File Structure

```
dark_matter_project/
├── README.md
├── requirements.txt
├── scripts/
│   ├── 01_download_fermi.py    — Data acquisition
│   ├── 02_fermi_residual.py    — Gamma-ray analysis
│   ├── 03_gaia_halo.py         — Gaia halo reconstruction
│   └── 04_cross_correlation.py — Statistical analysis
├── data/                        — Downloaded data (gitignore)
│   ├── fermi/
│   └── gaia/
└── results/                     — Output maps and plots
```

---

## Limitations & Caveats

- The simplified Fermi analysis (no fermipy) is a rough approximation.
  For publication-quality results, use the full fermipy pipeline.
- The Jeans equation approach for the DM halo is simplified. A full
  analysis would use galpy or AGAMA for orbit-based modeling.
- Systematic uncertainties in the Galactic diffuse model are the
  dominant source of error. The cross-correlation with Gaia is
  specifically designed to be independent of this model.
- The Gaia stellar sample at large distances (>15 kpc) is incomplete.
  This limits the halo map resolution in outer regions.

---

## License

This pipeline is released into the public domain. Use it for anything.
If you publish results, cite Totani (2025) and this pipeline.
