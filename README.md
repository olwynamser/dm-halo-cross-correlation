Dark Matter Cross-Correlation Pipeline
Goal
Cross-correlate Fermi LAT gamma-ray data with Gaia DR3 stellar halo kinematics to test for spatial coincidence between the gamma-ray halo component and the stellar mass distribution of the Milky Way halo.

If the gamma-ray halo component correlates spatially with the stellar halo → the signal is tied to the halo mass distribution. This is consistent with dark matter annihilation, but baryonic sources (millisecond pulsars, inverse Compton) cannot yet be excluded.


Results Summary
v1.0 — Full-sky GALPROP residual (10–50 GeV)
Metric
Value
Spearman ρ
+0.080
Significance
3.3σ
Bootstrap 95% CI
[0.041, 0.122]
Photons
687,491
Gaia stars
500,000


Validation: Sculptor dwarf spheroidal σ_r = 9.3 km/s (literature: 9.2 km/s).
v2.0 — Totani 21 GeV halo template
Using a professionally decomposed halo component map (Totani, priv. comm.; southern sky patch b ∈ [−60°, −28°]):

Metric
Value
Spearman ρ
+0.194
Significance
11.6σ
Bootstrap 95% CI
[0.161, 0.227]
Valid pixels
3,394


The increase relative to v1.0 reflects both the improved template and the different sky region.
Selection stringency test
Selection
Criteria
ρ
σ
Loose
pm > 3, plx < 1.0
+0.174
10.4
Standard
pm > 5, plx < 0.5
+0.194
11.3
Strict
pm > 8, plx < 0.3
+0.226
13.9
Very strict
pm > 12, plx < 0.2
+0.306
19.1
Disk control
pm < 3, plx > 2.0
+0.017
0.0

Null tests
Test
Result
Rotation 180°/270°
ρ → 0 (signal is spatial)
Pixel shuffle (1000 iter)
11.2σ
Latitude detrending
ρ = 0.141, 8.2σ (persists)

Energy dependence (GALPROP)
Band
ρ
Note
1–10 GeV
−0.075
Sign inversion (background-dominated)
10–50 GeV
+0.080
Positive (3.3σ)



Limitations
Totani halo map covers only b ∈ [−60°, −28°] — partial sky
Single energy band (21 GeV) — cannot distinguish DM from baryonic halo sources (MSPs, inverse Compton)
Correlation may trace gravitational potential without requiring DM
Bootstrap CIs assume pixel independence — spatial autocorrelation not corrected
Gaia queries use random 500K subsamples, not full catalog


Quick Start
# 1. Install dependencies

pip install -r requirements.txt

# 2. Download Fermi data (follow manual instructions in output)

python 01_download_fermi.py

# 3. Build gamma-ray residual map

python 02_fermi_residual.py

# 4. Build dark matter halo map from Gaia

python 03_gaia_halo.py

# 5. Run cross-correlation

python 04_cross_correlation.py


Pipeline
Path A (v1.0):

  Fermi LAT (10-50 GeV) → HEALPix → GALPROP subtraction → Residual map ─┐

                                                                          ├→ Spearman ρ

  Gaia DR3 (halo stars) → HEALPix → Density map ────────────────────────┘

Path B (v2.0):

  Totani 21 GeV halo → HEALPix rebinning → Template map ────────────────┐

                                                                          ├→ Spearman ρ

  Gaia DR3 (halo stars) → HEALPix → Density map ────────────────────────┘


Files
File
Description
01_download_fermi.py
Fermi LAT data query and download
02_fermi_residual.py
GALPROP model subtraction
03_gaia_halo.py
Gaia DR3 halo star selection and mapping
04_cross_correlation.py
Spearman correlation with bootstrap
professional_pipeline.py
Full automated pipeline
Untitled0 (1).ipynb
Analysis notebook with all tests



Citation
Amser, O. (2026). Cross-correlation of Gaia DR3 Stellar Kinematics with Fermi LAT 

Gamma-ray Residual and Halo Template: Evidence for Spatial Coincidence up to 19.1σ. 

Zenodo. https://doi.org/10.5281/zenodo.19221429


Acknowledgements
Prof. Tomonori Totani (University of Tokyo) — 21 GeV halo component map
Prof. Oscar Macias (San Francisco State University) — facilitating contact with C. Eckner
ESA Gaia DR3, NASA Fermi LAT
