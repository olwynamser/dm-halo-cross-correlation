#!/usr/bin/env python3
"""
STEP A2: Reproduce Totani's Galactic Center gamma-ray residual map.

Uses fermipy to:
1. Select and bin Fermi LAT photon data (10-50 GeV)
2. Fit the standard diffuse model + 4FGL point sources
3. Subtract the model → residual map
4. Save the residual as a HEALPix map

Prerequisites:
  pip install fermipy astropy healpy numpy matplotlib
  Completed step 01 (data files in data/fermi/)

Output:
  results/fermi_residual_healpix.fits  — residual gamma-ray map
  results/fermi_residual_plot.png      — visualization
"""
import os
import yaml
import numpy as np
from pathlib import Path

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def create_fermipy_config():
    """Create fermipy configuration YAML for the Galactic Center analysis."""

    config = {
        "data": {
            "evfile": "data/fermi/photons.fits",
            "scfile": "data/fermi/spacecraft.fits",
            "ltcube": None,  # will be generated
        },
        "binning": {
            "roiwidth": 40.0,    # 40° × 40° ROI
            "binsz": 0.2,        # 0.2° pixel size
            "projtype": "WCS",
        },
        "selection": {
            "emin": 10000,       # 10 GeV in MeV
            "emax": 50000,       # 50 GeV in MeV
            "glat": 0.0,
            "glon": 0.0,
            "coordsys": "GAL",
            "zmax": 90,          # zenith cut
            "evclass": 128,      # SOURCE class
            "evtype": 3,         # FRONT+BACK
            "tmin": 239557417,   # mission start
            "tmax": 725846400,
        },
        "model": {
            "src_roiwidth": 45.0,
            "galdiff": "data/fermi/gll_iem_v07.fits",
            "isodiff": "data/fermi/iso_P8R3_SOURCE_V3_v1.txt",
            "catalogs": ["data/fermi/gll_psc_v35.fit"],
        },
        "ltcube": {
            "zmax": 90,
        },
        "fileio": {
            "outdir": "results/fermipy_gc",
            "logfile": "results/fermipy_gc/fermipy.log",
        },
        "optimizer": {
            "optimizer": "MINUIT",
        },
    }

    config_path = "config_gc.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False)
    print(f"  Created config: {config_path}")
    return config_path


def run_analysis(config_path):
    """Run the full fermipy analysis pipeline."""
    from fermipy.gtanalysis import GTAnalysis

    print("═" * 60)
    print("  STEP A2: Fermipy Analysis — Galactic Center")
    print("═" * 60)
    print("  This will take 30-120 minutes depending on your machine.")
    print()

    # Initialize
    print("  [1/5] Initializing GTAnalysis...")
    gta = GTAnalysis(config_path, logging={"verbosity": 3})

    # Setup (generates livetime cube, exposure map, source maps)
    print("  [2/5] Running setup (ltcube, expmap, srcmaps)...")
    gta.setup()

    # Free parameters of nearby bright sources and diffuse models
    print("  [3/5] Freeing source parameters...")
    gta.free_sources(minmax_ts=[100, None], pars="norm")
    gta.free_source("galdiff", pars=["Prefactor", "Index"], free=True)
    gta.free_source("isodiff", pars=["Normalization"], free=True)

    # Fit
    print("  [4/5] Running likelihood fit...")
    fit_results = gta.fit()
    print(f"       Fit quality: {fit_results['fit_quality']}")
    print(f"       Log-likelihood: {fit_results['loglike']:.1f}")

    # Generate residual map
    print("  [5/5] Generating residual map...")
    resid = gta.residmap("gc_residual")

    # Save results
    gta.write_roi("gc_final")
    print()
    print(f"  ✅ Analysis complete. Results in: results/fermipy_gc/")
    return gta, resid


def extract_residual_healpix(gta):
    """Convert the WCS residual map to HEALPix format."""
    import healpy as hp
    from astropy.io import fits
    from astropy.wcs import WCS

    print("  Converting residual to HEALPix...")

    # Load the residual counts map
    resid_file = "results/fermipy_gc/gc_residual_residmap.fits"
    if not os.path.exists(resid_file):
        # Try alternate naming
        resid_file = "results/fermipy_gc/gc_final_residmap.fits"

    with fits.open(resid_file) as hdul:
        data = hdul[0].data
        wcs = WCS(hdul[0].header)

    # Create HEALPix map
    nside = 64  # ~0.9° resolution
    npix = hp.nside2npix(nside)
    hpx_map = np.zeros(npix)
    counts_map = np.zeros(npix)

    # Project WCS pixels to HEALPix
    ny, nx = data.shape[-2:]
    for iy in range(ny):
        for ix in range(nx):
            val = data[iy, ix] if data.ndim == 2 else data[0, iy, ix]
            if np.isnan(val) or val == 0:
                continue
            coords = wcs.pixel_to_world(ix, iy)
            l_deg = coords.galactic.l.deg
            b_deg = coords.galactic.b.deg
            theta = np.radians(90.0 - b_deg)
            phi = np.radians(l_deg)
            ipix = hp.ang2pix(nside, theta, phi)
            hpx_map[ipix] += val
            counts_map[ipix] += 1

    # Average where multiple WCS pixels map to same HEALPix pixel
    mask = counts_map > 0
    hpx_map[mask] /= counts_map[mask]

    # Save
    outpath = RESULTS_DIR / "fermi_residual_healpix.fits"
    hp.write_map(str(outpath), hpx_map, overwrite=True, coord="G")
    print(f"  → Saved: {outpath}")

    # Plot
    plot_path = RESULTS_DIR / "fermi_residual_plot.png"
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    hp.mollview(hpx_map, coord="G", title="Fermi LAT Residual (10-50 GeV)",
                unit="counts", cmap="inferno", min=-2, max=5)
    hp.graticule()
    plt.savefig(str(plot_path), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → Plot saved: {plot_path}")

    return hpx_map


def run_simplified():
    """
    Simplified version: if fermipy is not available or too complex,
    perform a basic photon count analysis.
    """
    import healpy as hp
    from astropy.io import fits
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    print("═" * 60)
    print("  STEP A2 (SIMPLIFIED): Direct photon binning")
    print("═" * 60)
    print("  Note: This is a simplified version without full likelihood.")
    print("  It bins photons directly into HEALPix and shows the raw map.")
    print()

    photon_file = "data/fermi/photons.fits"
    if not os.path.exists(photon_file):
        print(f"  ❌ File not found: {photon_file}")
        print("     Run 01_download_fermi.py first.")
        return None

    print("  Loading photon data...")
    with fits.open(photon_file) as hdul:
        events = hdul[1].data
        ra = events["RA"]
        dec = events["DEC"]
        energy = events["ENERGY"]  # MeV

    # Filter 10-50 GeV
    mask = (energy >= 10000) & (energy <= 50000)
    ra_sel = ra[mask]
    dec_sel = dec[mask]
    print(f"  Selected {np.sum(mask)} photons in 10-50 GeV range")

    # Convert to Galactic coordinates
    coords = SkyCoord(ra=ra_sel * u.deg, dec=dec_sel * u.deg, frame="icrs")
    l_gal = coords.galactic.l.deg
    b_gal = coords.galactic.b.deg

    # Bin into HEALPix
    nside = 64
    npix = hp.nside2npix(nside)
    theta = np.radians(90.0 - b_gal)
    phi = np.radians(l_gal)
    ipix = hp.ang2pix(nside, theta, phi)

    hpx_map = np.zeros(npix)
    for p in ipix:
        hpx_map[p] += 1

    # Simple background: smooth the map heavily and subtract
    # (this is a crude approximation — full analysis uses likelihood)
    smoothed = hp.smoothing(hpx_map, sigma=np.radians(5.0))
    residual = hpx_map - smoothed

    # Save
    outpath = RESULTS_DIR / "fermi_residual_healpix.fits"
    hp.write_map(str(outpath), residual, overwrite=True, coord="G")
    print(f"  → Saved: {outpath}")

    # Plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.figure(figsize=(16, 5))
    hp.mollview(hpx_map, coord="G", title="Raw photon counts (10-50 GeV)",
                cmap="inferno", sub=121)
    hp.graticule()
    hp.mollview(residual, coord="G", title="Residual (data - smoothed background)",
                cmap="RdBu_r", sub=122, min=-3, max=3)
    hp.graticule()
    plt.savefig(str(RESULTS_DIR / "fermi_residual_plot.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → Plot saved: {RESULTS_DIR / 'fermi_residual_plot.png'}")

    return residual


if __name__ == "__main__":
    # Try full fermipy analysis first; fall back to simplified
    try:
        import fermipy
        config_path = create_fermipy_config()
        gta, resid = run_analysis(config_path)
        hpx = extract_residual_healpix(gta)
    except ImportError:
        print("  ⚠ fermipy not installed. Running simplified analysis.")
        print("    (For full analysis: conda install -c conda-forge fermipy)")
        print()
        hpx = run_simplified()
    except Exception as e:
        print(f"  ⚠ fermipy failed ({e}). Running simplified analysis.")
        hpx = run_simplified()

    if hpx is not None:
        print()
        print("  ✅ Step A2 complete. Next: run 03_gaia_halo.py")
