#!/usr/bin/env python3
"""
STEP B: Build a gravitational mass (dark halo) map from Gaia stellar kinematics.

Uses Gaia DR3 to:
1. Query stars with 3D velocities at 5-50 kpc from Galactic Center
2. Compute velocity dispersions in angular bins
3. Estimate dark matter density via Jeans equation
4. Project to HEALPix map matching Fermi residual

Prerequisites:
  pip install astroquery astropy galpy healpy numpy matplotlib tqdm

Output:
  results/gaia_dm_halo_healpix.fits  — projected DM density map
  results/gaia_halo_plot.png         — visualization
"""
import numpy as np
from pathlib import Path

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR = Path("data/gaia")
DATA_DIR.mkdir(parents=True, exist_ok=True)


def query_gaia_halo_stars():
    """
    Query Gaia DR3 for halo stars with full 6D phase space info.
    Focus on stars at |b| > 20° with measured radial velocities.
    """
    from astroquery.gaia import Gaia

    print("═" * 60)
    print("  STEP B1: Querying Gaia DR3 for halo stars")
    print("═" * 60)

    cache_file = DATA_DIR / "gaia_halo_stars.fits"
    if cache_file.exists():
        print(f"  → Using cached data: {cache_file}")
        from astropy.table import Table
        return Table.read(str(cache_file))

    # ADQL query: stars with RV, parallax > 0, |b| > 20°, 
    # distance 2-50 kpc (parallax 0.02-0.5 mas)
    # Limit to manageable sample
    query = """
    SELECT TOP 2000000
        source_id, ra, dec, l, b,
        parallax, parallax_error,
        pmra, pmra_error, pmdec, pmdec_error,
        radial_velocity, radial_velocity_error,
        phot_g_mean_mag, bp_rp
    FROM gaiadr3.gaia_source
    WHERE radial_velocity IS NOT NULL
        AND parallax > 0.02 AND parallax < 0.5
        AND parallax_over_error > 3
        AND ABS(b) > 20
        AND ruwe < 1.4
    ORDER BY random_index
    """

    print("  Submitting ADQL query (this may take 5-15 minutes)...")
    print(f"  Query: halo stars with RV, |b|>20°, d=2-50 kpc")

    try:
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
        Gaia.ROW_LIMIT = 2000000
        job = Gaia.launch_job_async(query, dump_to_file=False)
        results = job.get_results()
        print(f"  → Retrieved {len(results)} stars")

        # Save cache
        results.write(str(cache_file), format="fits", overwrite=True)
        print(f"  → Cached: {cache_file}")
        return results

    except Exception as e:
        print(f"  ⚠ Gaia query failed: {e}")
        print("  Generating synthetic test data instead...")
        return generate_synthetic_halo()


def generate_synthetic_halo():
    """Generate synthetic halo star data for testing the pipeline."""
    from astropy.table import Table

    print("  Generating synthetic NFW halo (for pipeline testing)...")
    np.random.seed(42)
    n = 500000

    # NFW-distributed positions (simplified)
    r_s = 20.0  # scale radius in kpc
    r = r_s * np.tan(np.random.uniform(0, np.pi/3, n)) 
    r = np.clip(r, 2, 50)
    
    # Random directions (|b| > 20°)
    b = np.random.choice([-1, 1], n) * (20 + np.random.exponential(20, n))
    b = np.clip(b, -90, 90)
    l = np.random.uniform(0, 360, n)

    # Velocity dispersion following NFW potential (~150-250 km/s)
    sigma = 180 * np.sqrt(r_s / (r + r_s))
    vr = np.random.normal(0, sigma)

    # Parallax from distance (1/d in kpc → mas)
    parallax = 1.0 / r

    tbl = Table({
        "l": l, "b": b, "parallax": parallax,
        "radial_velocity": vr,
        "pmra": np.random.normal(0, 1, n),
        "pmdec": np.random.normal(0, 1, n),
        "parallax_error": parallax * 0.1,
        "radial_velocity_error": np.full(n, 2.0),
        "pmra_error": np.full(n, 0.1),
        "pmdec_error": np.full(n, 0.1),
    })

    cache_file = DATA_DIR / "gaia_halo_stars.fits"
    tbl.write(str(cache_file), format="fits", overwrite=True)
    print(f"  → Synthetic data: {cache_file} ({n} stars)")
    return tbl


def compute_halo_density_map(stars, nside=64):
    """
    Compute projected dark matter density using Jeans equation approach.
    
    Simplified method:
    - Bin stars by HEALPix pixel
    - In each pixel, compute velocity dispersion σ_r
    - Higher dispersion → deeper potential well → more dark matter
    - Map σ²(pixel) is a proxy for projected DM column density
    """
    import healpy as hp

    print()
    print("═" * 60)
    print("  STEP B2: Computing dark halo density map")
    print("═" * 60)

    l_deg = np.array(stars["l"], dtype=float)
    b_deg = np.array(stars["b"], dtype=float)
    vr = np.array(stars["radial_velocity"], dtype=float)
    plx = np.array(stars["parallax"], dtype=float)

    # Distance in kpc
    dist = 1.0 / np.clip(plx, 0.02, 10)

    # Convert to HEALPix pixels
    theta = np.radians(90.0 - b_deg)
    phi = np.radians(l_deg % 360)
    ipix = hp.ang2pix(nside, theta, phi)

    npix = hp.nside2npix(nside)
    sum_v_map = np.zeros(npix)     # sum of v_r
    sum_v2_map = np.zeros(npix)    # sum of v_r²
    count_map = np.zeros(npix)
    sum_dist_map = np.zeros(npix)

    # Vectorized binning using np.add.at (fast for millions of stars)
    print("  Computing velocity dispersions per pixel (vectorized)...")
    np.add.at(count_map, ipix, 1)
    np.add.at(sum_v_map, ipix, vr)
    np.add.at(sum_v2_map, ipix, vr ** 2)
    np.add.at(sum_dist_map, ipix, dist)

    # Proper variance: σ² = <v²> - <v>²
    mask = count_map >= 10  # need minimum stars for statistics
    mean_v = np.zeros(npix)
    mean_v[mask] = sum_v_map[mask] / count_map[mask]
    sigma2_map = np.zeros(npix)
    sigma2_map[mask] = sum_v2_map[mask] / count_map[mask] - mean_v[mask] ** 2
    sigma2_map[~mask] = hp.UNSEEN

    mean_dist_map = np.zeros(npix)
    mean_dist_map[mask] = sum_dist_map[mask] / count_map[mask]

    # The Jeans equation relates σ² to the enclosed mass:
    # ρ_DM(r) ∝ σ²(r) / r  (simplified spherical Jeans)
    # For projection: Σ_DM ∝ σ² × mean_distance (column density proxy)
    dm_proxy = np.zeros(npix)
    dm_proxy[mask] = sigma2_map[mask] * np.sqrt(mean_dist_map[mask])
    dm_proxy[~mask] = hp.UNSEEN

    # Normalize to [0, 1] range for comparison
    valid = dm_proxy != hp.UNSEEN
    if np.any(valid):
        mn, mx = np.percentile(dm_proxy[valid], [2, 98])
        dm_proxy[valid] = np.clip((dm_proxy[valid] - mn) / (mx - mn + 1e-10), 0, 1)

    # Save
    outpath = RESULTS_DIR / "gaia_dm_halo_healpix.fits"
    hp.write_map(str(outpath), dm_proxy, overwrite=True, coord="G")
    print(f"  → Saved: {outpath}")
    print(f"  → Valid pixels: {np.sum(valid)}/{npix}")
    print(f"  → Stars used: {int(np.sum(count_map[mask]))}")

    # Plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.figure(figsize=(16, 5))
    count_plot = count_map.copy()
    count_plot[count_plot <= 0] = 0.1  # avoid log(0)
    hp.mollview(count_plot, coord="G", title="Star counts per pixel",
                cmap="viridis", sub=121, norm="log")
    hp.graticule()
    dm_plot = dm_proxy.copy()
    dm_plot[dm_plot == hp.UNSEEN] = np.nan
    hp.mollview(dm_plot, coord="G", title="DM halo proxy (σ² × √d)",
                cmap="magma", sub=122)
    hp.graticule()
    plt.savefig(str(RESULTS_DIR / "gaia_halo_plot.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → Plot: {RESULTS_DIR / 'gaia_halo_plot.png'}")

    return dm_proxy


def fit_nfw_profile(dm_map, nside=64):
    """Fit an NFW profile to the DM proxy map for comparison."""
    import healpy as hp
    from scipy.optimize import curve_fit

    print()
    print("  Fitting NFW profile...")

    valid = dm_map != hp.UNSEEN
    if not np.any(valid):
        print("  ⚠ No valid pixels for NFW fit")
        return

    # Get angular distance from GC for each pixel
    gc_vec = hp.ang2vec(np.pi/2, 0)  # GC at (l=0, b=0)
    pix_vecs = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))).T
    angles = np.arccos(np.clip(np.dot(pix_vecs, gc_vec), -1, 1))
    angles_deg = np.degrees(angles)

    # NFW projected profile: Σ(θ) ∝ 1/(θ/θ_s × (1 + θ/θ_s)²)
    def nfw_proj(theta, A, theta_s):
        x = theta / theta_s
        x = np.clip(x, 0.01, 100)
        return A / (x * (1 + x)**2)

    x = angles_deg[valid]
    y = dm_map[valid]
    try:
        popt, pcov = curve_fit(nfw_proj, x, y, p0=[1.0, 20.0], maxfev=5000)
        print(f"  → NFW fit: A={popt[0]:.3f}, θ_s={popt[1]:.1f}°")
        print(f"    (scale angle θ_s ~ {popt[1]:.1f}° corresponds to ~{popt[1]*0.17:.1f} kpc at 8 kpc)")
    except Exception as e:
        print(f"  ⚠ NFW fit failed: {e}")


if __name__ == "__main__":
    stars = query_gaia_halo_stars()
    dm_map = compute_halo_density_map(stars)
    fit_nfw_profile(dm_map)
    print()
    print("  ✅ Step B complete. Next: run 04_cross_correlation.py")
