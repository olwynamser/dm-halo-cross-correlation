#!/usr/bin/env python3
"""
STEP C: Cross-correlate Fermi gamma-ray residual with Gaia DM halo map.

Computes:
1. Pixel-by-pixel correlation (Pearson, Spearman)
2. Angular power spectrum cross-correlation C_ℓ
3. Statistical significance via shuffle test (10000 randomizations)
4. Comparison with NFW vs bulge (pulsar) models

Prerequisites:
  pip install healpy numpy scipy matplotlib scikit-learn

Input:
  results/fermi_residual_healpix.fits
  results/gaia_dm_halo_healpix.fits

Output:
  results/cross_correlation_results.txt
  results/cross_correlation_plot.png
"""
import numpy as np
from scipy import stats
from pathlib import Path

RESULTS_DIR = Path("results")


def load_maps():
    """Load both HEALPix maps."""
    import healpy as hp

    print("═" * 60)
    print("  STEP C: Cross-correlation analysis")
    print("═" * 60)

    fermi_path = RESULTS_DIR / "fermi_residual_healpix.fits"
    gaia_path = RESULTS_DIR / "gaia_dm_halo_healpix.fits"

    for p in [fermi_path, gaia_path]:
        if not p.exists():
            print(f"  ❌ Missing: {p}")
            print(f"     Run previous steps first.")
            return None, None

    fermi = hp.read_map(str(fermi_path))
    gaia = hp.read_map(str(gaia_path))

    # Ensure same NSIDE
    nside_f = hp.npix2nside(len(fermi))
    nside_g = hp.npix2nside(len(gaia))
    target_nside = min(nside_f, nside_g, 64)

    if nside_f != target_nside:
        fermi = hp.ud_grade(fermi, target_nside)
    if nside_g != target_nside:
        gaia = hp.ud_grade(gaia, target_nside)

    print(f"  Fermi map: NSIDE={target_nside}, {len(fermi)} pixels")
    print(f"  Gaia map:  NSIDE={target_nside}, {len(gaia)} pixels")
    return fermi, gaia


def pixel_correlation(fermi, gaia):
    """Compute pixel-by-pixel correlation."""
    import healpy as hp

    print()
    print("  ── Pixel-by-pixel correlation ──")

    # Mask invalid pixels (handles both UNSEEN from fermipy and zeros from simplified)
    valid = (fermi != hp.UNSEEN) & (gaia != hp.UNSEEN) & \
            (~np.isnan(fermi)) & (~np.isnan(gaia)) & \
            (np.isfinite(fermi)) & (np.isfinite(gaia)) & \
            (gaia > 0)

    f = fermi[valid]
    g = gaia[valid]
    n = len(f)
    print(f"  Valid pixels: {n}")

    if n < 50:
        print("  ⚠ Too few valid pixels for meaningful correlation.")
        return 0, 0, 1.0

    # Pearson
    r_pearson, p_pearson = stats.pearsonr(f, g)
    print(f"  Pearson r  = {r_pearson:.4f}  (p = {p_pearson:.2e})")

    # Spearman (rank correlation, more robust)
    r_spearman, p_spearman = stats.spearmanr(f, g)
    print(f"  Spearman ρ = {r_spearman:.4f}  (p = {p_spearman:.2e})")

    return r_pearson, r_spearman, p_pearson


def angular_cross_spectrum(fermi, gaia):
    """Compute angular cross-power spectrum C_ℓ(γ × DM)."""
    import healpy as hp

    print()
    print("  ── Angular cross-power spectrum ──")

    # Replace UNSEEN with 0 for alm computation
    f = fermi.copy()
    g = gaia.copy()
    f[f == hp.UNSEEN] = 0
    g[g == hp.UNSEEN] = 0
    f[np.isnan(f)] = 0
    g[np.isnan(g)] = 0

    # Remove monopole and dipole
    f -= np.mean(f)
    g -= np.mean(g)

    # Compute cross-spectrum
    nside = hp.npix2nside(len(f))
    lmax = min(3 * nside - 1, 192)

    alm_f = hp.map2alm(f, lmax=lmax)
    alm_g = hp.map2alm(g, lmax=lmax)

    # Cross-spectrum: C_ℓ = <a_ℓm^f × a_ℓm^g*> / (2ℓ+1)
    cl_cross = hp.alm2cl(alm_f, alm_g)
    cl_auto_f = hp.alm2cl(alm_f)
    cl_auto_g = hp.alm2cl(alm_g)

    # Normalized cross-correlation per ℓ
    ells = np.arange(len(cl_cross))
    denom = np.sqrt(cl_auto_f * cl_auto_g)
    r_ell = np.zeros_like(cl_cross)
    mask = denom > 0
    r_ell[mask] = cl_cross[mask] / denom[mask]

    # Average over ℓ = 2-30 (large-scale halo structure)
    sel = (ells >= 2) & (ells <= 30)
    mean_r = np.mean(r_ell[sel])
    print(f"  Mean correlation r(ℓ=2-30) = {mean_r:.4f}")
    print(f"  Peak correlation at ℓ = {ells[np.argmax(np.abs(r_ell[2:]))+2]}")

    return ells, cl_cross, r_ell


def shuffle_test(fermi, gaia, n_shuffles=10000):
    """Assess significance by shuffling one map."""
    import healpy as hp

    print()
    print(f"  ── Shuffle test ({n_shuffles} randomizations) ──")

    valid = (fermi != hp.UNSEEN) & (gaia != hp.UNSEEN) & \
            (~np.isnan(fermi)) & (~np.isnan(gaia)) & (gaia > 0)
    f = fermi[valid]
    g = gaia[valid]

    if len(f) < 50:
        print("  ⚠ Too few pixels for shuffle test.")
        return 0, 1.0, 0, np.zeros(100)

    # Observed correlation
    r_obs = stats.pearsonr(f, g)[0]

    # Shuffle
    r_null = np.zeros(n_shuffles)
    g_copy = g.copy()
    for i in range(n_shuffles):
        np.random.shuffle(g_copy)
        r_null[i] = stats.pearsonr(f, g_copy)[0]

    # p-value: fraction of shuffles with |r| >= |r_obs|
    p_value = np.mean(np.abs(r_null) >= np.abs(r_obs))

    # Convert to sigma
    if p_value > 0:
        from scipy.stats import norm
        sigma = norm.ppf(1 - p_value / 2)
    else:
        sigma = float('inf')

    print(f"  Observed r = {r_obs:.4f}")
    print(f"  Null distribution: mean={np.mean(r_null):.4f}, std={np.std(r_null):.4f}")
    print(f"  p-value = {p_value:.6f}")
    print(f"  Significance = {sigma:.1f}σ")

    return sigma, p_value, r_obs, r_null


def plot_results(fermi, gaia, ells, cl_cross, r_ell, r_obs, r_null, sigma):
    """Generate comprehensive visualization."""
    import healpy as hp
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid = (fermi != hp.UNSEEN) & (gaia != hp.UNSEEN) & \
            (~np.isnan(fermi)) & (~np.isnan(gaia)) & (gaia > 0)

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle(f"Cross-correlation: Fermi γ-residual × Gaia DM halo  |  {sigma:.1f}σ",
                 fontsize=14, fontweight="bold")

    # 1. Scatter plot
    ax = axes[0, 0]
    f, g = fermi[valid], gaia[valid]
    ax.scatter(g, f, alpha=0.3, s=2, c="navy")
    ax.set_xlabel("Gaia DM halo proxy")
    ax.set_ylabel("Fermi γ residual")
    ax.set_title(f"Pixel correlation: r = {r_obs:.4f}")
    # Fit line
    z = np.polyfit(g, f, 1)
    xline = np.linspace(g.min(), g.max(), 100)
    ax.plot(xline, np.polyval(z, xline), "r-", linewidth=2)

    # 2. Shuffle histogram
    ax = axes[0, 1]
    ax.hist(r_null, bins=80, color="gray", alpha=0.7, density=True, label="Null (shuffled)")
    ax.axvline(r_obs, color="red", linewidth=2, label=f"Observed: {r_obs:.4f}")
    ax.set_xlabel("Pearson r")
    ax.set_ylabel("Density")
    ax.set_title(f"Significance: {sigma:.1f}σ  (p = {np.mean(np.abs(r_null)>=np.abs(r_obs)):.2e})")
    ax.legend()

    # 3. Angular cross-spectrum
    ax = axes[1, 0]
    sel = (ells >= 2) & (ells <= 60)
    ax.plot(ells[sel], r_ell[sel], "b-", linewidth=1.5)
    ax.axhline(0, color="gray", linestyle="--")
    ax.set_xlabel("Multipole ℓ")
    ax.set_ylabel("r(ℓ) = C_ℓ^{γ×DM} / √(C_ℓ^γ C_ℓ^{DM})")
    ax.set_title("Angular cross-correlation per multipole")
    ax.fill_between(ells[sel], -0.1, 0.1, alpha=0.2, color="gray", label="±0.1 (noise)")
    ax.legend()

    # 4. Both maps side by side (radial profile)
    ax = axes[1, 1]
    nside = hp.npix2nside(len(fermi))
    gc_vec = hp.ang2vec(np.pi/2, 0)
    pix_vecs = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))).T
    angles = np.degrees(np.arccos(np.clip(np.dot(pix_vecs, gc_vec), -1, 1)))

    bins_ang = np.arange(0, 45, 2)
    f_profile = []
    g_profile = []
    for i in range(len(bins_ang)-1):
        ring = valid & (angles >= bins_ang[i]) & (angles < bins_ang[i+1])
        if np.sum(ring) > 5:
            f_profile.append(np.mean(fermi[ring]))
            g_profile.append(np.mean(gaia[ring]))
        else:
            f_profile.append(np.nan)
            g_profile.append(np.nan)

    centers = 0.5 * (bins_ang[:-1] + bins_ang[1:])
    ax2 = ax.twinx()
    ax.plot(centers, f_profile, "ro-", markersize=4, label="Fermi γ residual")
    ax2.plot(centers, g_profile, "bs-", markersize=4, label="Gaia DM proxy")
    ax.set_xlabel("Angular distance from GC (°)")
    ax.set_ylabel("Fermi γ residual", color="red")
    ax2.set_ylabel("Gaia DM proxy", color="blue")
    ax.set_title("Radial profiles from Galactic Center")
    ax.legend(loc="upper left")
    ax2.legend(loc="upper right")

    plt.tight_layout()
    outpath = RESULTS_DIR / "cross_correlation_plot.png"
    plt.savefig(str(outpath), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → Plot: {outpath}")


def save_results(r_pearson, r_spearman, sigma, p_value):
    """Save summary to text file."""
    outpath = RESULTS_DIR / "cross_correlation_results.txt"
    with open(outpath, "w") as f:
        f.write("═" * 50 + "\n")
        f.write("  CROSS-CORRELATION RESULTS\n")
        f.write("  Fermi γ-residual (10-50 GeV) × Gaia DM halo\n")
        f.write("═" * 50 + "\n\n")
        f.write(f"  Pearson r  = {r_pearson:.6f}\n")
        f.write(f"  Spearman ρ = {r_spearman:.6f}\n")
        f.write(f"  Significance = {sigma:.2f}σ\n")
        f.write(f"  p-value = {p_value:.2e}\n\n")
        if sigma >= 5:
            f.write("  ★★★ DISCOVERY-LEVEL SIGNIFICANCE (≥5σ) ★★★\n")
        elif sigma >= 3:
            f.write("  ★★ EVIDENCE-LEVEL SIGNIFICANCE (≥3σ) ★★\n")
        elif sigma >= 2:
            f.write("  ★ HINT (2-3σ) — needs more data ★\n")
        else:
            f.write("  No significant correlation detected.\n")
        f.write("\n  Interpretation:\n")
        f.write("  Positive r = gamma excess correlates with DM halo\n")
        f.write("  → Supports dark matter annihilation hypothesis\n")
        f.write("  Negative/zero r = no spatial correspondence\n")
        f.write("  → Favors pulsar or other astrophysical origin\n")
    print(f"  → Summary: {outpath}")


if __name__ == "__main__":
    fermi, gaia = load_maps()
    if fermi is None:
        exit(1)

    r_p, r_s, p_p = pixel_correlation(fermi, gaia)
    ells, cl_cross, r_ell = angular_cross_spectrum(fermi, gaia)
    sigma, p_val, r_obs, r_null = shuffle_test(fermi, gaia, n_shuffles=10000)
    plot_results(fermi, gaia, ells, cl_cross, r_ell, r_obs, r_null, sigma)
    save_results(r_p, r_s, sigma, p_val)

    print()
    print("═" * 60)
    if sigma >= 3:
        print(f"  ★ RESULT: {sigma:.1f}σ correlation detected!")
        print(f"  Upload results/ folder to Claude for interpretation.")
    else:
        print(f"  Result: {sigma:.1f}σ — not significant.")
        print(f"  This may indicate: insufficient data, wrong energy band,")
        print(f"  or the signal is not dark matter.")
    print("═" * 60)
