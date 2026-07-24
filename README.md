Gaia DR3 × Totani (2025) — kinematic halo cross-correlation

Show Image

The result this repository reported is withdrawn. Versions 1.0 through 5.0 reported a positive Gaia × Fermi spatial cross-correlation and presented it as evidence for dark-matter annihilation, at significances quoted up to 9.4σ and, in v5.0, 6.1σ from a "frozen pipeline" together with a Bayes factor of 2×10³⁴. None of it survives correction. See below, and the corrected analysis linked at the end.

What was wrong

The working pixel set was conditioned on the target. The selection used throughout was

python
ok_r = (totani_64[idx] > 0) & ring & (g_vs > 0)

The first condition is an explicit positivity cut on the map being correlated. Of the pixels the stated coordinate rule defines, every one where the target flux is negative was excluded — 1782 of them, leaving a support containing no negative-target pixels at all. Every significance computed on that support, including the 6.1σ, is post-selection rather than calibrated inference.

The pipeline was not frozen. The Known Limitations section of this README already stated that there was no pre-registration and that five post-hoc decisions had been made. The Zenodo abstract described the same pipeline as frozen. The README was right.

The three "messengers" are not independent. The γ-ray residual, the 408 MHz synchrotron map and the stellar counts are all centrally concentrated and all share the Galactic latitude gradient; the synchrotron map alone has a rank correlation of 0.74 with |b|. Their mutual correlation is one degeneracy seen three times.

The low-energy bins do not exclude millisecond pulsars. Totani's fit assigns the halo component a negative amplitude at 1.5–2.5 GeV. The vanishing correlation there follows from that amplitude, not from anything about the source population.

The Bayes factor is withdrawn without replacement. A factor of 10³⁴ from nine spectral points indicates underestimated uncertainties, not evidence.

What replaces it

On a support generated from pixel coordinates alone — 7596 pixels, no map values consulted — the latitude-controlled partial rank correlation is +0.182, significant under a null whose type-I error is measured rather than assumed (p = 0.014). Controlling additionally for the angular distance from the Galactic Centre leaves +0.016 (p = 0.44), and the same collapse occurs for every independent stellar tracer tested. Across all nine of Totani's per-energy maps, no bin is significant after that control.

The replacement is not a smaller version of the same claim. It is a structural result about the method, which holds independently of these data:

For a spherically symmetric halo the line-of-sight annihilation (J-factor) morphology depends on direction only through the angle from the Galactic Centre, and decreases monotonically with it. Rank statistics are invariant under monotone transformation. A computed NFW² J-factor map therefore has Spearman correlation −1.000000 with that angle.

A rank cross-correlation between any tracer and a spherically symmetric γ-ray template is mathematically a test of central concentration, and cannot in principle isolate dark-matter-specific morphology.

Status of this code

Kept for provenance. It implements the withdrawn analysis and should not be used as-is. The corrected pipeline builds its mask from pixel coordinates alone, never touching map values, and calibrates its null empirically.

On the γ-ray excess itself

Nothing here bears on whether the Totani (2025) excess is real. It has since been independently reproduced at pixel level by Stenhouse, Ghag & Deppisch (arXiv:2607.08552). Whether it is dark matter remains contested on other grounds — a factor-of-a-few tension with dwarf-spheroidal limits, and a more severe one with AMS-02 antiproton data (Wang & Duan, JCAP 05 (2026) 087).

Corrected analysis

Paper, code, input maps, reproduction harness, and the responses to the two external audits that found the defects:

Zenodo: https://doi.org/10.5281/zenodo.21518911

Version history
Version	Date	DOI	Status
v1.0	2026-03-25	10.5281/zenodo.19221430	withdrawn
v2.0	2026-03-26	10.5281/zenodo.19236371	withdrawn
v3.0	2026-03-27	10.5281/zenodo.19244340	withdrawn
v4.0	2026-03-28	10.5281/zenodo.19268430	withdrawn
v5.0	2026-03-30	10.5281/zenodo.19325130	withdrawn
v6.0	2026-07-23	10.5281/zenodo.21518911	current — corrects all of the above
Citation

Cite the corrected version, not the withdrawn ones:

@misc{amser2026corrected,
  author    = {Amser, Olwyn},
  title     = {{Why a rank cross-correlation with a spherical template cannot
                test dark matter: a corrected Gaia x Fermi analysis}},
  year      = {2026},
  publisher = {Zenodo},
  version   = {v6.0},
  doi       = {10.5281/zenodo.21518911}
}
