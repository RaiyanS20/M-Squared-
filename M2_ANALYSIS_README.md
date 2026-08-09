# M² analysis of the fibre speckle image sets

`m2_analysis.py` computes the beam propagation ratio M² from the caustic z-scans
in `M-Squared Results/`, and compares the two step-index fibres (50 µm and
105 µm) under the two magnet configurations (Normal / non-alternating and
Alternating) at 0, 25 and 50 magnets.

## Running it

```bash
pip install numpy scipy matplotlib pillow
python3 m2_analysis.py                       # writes ./analysis_output
```

Useful flags:

| flag | meaning |
|---|---|
| `--wavelength-nm 532` | laser wavelength |
| `--z-step-mm 10` | camera translation between consecutive frames |
| `--pixel-size-um 3.45` | sensor pixel pitch |
| `--root "M-Squared Results"` | where the image folders live |
| `--out analysis_output` | output directory |
| `--jobs N` | parallel workers for the image processing |
| `--step-index-only` | skip the earlier runs that do not name the fibre |
| `--no-plots` | numbers only |
| `--no-cache` | ignore the cached per-frame widths and recompute |

The per-frame width extraction is the slow part (~4 min for all 1641 frames) and
is cached in `analysis_output/.widths_cache.json`, keyed by file mtime and by the
extraction constants. Re-running with a different wavelength or step size, or
just to redraw the plots, is therefore instant — only a change to the extraction
settings themselves forces a full recompute.

## Calibration constants — read this first

Three constants are **not** recoverable from the images and are supplied at the
top of the script:

```python
DEFAULT_WAVELENGTH_NM = 532.0    # laser wavelength
DEFAULT_Z_STEP_MM     = 10.0     # camera translation between consecutive frames
DEFAULT_PIXEL_SIZE_UM = 3.45     # Blackfly S / Sony IMX273 pixel pitch
```

The pixel pitch is well determined — the 1440×1080 and 1456×1088 frame sizes and
the "Blackfly" note in `Msquared script.mlx` identify a Sony IMX273 sensor.
The wavelength and the stage step were supplied by hand.

**M² scales linearly with the assumed step and inversely with the assumed
wavelength.** If either is wrong, every absolute M² moves by the same factor.
Because those constants are common to all scans, every *comparison* in the
report — fibre vs fibre, configuration vs configuration, magnet count vs magnet
count — is unaffected. Only the absolute scale would need correcting.

## What the script does

1. **Discovers z-scans.** Every leaf folder holding at least 5 numbered PNGs is
   one scan. Fibre, magnet count, configuration, trial number and polarisation
   arm are parsed from the folder names (including the spelling variants in the
   data: `Norrmal`, `Alternting`, `Magents`, `Norma`). Folders that do not name
   the fibre — the earlier `Normal Config … LP_n`, `Wollaston`, `singleMode` and
   `mutliMode_50um` runs — are still measured and written to the CSVs, but are
   excluded from the headline comparison and tagged `family=legacy`.

2. **Measures each frame** (ISO 11146 second moments):
   - background offset from the median of a border ring, subtracted;
   - the offset-corrected frame is **not** clipped at zero — rectifying the
     negative half of the noise puts a positive pedestal across the sensor, and
     since the second moment weights by distance squared that pedestal dominates
     for a beam this large (it also biases the wider sensor axis more than the
     narrow one, which appears as an x/y asymmetry that is not in the beam);
   - **secondary bright lobes are masked out.** Many frames contain a second
     spot — the other Wollaston arm, or a reflection near the border. A lobe far
     from the centroid dominates the second moment even when it carries little
     power. The frame is block-averaged 4×, smoothed, thresholded and labelled;
     only the brightest connected lobe is kept. This affects 320 of 1641 frames
     and was the single largest source of error in the first version of this
     analysis (it inflated some M² values by a factor of two to three);
   - moments are integrated over an elliptical aperture of 3× the D4σ diameters,
     iterated to convergence. When the beam is large enough that this aperture
     would run off the sensor, the aperture shrinks to the largest factor that
     still fits symmetrically about the centroid, so the integration region is
     never truncated on one side only. The factor actually used is recorded per
     frame as `aperture_factor`.

3. **Fits the caustic** d²(z) = a + bz + cz² per axis, and reports
   M² = (π/8λ)·√(4ac − b²), together with the waist diameter, waist position,
   Rayleigh length and far-field divergence.

   The fit uses an iterated 3σ robust residual cut (MAD-based, capped at 25% of
   the points). About 2.7% of frames are dropped this way. Every rejection is
   recorded per frame (`used_in_fit_x`, `used_in_fit_y`) and drawn as ✕ on the
   caustic plots, so nothing is silently discarded.

4. **Aggregates** over trials and over the H and V arms, and writes the report.

## Outputs

| file | contents |
|---|---|
| `per_image_widths.csv` | one row per frame: centroid, D4σ widths in px and µm, principal axes, aperture factor, stray-lobe count, rejection flags |
| `per_scan_fits.csv` | one row per z-scan: M²x, M²y with bootstrap errors, waist, Rayleigh length, divergence, R², ISO sampling flag |
| `summary_by_condition.csv` | trial-averaged M² per fibre / magnets / configuration / arm |
| `caustics_50um.png`, `caustics_105um.png` | measured caustics and fits |
| `m2_vs_magnets.png` | headline result |
| `m2_summary_bars.png` | trial-averaged M² per fibre and configuration |
| `waist_divergence.png` | waist size and divergence behind the M² |
| `example_beams.png` | near-waist profile for each condition |
| `M2_REPORT.md` | generated report and conclusion |

## Known limitation of the data

None of the scans fully meets the ISO 11146 sampling recommendation (≥10 points,
with ≥5 within one Rayleigh length of the waist **and** ≥5 beyond two Rayleigh
lengths). The waist sits near the start of the travel in every scan, roughly
20–40 mm in, while the Rayleigh length is around 40 mm — so the scans reach only
about 2 z_R on the far side and have almost nothing on the near side. The
hyperbola is therefore constrained mostly from one wing.

The fits are excellent (median R² = 0.999) and the trial-to-trial reproducibility
is good, so the *relative* comparison is sound. But the absolute M² carries a
systematic uncertainty larger than the quoted bootstrap error. **If you can
re-take data, extend the travel to roughly 2× the current range on the far side
and add three or four frames on the near side of the waist** — that alone would
tighten the absolute numbers considerably.
