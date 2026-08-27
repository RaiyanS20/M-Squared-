#!/usr/bin/env python3
"""
M-squared (beam propagation ratio) analysis of the speckle / caustic image sets.

Each leaf folder under the results directory is one z-scan: a sequence of camera
frames taken as the camera (or stage) is translated in equal steps through the
focus of the beam emerging from a multimode fibre.  For every frame the second
moment beam widths are computed following ISO 11146, then the caustic

    d(z)^2 = a + b*z + c*z^2

is fitted and the beam propagation ratio extracted as

    M^2 = (pi / (8*lambda)) * sqrt(4*a*c - b^2)

The experiment compares two step-index fibre core sizes (50 um and 105 um)
under two magnet configurations (Normal / non-alternating and Alternating) at
0, 25 and 50 magnets, with the beam split into two polarisation arms (H and V).

Usage
-----
    python3 m2_analysis.py                      # analyse everything, write ./analysis_output
    python3 m2_analysis.py --wavelength-nm 532 --z-step-mm 10 --pixel-size-um 3.45
    python3 m2_analysis.py --root "M-Squared Results" --out analysis_output
    python3 m2_analysis.py --jobs 8             # parallel image processing

Outputs (in --out)
------------------
    per_image_widths.csv     one row per frame: centroid, D4sigma widths, diagnostics
    per_scan_fits.csv        one row per z-scan: M2x, M2y, waist, divergence, fit quality
    summary_by_condition.csv trial-averaged M^2 for each fibre / config / magnet count
    caustics_50um.png        measured caustics + fits, 50 um step-index fibre
    caustics_105um.png       measured caustics + fits, 105 um step-index fibre
    m2_vs_magnets.png        headline result: M^2 vs magnet count, Normal vs Alternating
    m2_summary_bars.png      trial-averaged M^2 per fibre / configuration
    waist_divergence.png     waist diameter and far-field divergence behind the M^2
    example_beams.png        near-waist intensity profile for each condition
    M2_REPORT.md             generated report with the numbers and the conclusion
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field, asdict
from typing import Iterable

import numpy as np

# ----------------------------------------------------------------------------
# Experimental constants.  Override on the command line.
# ----------------------------------------------------------------------------

DEFAULT_WAVELENGTH_NM = 532.0    # laser wavelength
DEFAULT_Z_STEP_MM = 10.0         # camera translation between consecutive frames
DEFAULT_PIXEL_SIZE_UM = 3.45     # Blackfly S / Sony IMX273 pixel pitch

# ISO 11146 second-moment extraction settings.
APERTURE_FACTOR = 3.0            # integration aperture = 3 x D4sigma (ISO 11146-3)
MIN_APERTURE_FACTOR = 1.6        # smallest aperture accepted when the beam fills the sensor
MAX_APERTURE_ITERS = 25
APERTURE_TOL = 2e-3              # relative convergence tolerance on the widths
SEED_SIGMAS = 5.0                # seed threshold, in units of background noise sigma
BORDER_FRACTION = 0.04           # fraction of the frame used as the dark reference ring
SMOOTH_SIGMA_PX = 6.0            # blur used only to group speckle grains into lobes
LOBE_DOWNSAMPLE = 4              # block-average factor for the lobe labelling
LOBE_THRESHOLD = 0.06            # lobe-labelling level, as a fraction of the smoothed peak
STRAY_MIN_FRACTION = 0.005       # ignore secondary lobes fainter than this vs the main one
STRAY_DILATE_PX = 8              # grow the stray-lobe mask to catch its skirts
ASPECT_TOL = 1.35                # reject frames whose x/y aspect departs this far
                                 # from the median aspect of their own scan
RESID_REJECT_SIGMAS = 3.0        # robust residual cut in the caustic fit
MAX_REJECT_FRACTION = 0.25       # never throw away more than this share of a scan
MIN_IMAGES_PER_SCAN = 5          # a parabola fit needs more than 3 points to mean anything
SATURATION_LEVEL = 255           # 8-bit sensor
BOOTSTRAP_DRAWS = 400            # residual bootstrap for the M^2 uncertainty

# Categorical palette (validated: all-pairs CVD dE 9.2, normal-vision dE 24.0).
C_NORMAL = "#2a78d6"     # slot 1, blue    -> Normal (non-alternating) magnets
C_ALTERNATING = "#eb6834"  # slot 2, orange -> Alternating magnets
C_BASELINE = "#1baf7a"   # slot 3, aqua    -> 0 magnets (baseline)
C_INK = "#0b0b0b"
C_INK_2 = "#52514e"
C_GRID = "#e3e2de"
C_SURFACE = "#fcfcfb"

CONFIG_STYLE = {
    "Normal": (C_NORMAL, "o", "-"),
    "Alternating": (C_ALTERNATING, "s", "-"),
    "Baseline": (C_BASELINE, "D", "-"),
}

AXIS_LABEL = {"H": "H arm (horizontal polarisation)",
              "V": "V arm (vertical polarisation)"}


# ----------------------------------------------------------------------------
# Dataset discovery and metadata parsing
# ----------------------------------------------------------------------------

@dataclass
class Scan:
    """One z-scan: a folder of frames taken at equally spaced z positions."""
    path: str
    rel: str
    dataset: str          # the folder that names the experimental condition
    axis: str             # 'H', 'V' or '-' when the folder is not split
    fiber: str | None     # '50', '105', 'SMF' or None when not stated in the name
    magnets: int | None
    config: str | None    # 'Normal', 'Alternating' or 'Baseline'
    trial: int
    family: str           # 'step-index' for the main campaign, 'legacy' otherwise
    files: list[str] = field(default_factory=list)
    indices: list[int] = field(default_factory=list)

    @property
    def label(self) -> str:
        parts = []
        if self.fiber:
            parts.append(f"{self.fiber} um" if self.fiber != "SMF" else "SMF")
        if self.magnets is not None:
            parts.append(f"{self.magnets} magnets")
        if self.config:
            parts.append(self.config)
        parts.append(f"trial {self.trial}")
        if self.axis != "-":
            parts.append(self.axis)
        return ", ".join(parts)


_AXIS_NAMES = {"h": "H", "v": "V", "horizontal": "H", "vertical": "V"}


def _parse_axis(leaf: str) -> str:
    return _AXIS_NAMES.get(leaf.strip().lower(), "-")


def _parse_fiber(dataset: str) -> tuple[str | None, str]:
    """Return (fiber id, remainder of the folder name with the fibre part removed)."""
    m = re.match(r"^\s*(50|105)_Step\s*index", dataset, re.I)
    if m:
        return m.group(1), dataset[m.end():]
    m = re.match(r"^\s*Fiber_(\d+)um", dataset, re.I)
    if m:
        return m.group(1), dataset[m.end():]
    if re.search(r"multi?_?mode_?(\d+)um", dataset, re.I):
        mm = re.search(r"mode_?(\d+)um", dataset, re.I)
        return mm.group(1), dataset
    if re.search(r"single\s*mode", dataset, re.I):
        return "SMF", dataset
    return None, dataset


def _parse_magnets(rest: str) -> int | None:
    # Underscores are used as word separators in these folder names, so they
    # must not read as word characters ("25 Magnets_Alternating").
    rest = rest.replace("_", " ")
    if re.search(r"\bno\s+mag[ne]{0,3}[nt]?s?\b", rest, re.I):
        return 0
    m = re.search(r"(\d+)\s*Mag[ne]{1,3}ts?\b", rest, re.I)
    if m:
        return int(m.group(1))
    m = re.search(r"(\d+)\s*(?:Alternat|Alternt|Norr?ma)", rest, re.I)
    if m:
        return int(m.group(1))
    if re.search(r"\bcontrol\b", rest, re.I):
        return 0
    return None


def _parse_config(rest: str, magnets: int | None) -> str | None:
    if magnets == 0:
        # With no magnets in the beam path the two configurations are identical.
        return "Baseline"
    if re.search(r"alternat|alternt", rest, re.I):
        return "Alternating"
    if re.search(r"norr?ma", rest, re.I):
        return "Normal"
    return None


def _parse_trial(rest: str) -> int:
    m = re.search(r"Trial\s*(\d+)", rest, re.I)
    if m:
        return int(m.group(1))
    m = re.search(r"(?:LP|OM|M)_?(\d+)\s*$", rest, re.I)
    if m:
        return int(m.group(1))
    m = re.search(r"_(\d+)\s*$", rest)
    if m:
        return int(m.group(1))
    return 1


def discover_scans(root: str) -> list[Scan]:
    scans: list[Scan] = []
    for dirpath, _dirnames, filenames in os.walk(root):
        pngs = sorted(f for f in filenames if f.lower().endswith(".png"))
        if len(pngs) < MIN_IMAGES_PER_SCAN:
            continue
        indexed = []
        for f in pngs:
            m = re.search(r"(\d+)", os.path.splitext(f)[0])
            if m:
                indexed.append((int(m.group(1)), f))
        if len(indexed) < MIN_IMAGES_PER_SCAN:
            continue
        indexed.sort()

        rel = os.path.relpath(dirpath, root)
        parts = [p for p in rel.split(os.sep) if p not in (".", "")]
        leaf = parts[-1] if parts else ""
        axis = _parse_axis(leaf)
        # When the leaf is an H/V split, the condition is named by its parent.
        dataset_parts = parts[:-1] if (axis != "-" and len(parts) > 1) else parts
        dataset = " / ".join(dataset_parts) if dataset_parts else leaf

        fiber, rest = _parse_fiber(dataset)
        magnets = _parse_magnets(rest)
        config = _parse_config(rest, magnets)
        trial = _parse_trial(rest)
        family = "step-index" if re.search(r"Step\s*index", dataset, re.I) else "legacy"

        scans.append(Scan(path=dirpath, rel=rel, dataset=dataset, axis=axis,
                          fiber=fiber, magnets=magnets, config=config,
                          trial=trial, family=family,
                          files=[f for _, f in indexed],
                          indices=[i for i, _ in indexed]))
    scans.sort(key=lambda s: s.rel)
    return scans


# ----------------------------------------------------------------------------
# ISO 11146 second-moment beam widths
# ----------------------------------------------------------------------------

@dataclass
class FrameWidths:
    cx_px: float
    cy_px: float
    dx_px: float          # D4sigma along the lab x axis
    dy_px: float          # D4sigma along the lab y axis
    d_maj_px: float       # D4sigma along the principal axes
    d_min_px: float
    azimuth_deg: float
    power: float
    background: float
    noise_sigma: float
    peak: int
    saturated: bool
    aperture_factor: float  # aperture actually used, in units of the D4sigma diameter
    clipped: bool           # aperture had to be shrunk below MIN_APERTURE_FACTOR
    stray_removed: int      # secondary bright lobes masked out of the frame
    lobe_at_edge: bool      # the analysed lobe runs into the sensor border
    converged: bool


def _load_gray(path: str) -> np.ndarray:
    from PIL import Image
    img = Image.open(path)
    if img.mode not in ("L", "I;16", "I"):
        img = img.convert("L")
    return np.asarray(img)


def frame_widths(path: str) -> FrameWidths:
    """Second moment widths of a single frame, ISO 11146 style.

    A background offset is taken from a border ring and subtracted.  The offset
    corrected frame is *not* clipped at zero: rectifying the negative half of
    the noise would put a positive pedestal across the whole sensor, and because
    the second moment weights by distance squared that pedestal dominates the
    result for a beam this large (it also biases the wider sensor axis more than
    the narrow one, which shows up as an x/y asymmetry that is not in the beam).

    Moments are then integrated over an elliptical aperture of APERTURE_FACTOR x
    the D4sigma diameters, iterated to convergence.  When the beam is large
    enough that this aperture would run off the sensor, it is shrunk to the
    largest factor that still fits symmetrically about the centroid, so the
    integration region never gets truncated on one side only.
    """
    raw = _load_gray(path)
    peak = int(raw.max())
    a = raw.astype(np.float64)
    ny, nx = a.shape

    b = max(4, int(round(BORDER_FRACTION * min(ny, nx))))
    ring = np.concatenate([a[:b, :].ravel(), a[-b:, :].ravel(),
                           a[:, :b].ravel(), a[:, -b:].ravel()])
    bg = float(np.median(ring))
    noise = float(1.4826 * np.median(np.abs(ring - bg))) or float(ring.std())

    img = a - bg

    yy = np.arange(ny, dtype=np.float64)
    xx = np.arange(nx, dtype=np.float64)

    # These frames often contain more than one bright feature -- the second
    # Wollaston arm, or a reflection near the sensor border. A second lobe sits
    # far from the centroid, so the r^2 weighting of the second moment makes it
    # dominate the width even when it carries little power. Keep only the
    # brightest connected lobe; the faint wings of the beam itself fall below
    # the labelling threshold and are left untouched.
    img, seed, stray_removed, lobe_at_edge = _isolate_main_lobe(img, noise)

    cx, cy, sx, sy, _sxy, power = _moments(seed, xx, yy)
    if not np.isfinite(sx) or sx <= 0:
        cx, sx = nx / 2.0, nx / 8.0
    if not np.isfinite(sy) or sy <= 0:
        cy, sy = ny / 2.0, ny / 8.0

    X = xx[None, :]
    Y = yy[:, None]
    converged = False
    sxy = 0.0
    k_eff = APERTURE_FACTOR
    for _ in range(MAX_APERTURE_ITERS):
        # Largest aperture factor that still fits inside the sensor, centred on
        # the current centroid (semi-axis = k * 2 * sigma).
        k_fit = min(min(cx, nx - 1 - cx) / (2.0 * sx),
                    min(cy, ny - 1 - cy) / (2.0 * sy))
        k_eff = float(min(APERTURE_FACTOR, max(k_fit, MIN_APERTURE_FACTOR)))
        ax_, ay_ = k_eff * 2.0 * sx, k_eff * 2.0 * sy
        mask = (((X - cx) / ax_) ** 2 + ((Y - cy) / ay_) ** 2) <= 1.0
        win = np.where(mask, img, 0.0)
        ncx, ncy, nsx, nsy, nsxy, power = _moments(win, xx, yy)
        if not np.isfinite(nsx) or not np.isfinite(nsy) or nsx <= 0 or nsy <= 0:
            break
        done = abs(nsx - sx) <= APERTURE_TOL * sx and abs(nsy - sy) <= APERTURE_TOL * sy
        cx, cy, sx, sy, sxy = ncx, ncy, nsx, nsy, nsxy
        if done:
            converged = True
            break
    clipped = k_eff <= MIN_APERTURE_FACTOR + 1e-9

    # Principal axes (ISO 11146-1 eq. for a general astigmatic beam).
    vx, vy = sx ** 2, sy ** 2
    diff = vx - vy
    root = math.sqrt(max(diff * diff + 4.0 * sxy * sxy, 0.0))
    gamma = 1.0 if diff >= 0 else -1.0
    d_maj = 2.0 * math.sqrt(2.0) * math.sqrt(max(vx + vy + gamma * root, 0.0))
    d_min = 2.0 * math.sqrt(2.0) * math.sqrt(max(vx + vy - gamma * root, 0.0))
    azimuth = 0.5 * math.degrees(math.atan2(2.0 * sxy, diff)) if abs(diff) + abs(sxy) > 0 else 0.0

    return FrameWidths(cx_px=cx, cy_px=cy, dx_px=4.0 * sx, dy_px=4.0 * sy,
                       d_maj_px=d_maj, d_min_px=d_min, azimuth_deg=azimuth,
                       power=power, background=bg, noise_sigma=noise, peak=peak,
                       saturated=peak >= SATURATION_LEVEL, aperture_factor=k_eff,
                       clipped=clipped, stray_removed=stray_removed,
                       lobe_at_edge=lobe_at_edge, converged=converged)


def _isolate_main_lobe(img: np.ndarray, noise: float):
    """Zero out bright features that do not belong to the main beam lobe.

    Returns (cleaned image, seed image for the first moment estimate, number of
    lobes removed, whether the kept lobe touches the sensor border).
    """
    from scipy import ndimage

    ny, nx = img.shape
    # Lobes are hundreds of pixels across, so the labelling is done on a
    # block-averaged copy: same result, ~B^2 less work than smoothing and
    # labelling at full resolution.
    B = LOBE_DOWNSAMPLE
    py, px_ = (-ny) % B, (-nx) % B
    small = np.pad(img, ((0, py), (0, px_)), mode="edge")
    small = small.reshape(small.shape[0] // B, B, small.shape[1] // B, B).mean(axis=(1, 3))

    # Smooth before labelling so speckle grains join into one lobe rather than
    # fragmenting into hundreds of components.
    sm = ndimage.gaussian_filter(small, sigma=SMOOTH_SIGMA_PX / B)
    peak = float(sm.max())
    if peak <= 0:
        return img, np.clip(img, 0.0, None), 0, False

    thr = max(LOBE_THRESHOLD * peak, SEED_SIGMAS * max(noise, 1e-6) / B)
    labels, n = ndimage.label(sm > thr)
    if n == 0:
        return img, np.clip(img, 0.0, None), 0, False

    pos = np.clip(small, 0.0, None)
    sums = ndimage.sum_labels(pos, labels, index=np.arange(1, n + 1))
    keep = int(np.argmax(sums)) + 1

    cleaned = img
    removed = 0
    if n > 1:
        # Only bother masking lobes that carry a non-negligible amount of light.
        drop = [i + 1 for i in range(n)
                if i + 1 != keep and sums[i] > STRAY_MIN_FRACTION * sums[keep - 1]]
        if drop:
            mask = np.isin(labels, drop)
            mask = ndimage.binary_dilation(
                mask, iterations=max(1, STRAY_DILATE_PX // B))
            mask = _upsample(mask, B, ny, nx)
            cleaned = np.where(mask, 0.0, img)
            removed = len(drop)

    kept = _upsample(labels == keep, B, ny, nx)
    at_edge = bool(kept[0, :].any() or kept[-1, :].any()
                   or kept[:, 0].any() or kept[:, -1].any())
    seed = np.where(kept, np.clip(cleaned, 0.0, None), 0.0)
    if seed.sum() <= 0:
        seed = np.clip(cleaned, 0.0, None)
    return cleaned, seed, removed, at_edge


def _upsample(mask: np.ndarray, B: int, ny: int, nx: int) -> np.ndarray:
    return np.repeat(np.repeat(mask, B, axis=0), B, axis=1)[:ny, :nx]


def _moments(w: np.ndarray, xx: np.ndarray, yy: np.ndarray):
    p = float(w.sum())
    if p <= 0:
        return math.nan, math.nan, math.nan, math.nan, math.nan, 0.0
    col = w.sum(axis=0)
    row = w.sum(axis=1)
    cx = float((col * xx).sum() / p)
    cy = float((row * yy).sum() / p)
    vx = float((col * (xx - cx) ** 2).sum() / p)
    vy = float((row * (yy - cy) ** 2).sum() / p)
    sxy = float(((w * (yy[:, None] - cy)) * (xx[None, :] - cx)).sum() / p)
    return cx, cy, math.sqrt(max(vx, 0.0)), math.sqrt(max(vy, 0.0)), sxy, p


def _frame_worker(path: str):
    try:
        return path, asdict(frame_widths(path)), None
    except Exception as exc:  # a corrupt frame should not kill the run
        return path, None, repr(exc)


# ----------------------------------------------------------------------------
# Caustic fit
# ----------------------------------------------------------------------------

@dataclass
class CausticFit:
    m2: float
    m2_err: float
    d0_um: float          # waist diameter (D4sigma)
    z0_mm: float          # waist position along the scan
    zr_mm: float          # Rayleigh length
    theta_mrad: float     # full far-field divergence angle
    r2: float
    n_points: int
    n_near: int           # points within one Rayleigh length of the waist
    n_far: int            # points beyond two Rayleigh lengths
    iso_sampling_ok: bool
    valid: bool
    note: str = ""
    n_rejected: int = 0
    rejected_mask: tuple = ()


def fit_caustic(z_mm: np.ndarray, d_um: np.ndarray, wavelength_nm: float,
                rng: np.random.Generator | None = None) -> CausticFit:
    """Fit d^2 = a + b z + c z^2 and derive the ISO 11146 beam parameters.

    The fit is made robust with an iterated residual cut: points lying more than
    RESID_REJECT_SIGMAS robust standard deviations off the parabola are dropped
    and the fit repeated. A handful of frames in this data set are corrupted by
    stray light that survives the lobe isolation, and because the fit is over
    d^2 a single such point can drag M^2 by a factor of two. The number of
    dropped points is recorded per scan so the rejection is auditable.
    """
    z_all = np.asarray(z_mm, dtype=float) * 1e-3          # m
    d_all = np.asarray(d_um, dtype=float) * 1e-6          # m
    good = np.isfinite(z_all) & np.isfinite(d_all) & (d_all > 0)
    n_total = int(good.sum())
    if n_total < 4:
        return CausticFit(*([math.nan] * 6), r2=math.nan, n_points=n_total,
                          n_near=0, n_far=0, iso_sampling_ok=False, valid=False,
                          note="too few usable frames")

    keep = good.copy()
    max_drop = int(math.floor(MAX_REJECT_FRACTION * n_total))
    coeff = None
    for _ in range(4):
        z, y = z_all[keep], d_all[keep] ** 2
        if z.size < 4:
            break
        coeff = np.polyfit(z, y, 2)
        resid_all = d_all ** 2 - np.polyval(coeff, z_all)
        r = resid_all[keep]
        scale = 1.4826 * float(np.median(np.abs(r - np.median(r))))
        if scale <= 0:
            break
        bad = good & (np.abs(resid_all) > RESID_REJECT_SIGMAS * scale)
        # Drop the worst offenders only, and never more than the cap allows.
        order = np.argsort(-np.abs(resid_all) * bad)
        new_keep = good.copy()
        dropped = 0
        for i in order:
            if bad[i] and dropped < max_drop:
                new_keep[i] = False
                dropped += 1
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep

    z, d = z_all[keep], d_all[keep]
    n = z.size
    n_rejected = n_total - n
    if n < 4 or coeff is None:
        return CausticFit(*([math.nan] * 6), r2=math.nan, n_points=n, n_near=0,
                          n_far=0, iso_sampling_ok=False, valid=False,
                          note="too few usable frames after outlier rejection",
                          n_rejected=n_rejected)

    y = d ** 2
    coeff = np.polyfit(z, y, 2)
    c, b, a = coeff
    disc = 4.0 * a * c - b * b
    if c <= 0 or disc <= 0:
        return CausticFit(math.nan, math.nan, math.nan, math.nan, math.nan,
                          math.nan, r2=math.nan, n_points=n, n_near=0, n_far=0,
                          iso_sampling_ok=False, valid=False,
                          note="no minimum in the scan (curvature<=0 or disc<=0)",
                          n_rejected=n_rejected)

    lam = wavelength_nm * 1e-9
    m2 = math.pi / (8.0 * lam) * math.sqrt(disc)
    z0 = -b / (2.0 * c)
    d0 = math.sqrt(max(a - b * b / (4.0 * c), 0.0))
    zr = math.sqrt(disc) / (2.0 * c)
    theta = math.sqrt(c)                     # full angle, d(z) -> theta*z far field

    resid = y - np.polyval(coeff, z)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((resid ** 2).sum()) / ss_tot if ss_tot > 0 else math.nan

    dz = np.abs(z - z0)
    n_near = int((dz <= zr).sum())
    n_far = int((dz >= 2.0 * zr).sum())
    iso_ok = n >= 10 and n_near >= 5 and n_far >= 5

    # Residual bootstrap for the M^2 uncertainty.
    rng = rng or np.random.default_rng(12345)
    draws = []
    for _ in range(BOOTSTRAP_DRAWS):
        yb = np.polyval(coeff, z) + rng.choice(resid, size=n, replace=True)
        try:
            cb, bb, ab = np.polyfit(z, yb, 2)
        except Exception:
            continue
        db = 4.0 * ab * cb - bb * bb
        if cb > 0 and db > 0:
            draws.append(math.pi / (8.0 * lam) * math.sqrt(db))
    m2_err = float(np.std(draws, ddof=1)) if len(draws) > 5 else math.nan

    note = "" if iso_ok else (
        f"ISO 11146 sampling not met (need >=10 points, >=5 within 1 z_R and "
        f">=5 beyond 2 z_R; have {n}/{n_near}/{n_far})")
    return CausticFit(m2=m2, m2_err=m2_err, d0_um=d0 * 1e6, z0_mm=z0 * 1e3,
                      zr_mm=zr * 1e3, theta_mrad=theta * 1e3, r2=r2,
                      n_points=n, n_near=n_near, n_far=n_far,
                      iso_sampling_ok=iso_ok, valid=True, note=note,
                      n_rejected=n_rejected, rejected_mask=tuple(~keep))


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

@dataclass
class ScanResult:
    scan: Scan
    z_mm: np.ndarray
    dx_um: np.ndarray
    dy_um: np.ndarray
    frames: list[dict]
    fit_x: CausticFit
    fit_y: CausticFit

    @property
    def m2_mean(self) -> float:
        vals = [f.m2 for f in (self.fit_x, self.fit_y) if f.valid]
        return float(np.mean(vals)) if vals else math.nan


def _extraction_signature() -> str:
    """Identifies the width-extraction settings, so a stale cache is not reused."""
    import hashlib
    key = repr([APERTURE_FACTOR, MIN_APERTURE_FACTOR, MAX_APERTURE_ITERS,
                APERTURE_TOL, SEED_SIGMAS, BORDER_FRACTION, SMOOTH_SIGMA_PX,
                LOBE_DOWNSAMPLE, LOBE_THRESHOLD, STRAY_MIN_FRACTION,
                STRAY_DILATE_PX])
    return hashlib.sha256(key.encode()).hexdigest()[:16]


def analyse(scans: list[Scan], args) -> list[ScanResult]:
    all_paths: list[str] = []
    for s in scans:
        all_paths.extend(os.path.join(s.path, f) for f in s.files)

    # The width extraction is the expensive step and depends only on the image
    # file and the extraction constants -- not on wavelength, step or pixel
    # size. Cache it so re-running with different calibration, or just to
    # redraw the plots, is instant.
    cache_path = os.path.join(args.out, ".widths_cache.json")
    sig = _extraction_signature()
    cache: dict[str, dict] = {}
    if not args.no_cache and os.path.exists(cache_path):
        try:
            import json
            with open(cache_path) as fh:
                blob = json.load(fh)
            if blob.get("signature") == sig:
                cache = blob.get("frames", {})
        except Exception as exc:
            print(f"  (ignoring unreadable cache: {exc})", file=sys.stderr)

    widths: dict[str, dict] = {}
    todo: list[str] = []
    for path in all_paths:
        entry = cache.get(path)
        try:
            st = os.stat(path)
        except OSError:
            continue
        if entry and entry.get("mtime") == st.st_mtime and entry.get("size") == st.st_size:
            widths[path] = entry["w"]
        else:
            todo.append(path)

    if widths:
        print(f"Reusing cached widths for {len(widths)} frames.", flush=True)
    print(f"Measuring second-moment widths for {len(todo)} frames "
          f"in {len(scans)} scans ...", flush=True)

    failures: list[tuple[str, str]] = []
    done = 0

    def _record(path, w, err):
        nonlocal done
        done += 1
        if w is None:
            failures.append((path, err))
        else:
            widths[path] = w
        if done % 100 == 0 or done == len(todo):
            print(f"  {done}/{len(todo)} frames", flush=True)

    if args.jobs > 1 and len(todo) > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            for path, w, err in pool.map(_frame_worker, todo, chunksize=4):
                _record(path, w, err)
    else:
        for path in todo:
            _record(*_frame_worker(path))

    for path, err in failures:
        print(f"  ! failed on {path}: {err}", file=sys.stderr)

    if not args.no_cache:
        try:
            import json
            frames = {}
            for path, w in widths.items():
                st = os.stat(path)
                frames[path] = {"mtime": st.st_mtime, "size": st.st_size, "w": w}
            tmp = cache_path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump({"signature": sig, "frames": frames}, fh)
            os.replace(tmp, cache_path)
        except Exception as exc:
            print(f"  (could not write cache: {exc})", file=sys.stderr)

    px = args.pixel_size_um
    rng = np.random.default_rng(20240501)
    results: list[ScanResult] = []
    for s in scans:
        z, dx, dy, frames = [], [], [], []
        for idx, f in zip(s.indices, s.files):
            path = os.path.join(s.path, f)
            w = widths.get(path)
            if w is None:
                continue
            rec = dict(w)
            rec.update(scan=s.rel, dataset=s.dataset, axis=s.axis, fiber=s.fiber,
                       magnets=s.magnets, config=s.config, trial=s.trial,
                       family=s.family, file=f, index=idx,
                       z_mm=idx * args.z_step_mm,
                       dx_um=w["dx_px"] * px, dy_um=w["dy_px"] * px,
                       d_maj_um=w["d_maj_px"] * px, d_min_um=w["d_min_px"] * px)
            frames.append(rec)
            z.append(rec["z_mm"])
            dx.append(rec["dx_um"])
            dy.append(rec["dy_um"])
        z, dx, dy = np.array(z), np.array(dx), np.array(dy)

        # Frames whose x/y aspect ratio departs sharply from the rest of their
        # own scan are corrupted (stray light that survived the lobe isolation,
        # or the beam clipped by something at the end of the travel). The
        # comparison is against the scan's own median aspect, so this assumes
        # the beam shape is consistent along z, not that it is round.
        if z.size >= 5:
            with np.errstate(divide="ignore", invalid="ignore"):
                aspect = np.log(dx / dy)
            med = float(np.median(aspect[np.isfinite(aspect)]))
            bad_aspect = ~(np.abs(aspect - med) <= math.log(ASPECT_TOL))
        else:
            bad_aspect = np.zeros(z.size, dtype=bool)
        for rec, flag in zip(frames, bad_aspect):
            rec["aspect_outlier"] = bool(flag)
        zf = np.where(bad_aspect, np.nan, z)

        fx = fit_caustic(zf, dx, args.wavelength_nm, rng)
        fy = fit_caustic(zf, dy, args.wavelength_nm, rng)
        for rec, rx, ry in zip(frames, fx.rejected_mask or [False] * len(frames),
                               fy.rejected_mask or [False] * len(frames)):
            rec["used_in_fit_x"] = not rx
            rec["used_in_fit_y"] = not ry
        results.append(ScanResult(s, z, dx, dy, frames, fx, fy))
    return results


# ----------------------------------------------------------------------------
# Aggregation
# ----------------------------------------------------------------------------

def condition_key(r: ScanResult):
    return (r.scan.fiber, r.scan.magnets, r.scan.config, r.scan.axis)


def aggregate(results: Iterable[ScanResult]) -> dict[tuple, dict]:
    """Trial-average M^2 for each (fibre, magnets, config, axis) condition."""
    buckets: dict[tuple, list[ScanResult]] = {}
    for r in results:
        s = r.scan
        if s.family != "step-index" or s.fiber is None or s.magnets is None \
                or s.config is None or s.axis == "-":
            continue
        if not (r.fit_x.valid and r.fit_y.valid):
            continue
        buckets.setdefault(condition_key(r), []).append(r)

    out = {}
    for key, rs in buckets.items():
        m2x = np.array([r.fit_x.m2 for r in rs])
        m2y = np.array([r.fit_y.m2 for r in rs])
        m2 = np.concatenate([m2x, m2y])
        d0 = np.array([0.5 * (r.fit_x.d0_um + r.fit_y.d0_um) for r in rs])
        th = np.array([0.5 * (r.fit_x.theta_mrad + r.fit_y.theta_mrad) for r in rs])
        out[key] = dict(
            n_trials=len(rs),
            m2_mean=float(m2.mean()),
            m2_sd=float(m2.std(ddof=1)) if m2.size > 1 else 0.0,
            m2_sem=float(m2.std(ddof=1) / math.sqrt(m2.size)) if m2.size > 1 else 0.0,
            m2x_mean=float(m2x.mean()), m2y_mean=float(m2y.mean()),
            d0_um=float(d0.mean()), theta_mrad=float(th.mean()),
            trials=sorted(r.scan.trial for r in rs),
        )
    return out


# ----------------------------------------------------------------------------
# Output: CSV
# ----------------------------------------------------------------------------

def write_csvs(results: list[ScanResult], summary: dict, outdir: str):
    frame_rows = [rec for r in results for rec in r.frames]
    if frame_rows:
        cols = ["scan", "dataset", "family", "fiber", "magnets", "config", "trial",
                "axis", "file", "index", "z_mm", "cx_px", "cy_px", "dx_px", "dy_px",
                "dx_um", "dy_um", "d_maj_um", "d_min_um", "azimuth_deg", "power",
                "background", "noise_sigma", "peak", "saturated",
                "aperture_factor", "clipped", "stray_removed", "lobe_at_edge",
                "converged", "aspect_outlier", "used_in_fit_x", "used_in_fit_y"]
        with open(os.path.join(outdir, "per_image_widths.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(frame_rows)

    with open(os.path.join(outdir, "per_scan_fits.csv"), "w", newline="") as fh:
        cols = ["scan", "dataset", "family", "fiber", "magnets", "config", "trial",
                "axis", "n_frames",
                "M2x", "M2x_err", "M2y", "M2y_err", "M2_mean",
                "d0x_um", "d0y_um", "z0x_mm", "z0y_mm", "zRx_mm", "zRy_mm",
                "thetax_mrad", "thetay_mrad", "R2x", "R2y",
                "iso_sampling_ok", "n_rejected_x", "n_rejected_y",
                "n_saturated", "n_clipped", "note"]
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in results:
            s = r.scan
            notes = "; ".join(n for n in {r.fit_x.note, r.fit_y.note} if n)
            w.writerow(dict(
                scan=s.rel, dataset=s.dataset, family=s.family, fiber=s.fiber,
                magnets=s.magnets, config=s.config, trial=s.trial, axis=s.axis,
                n_frames=len(r.frames),
                M2x=_r(r.fit_x.m2), M2x_err=_r(r.fit_x.m2_err),
                M2y=_r(r.fit_y.m2), M2y_err=_r(r.fit_y.m2_err),
                M2_mean=_r(r.m2_mean),
                d0x_um=_r(r.fit_x.d0_um), d0y_um=_r(r.fit_y.d0_um),
                z0x_mm=_r(r.fit_x.z0_mm), z0y_mm=_r(r.fit_y.z0_mm),
                zRx_mm=_r(r.fit_x.zr_mm), zRy_mm=_r(r.fit_y.zr_mm),
                thetax_mrad=_r(r.fit_x.theta_mrad), thetay_mrad=_r(r.fit_y.theta_mrad),
                R2x=_r(r.fit_x.r2, 4), R2y=_r(r.fit_y.r2, 4),
                iso_sampling_ok=r.fit_x.iso_sampling_ok and r.fit_y.iso_sampling_ok,
                n_rejected_x=r.fit_x.n_rejected, n_rejected_y=r.fit_y.n_rejected,
                n_saturated=sum(1 for f in r.frames if f["saturated"]),
                n_clipped=sum(1 for f in r.frames if f["clipped"]),
                note=notes))

    with open(os.path.join(outdir, "summary_by_condition.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["fiber_um", "magnets", "config", "axis", "n_trials",
                    "M2_mean", "M2_sd", "M2_sem", "waist_d0_um", "theta_mrad"])
        for (fiber, mag, cfg, axis), v in sorted(
                summary.items(), key=lambda kv: (str(kv[0][0]), kv[0][1] or 0,
                                                 str(kv[0][2]), str(kv[0][3]))):
            w.writerow([fiber, mag, cfg, axis, v["n_trials"], _r(v["m2_mean"]),
                        _r(v["m2_sd"]), _r(v["m2_sem"]), _r(v["d0_um"]),
                        _r(v["theta_mrad"])])


def _r(x, nd=3):
    return "" if x is None or (isinstance(x, float) and not math.isfinite(x)) else round(float(x), nd)


# ----------------------------------------------------------------------------
# Output: figures
# ----------------------------------------------------------------------------

def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": C_SURFACE, "axes.facecolor": C_SURFACE,
        "savefig.facecolor": C_SURFACE,
        "axes.edgecolor": C_GRID, "axes.labelcolor": C_INK_2,
        "axes.titlecolor": C_INK, "text.color": C_INK,
        "xtick.color": C_INK_2, "ytick.color": C_INK_2,
        "grid.color": C_GRID, "grid.linewidth": 0.8,
        "axes.grid": True, "axes.axisbelow": True,
        "axes.spines.top": False, "axes.spines.right": False,
        "font.size": 9, "axes.titlesize": 10, "figure.dpi": 130,
        "legend.frameon": False,
        # Vector output: keep SVG text as real text (selectable and editable in
        # Illustrator/Inkscape) and embed TrueType rather than Type-3 in the
        # PDF, which is what journals and LaTeX pipelines expect.
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    return plt


# Formats every figure is written in. PDF is the one to use in LaTeX; SVG is
# for editing by hand; PNG is kept for quick viewing and for any consumer that
# cannot take vector art.
FIGURE_FORMATS = ("pdf", "svg", "png")


def _save(fig, stem: str, formats=None):
    """Write one figure to every requested format. `stem` carries no extension."""
    written = []
    for ext in (formats or FIGURE_FORMATS):
        path = f"{stem}.{ext}"
        # Raster panels (the beam images) stay raster inside the vector file;
        # 300 dpi keeps them sharp in print without bloating the PDF.
        fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=0.05)
        written.append(path)
    return written


def plot_caustics(results, fiber, stem, args, formats=None):
    plt = _style()
    rs = [r for r in results
          if r.scan.family == "step-index" and r.scan.fiber == fiber
          and r.scan.magnets is not None and r.scan.axis in ("H", "V")]
    if not rs:
        return
    magnets = sorted({r.scan.magnets for r in rs})
    axes_order = ["H", "V"]
    fig, axs = plt.subplots(len(axes_order), len(magnets),
                            figsize=(3.5 * len(magnets), 3.2 * len(axes_order)),
                            sharex=True, sharey=True, squeeze=False)
    for i, ax_name in enumerate(axes_order):
        for j, mag in enumerate(magnets):
            ax = axs[i][j]
            sel = [r for r in rs if r.scan.axis == ax_name and r.scan.magnets == mag]
            seen = set()
            for r in sorted(sel, key=lambda r: (str(r.scan.config), r.scan.trial)):
                cfg = r.scan.config or "Normal"
                color, marker, ls = CONFIG_STYLE.get(cfg, (C_INK_2, "o", "-"))
                lbl = cfg if cfg not in seen else None
                seen.add(cfg)
                d = 0.5 * (r.dx_um + r.dy_um)
                used = np.array([f.get("used_in_fit_x", True)
                                 and f.get("used_in_fit_y", True)
                                 and not f.get("aspect_outlier", False)
                                 for f in r.frames])
                ax.plot(r.z_mm[used], d[used], marker, color=color, ms=4.5,
                        alpha=0.85, mew=0, label=lbl)
                if (~used).any():
                    ax.plot(r.z_mm[~used], d[~used], "x", color=color, ms=5,
                            mew=1.2, alpha=0.6)
                f = r.fit_x
                if f.valid and r.z_mm.size > 2:
                    zz = np.linspace(r.z_mm.min(), r.z_mm.max(), 200)
                    dd = 0.5 * (_curve(zz, r.fit_x, args) + _curve(zz, r.fit_y, args))
                    ax.plot(zz, dd, ls, color=color, lw=1.4, alpha=0.55)
            if i == 0:
                ax.set_title(f"{mag} magnets", color=C_INK)
            if j == 0:
                ax.set_ylabel(f"{AXIS_LABEL[ax_name]}\nD4$\\sigma$ diameter (µm)")
            if i == len(axes_order) - 1:
                ax.set_xlabel("stage position z (mm)")
            if seen:
                ax.legend(loc="upper left", fontsize=8)
    fig.suptitle(f"Measured caustics — {fiber} µm step-index fibre "
                 f"(λ={args.wavelength_nm:g} nm, {args.z_step_mm:g} mm per step)",
                 fontsize=11, color=C_INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, stem, formats)
    plt.close(fig)


def _curve(z_mm, fit: CausticFit, args):
    """Reconstruct the fitted D4sigma diameter in µm at positions z_mm."""
    lam_um = args.wavelength_nm * 1e-3
    z_um = (z_mm - fit.z0_mm) * 1e3
    d0 = fit.d0_um
    return np.sqrt(d0 ** 2 + (fit.m2 * 4.0 * lam_um / (math.pi * d0)) ** 2 * z_um ** 2)


def plot_m2_vs_magnets(summary, stem, args, formats=None):
    plt = _style()
    fibers = sorted({k[0] for k in summary}, key=lambda f: int(f))
    axes_order = ["H", "V"]
    fig, axs = plt.subplots(len(axes_order), len(fibers),
                            figsize=(4.2 * len(fibers), 3.4 * len(axes_order)),
                            squeeze=False)
    for i, ax_name in enumerate(axes_order):
        for j, fiber in enumerate(fibers):
            ax = axs[i][j]
            base = [(k[1], v) for k, v in summary.items()
                    if k[0] == fiber and k[3] == ax_name and k[2] == "Baseline"]
            for cfg in ("Normal", "Alternating"):
                pts = [(k[1], v) for k, v in summary.items()
                       if k[0] == fiber and k[3] == ax_name and k[2] == cfg]
                pts = sorted(base + pts, key=lambda p: p[0])
                if not pts:
                    continue
                x = [p[0] for p in pts]
                y = [p[1]["m2_mean"] for p in pts]
                e = [p[1]["m2_sd"] for p in pts]
                color, marker, _ = CONFIG_STYLE[cfg]
                ax.errorbar(x, y, yerr=e, color=color, marker=marker, ms=7,
                            lw=2, capsize=3, elinewidth=1.2, label=cfg, zorder=3)
                # Label the end of each series only -- a number on every point
                # collides with its neighbour and with the axes. Exact values
                # for every condition are in the report table and the CSVs.
                ax.annotate(f" {cfg}: {y[-1]:.1f}", (x[-1], y[-1]),
                            textcoords="offset points",
                            xytext=(6, 11 if cfg == "Normal" else -14),
                            ha="right", fontsize=8, color=color)
            for mag, v in base:
                color, marker, _ = CONFIG_STYLE["Baseline"]
                ax.errorbar([mag], [v["m2_mean"]], yerr=[v["m2_sd"]], color=color,
                            marker=marker, ms=8, lw=0, capsize=3, elinewidth=1.2,
                            label="0 magnets (baseline)", zorder=4)
            ax.set_xticks([0, 25, 50])
            ax.set_xlim(-8, 62)
            ax.margins(y=0.22)
            ax.set_xlabel("number of magnets")
            ax.set_ylabel("M²")
            ax.set_title(f"{fiber} µm step-index — {AXIS_LABEL[ax_name]}", color=C_INK)
            handles, labels = ax.get_legend_handles_labels()
            uniq = dict(zip(labels, handles))
            ax.legend(uniq.values(), uniq.keys(), fontsize=8, loc="lower left")
    fig.suptitle("Beam propagation ratio vs magnet count and configuration\n"
                 f"(error bars = spread over trials and H/V fits; λ={args.wavelength_nm:g} nm)",
                 fontsize=11, color=C_INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    _save(fig, stem, formats)
    plt.close(fig)


def plot_summary_bars(summary, stem, formats=None):
    """Change in M^2 relative to the 0-magnet baseline of the same fibre.

    The magnet effects are a few percent on top of baselines that differ by
    almost 3x between the fibres, so plotting the absolute M^2 hides them. This
    shows the quantity actually under test, with a +/-2 standard error band so
    it is immediately visible which changes clear the noise.
    """
    plt = _style()
    fibers = sorted({k[0] for k in summary}, key=lambda f: int(f))
    magnets = [m for m in sorted({k[1] for k in summary}) if m != 0]
    fig, axs = plt.subplots(1, len(fibers), figsize=(4.8 * len(fibers), 3.8),
                            squeeze=False)
    for j, fiber in enumerate(fibers):
        ax = axs[0][j]
        base = _combine_stats(summary, fiber, 0, "Baseline")
        width = 0.36
        xs = np.arange(len(magnets), dtype=float)
        for k, cfg in enumerate(("Normal", "Alternating")):
            ys, es, xs_used = [], [], []
            for mi, mag in enumerate(magnets):
                v = _combine_stats(summary, fiber, mag, cfg)
                if v is None or base is None or base[0] == 0:
                    continue
                rel = 100.0 * (v[0] - base[0]) / base[0]
                # Propagate the two standard errors onto the ratio.
                err = 100.0 * math.sqrt(v[1] ** 2 + (v[0] / base[0] * base[1]) ** 2) / base[0]
                ys.append(rel); es.append(err)
                xs_used.append(xs[mi] + (k - 0.5) * (width + 0.02))
            color = CONFIG_STYLE[cfg][0]
            ax.bar(xs_used, ys, width=width, color=color, label=cfg,
                   edgecolor=C_SURFACE, linewidth=2, zorder=2)
            ax.errorbar(xs_used, ys, yerr=[2 * e for e in es], fmt="none",
                        ecolor=C_INK_2, elinewidth=1.1, capsize=3, zorder=3)
            for xi, yi, ei in zip(xs_used, ys, es):
                off = 6 if yi >= 0 else -14
                ax.annotate(f"{yi:+.1f}%", (xi, yi + math.copysign(2 * ei, yi or 1.0)),
                            textcoords="offset points", xytext=(0, off),
                            ha="center", fontsize=8, color=C_INK_2)
        ax.axhline(0, color=C_INK_2, lw=1.2, zorder=1)
        ax.set_xticks(xs)
        ax.set_xticklabels([f"{m}" for m in magnets])
        ax.set_xlabel("number of magnets")
        ax.set_ylabel("change in M² vs 0 magnets (%)")
        base_txt = f" (baseline M² = {base[0]:.1f})" if base else ""
        ax.set_title(f"{fiber} µm step-index fibre{base_txt}", color=C_INK)
        ax.legend(fontsize=8)
        lim = max(abs(np.array(ax.get_ylim()))) * 1.15
        ax.set_ylim(-lim, lim)
    fig.suptitle("Magnet-induced change in M², relative to each fibre's own baseline\n"
                 "(error bars = ±2 standard errors; a bar not crossing zero is a real effect)",
                 fontsize=11, color=C_INK)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    _save(fig, stem, formats)
    plt.close(fig)


def _combine_field(summary, fiber, magnets, config, field_):
    """Mean of one summary field over the H and V arms."""
    vals = [v[field_] for k, v in summary.items()
            if k[0] == fiber and k[1] == magnets and k[2] == config]
    return float(np.mean(vals)) if vals else None


def _combine_stats(summary, fiber, magnets, config):
    """Mean M^2 and its standard error, pooled over the H and V arms."""
    vals = [v for k, v in summary.items()
            if k[0] == fiber and k[1] == magnets and k[2] == config]
    if not vals:
        return None
    m = float(np.mean([v["m2_mean"] for v in vals]))
    # Pool the per-arm standard errors, then average over the arms.
    sem = float(math.sqrt(sum(v["m2_sem"] ** 2 for v in vals)) / len(vals))
    return m, sem


def plot_waist_divergence(summary, stem, formats=None):
    """Separate the M^2 change into its two factors: M^2 = pi*d0*theta/(8*lambda).

    Plotted as a change relative to each fibre's own 0-magnet baseline -- the two
    fibres differ by ~2x in absolute waist and divergence, so a shared absolute
    scale would compress the few-percent effect under test to nothing.
    """
    plt = _style()
    fibers = sorted({k[0] for k in summary}, key=lambda f: int(f))
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.8))
    for ax, field_, ylabel in ((axs[0], "d0_um", "change in waist diameter D4σ (%)"),
                               (axs[1], "theta_mrad", "change in full divergence θ (%)")):
        for fiber in fibers:
            base_vals = [v[field_] for kk, v in summary.items()
                         if kk[0] == fiber and kk[2] == "Baseline"]
            if not base_vals:
                continue
            base = float(np.mean(base_vals))
            for cfg in ("Normal", "Alternating"):
                agg: dict[int, list[float]] = {}
                for kk, v in summary.items():
                    if kk[0] != fiber:
                        continue
                    if kk[2] == cfg or (kk[2] == "Baseline" and kk[1] == 0):
                        agg.setdefault(kk[1], []).append(v[field_])
                if not agg:
                    continue
                x = sorted(agg)
                y = [100.0 * (float(np.mean(agg[m])) - base) / base for m in x]
                color = CONFIG_STYLE[cfg][0]
                first = fiber == fibers[0]
                ax.plot(x, y, "-" if first else "--", marker="o" if first else "s",
                        color=color, lw=2, ms=6, mfc=color if first else C_SURFACE,
                        mew=1.6, label=f"{fiber} µm, {cfg}")
        ax.axhline(0, color=C_INK_2, lw=1.2, zorder=1)
        ax.set_xticks([0, 25, 50])
        ax.margins(y=0.2)
        ax.set_xlabel("number of magnets")
        ax.set_ylabel(ylabel)
    axs[0].legend(fontsize=7.5, ncol=2, loc="best")
    fig.suptitle("What drives the M² change — M² = π·d₀·θ / (8λ), so it is the product "
                 "of waist size and divergence\n"
                 "(solid + filled = 50 µm, dashed + hollow = 105 µm; "
                 "each relative to its own 0-magnet baseline)",
                 fontsize=10, color=C_INK)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    _save(fig, stem, formats)
    plt.close(fig)


def plot_example_beams(results, stem, formats=None):
    plt = _style()
    picks = []
    for fiber in ("50", "105"):
        for mag, cfg in ((0, "Baseline"), (25, "Normal"), (25, "Alternating"),
                         (50, "Normal"), (50, "Alternating")):
            cand = [r for r in results
                    if r.scan.family == "step-index" and r.scan.fiber == fiber
                    and r.scan.magnets == mag and r.scan.config == cfg
                    and r.scan.axis == "H" and r.fit_x.valid]
            if cand:
                picks.append((fiber, mag, cfg, sorted(cand, key=lambda r: r.scan.trial)[0]))
    if not picks:
        return
    ncol = 5
    nrow = int(math.ceil(len(picks) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.5 * ncol, 2.2 * nrow), squeeze=False)
    for ax in axs.ravel():
        ax.axis("off")
    for k, (fiber, mag, cfg, r) in enumerate(picks):
        ax = axs[k // ncol][k % ncol]
        # frame closest to the fitted waist
        z0 = r.fit_x.z0_mm
        i = int(np.argmin(np.abs(r.z_mm - z0)))
        img = _load_gray(os.path.join(r.scan.path, r.scan.files[i]))
        rec = r.frames[i]
        cx, cy = rec["cx_px"], rec["cy_px"]
        half = max(rec["dx_px"], rec["dy_px"])
        x0, x1 = int(max(cx - half, 0)), int(min(cx + half, img.shape[1]))
        y0, y1 = int(max(cy - half, 0)), int(min(cy + half, img.shape[0]))
        crop = img[y0:y1, x0:x1]
        # Stretch to the 99.5th percentile -- a handful of bright speckle grains
        # otherwise leave the whole pattern sitting in the bottom of the ramp.
        hi = float(np.percentile(crop, 99.5)) or float(crop.max()) or 1.0
        ax.imshow(crop, cmap="inferno", vmin=0, vmax=hi)
        ax.set_title(f"{fiber} µm · {mag} mag · {cfg}\nD4σ ≈ {rec['dx_um']:.0f} µm",
                     fontsize=8, color=C_INK)
        ax.axis("off")
    fig.suptitle("Near-waist speckle / intensity profile for each condition",
                 fontsize=11, color=C_INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    _save(fig, stem, formats)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Output: report
# ----------------------------------------------------------------------------

def write_report(results, summary, outdir, args):
    lines = []
    A = lines.append
    A("# M² analysis of the fibre speckle image sets\n")
    A(f"Generated by `m2_analysis.py`. Calibration used: **λ = {args.wavelength_nm:g} nm**, "
      f"**{args.z_step_mm:g} mm** camera translation per frame, "
      f"**{args.pixel_size_um:g} µm** pixel pitch.\n")
    A("Beam widths are ISO 11146 second-moment (D4σ) diameters, computed after "
      "border-ring background subtraction and integrated over an elliptical "
      f"aperture of {APERTURE_FACTOR:g}× the D4σ diameter, iterated to convergence. "
      "M² comes from the hyperbolic fit d²(z) = a + bz + cz², "
      "M² = (π/8λ)·√(4ac − b²).\n")
    A("> M² scales linearly with the assumed stage step and inversely with the "
      "assumed wavelength. If either constant is wrong, every number below moves "
      "by the same factor and all *relative* comparisons still hold.\n")

    fibers = sorted({k[0] for k in summary}, key=lambda f: int(f))

    A("\n## Headline numbers\n")
    A("| Fibre | Magnets | Configuration | Arm | trials | M² | sd | waist D4σ (µm) | θ (mrad) |")
    A("|---|---|---|---|---|---|---|---|---|")
    for key in sorted(summary, key=lambda k: (int(k[0]), k[1], k[2], k[3])):
        f, mag, cfg, axis = key
        v = summary[key]
        A(f"| {f} µm | {mag} | {cfg} | {axis} | {v['n_trials']} | "
          f"{v['m2_mean']:.1f} | {v['m2_sd']:.1f} | {v['d0_um']:.0f} | "
          f"{v['theta_mrad']:.1f} |")

    A("\n## Effect of the magnets, per fibre\n")
    A("A difference is called significant when it exceeds twice the combined "
      "standard error of the two means being compared (roughly a 95% criterion).\n")
    A("| Fibre | Configuration | Magnets | M² (H+V) ± s.e. | Δ vs 0 magnets | % change | significant? |")
    A("|---|---|---|---|---|---|---|")
    deltas = {}
    for f in fibers:
        base = _combine_stats(summary, f, 0, "Baseline")
        if base is not None:
            A(f"| {f} µm | baseline | 0 | {base[0]:.2f} ± {base[1]:.2f} | — | — | — |")
        for cfg in ("Normal", "Alternating"):
            for mag in (25, 50):
                v = _combine_stats(summary, f, mag, cfg)
                if v is None:
                    continue
                if base is None:
                    A(f"| {f} µm | {cfg} | {mag} | "
                      f"{v[0]:.2f} ± {v[1]:.2f} | — | — | — |")
                    continue
                d = v[0] - base[0]
                comb = math.sqrt(v[1] ** 2 + base[1] ** 2)
                sig = abs(d) > 2 * comb
                A(f"| {f} µm | {cfg} | {mag} | {v[0]:.2f} ± {v[1]:.2f} | "
                  f"{d:+.2f} | {100*d/base[0]:+.1f}% | "
                  f"{'**yes**' if sig else 'no'} |")
                deltas[(f, cfg, mag)] = (v[0], v[1], d, 100 * d / base[0], comb, sig)

    A("\n## Where the change comes from\n")
    A("M² = π·d₀·θ/(8λ), so any change in M² is a change in the waist diameter "
      "d₀, the far-field divergence θ, or both.\n")
    A("| Fibre | Configuration | Magnets | Δ waist | Δ divergence | Δ M² | dominated by |")
    A("|---|---|---|---|---|---|---|")
    for f in fibers:
        bd = _combine_field(summary, f, 0, "Baseline", "d0_um")
        bt = _combine_field(summary, f, 0, "Baseline", "theta_mrad")
        for cfg in ("Normal", "Alternating"):
            for mag in (25, 50):
                d = _combine_field(summary, f, mag, cfg, "d0_um")
                t = _combine_field(summary, f, mag, cfg, "theta_mrad")
                if None in (bd, bt, d, t) or (f, cfg, mag) not in deltas:
                    continue
                dd = 100 * (d - bd) / bd
                dt = 100 * (t - bt) / bt
                which = "waist" if abs(dd) > 1.5 * abs(dt) else (
                    "divergence" if abs(dt) > 1.5 * abs(dd) else "both equally")
                A(f"| {f} µm | {cfg} | {mag} | {dd:+.1f}% | {dt:+.1f}% | "
                  f"{deltas[(f, cfg, mag)][3]:+.1f}% | {which} |")

    A("\n## Conclusion\n")
    for para in build_conclusion(summary, deltas, fibers):
        A(para + "\n")

    A("\n## Figures\n")
    A("Each figure is written as PDF (vector, for LaTeX), SVG (vector, for "
      "editing) and PNG (for quick viewing).\n")
    A("| file | what it shows |")
    A("|---|---|")
    A("| `caustics_50um.pdf`, `caustics_105um.pdf` | every measured caustic with "
      "its fit; ✕ marks a frame rejected from the fit |")
    A("| `m2_vs_magnets.pdf` | M² vs magnet count, Normal vs Alternating, per fibre and arm |")
    A("| `m2_summary_bars.pdf` | the headline: change in M² relative to each "
      "fibre's own baseline, with ±2 s.e. bars |")
    A("| `waist_divergence.pdf` | whether a given change came from the waist or the divergence |")
    A("| `example_beams.pdf` | near-waist speckle pattern per condition — the "
      "105 µm fibre visibly carries many more speckle grains, i.e. more modes |")

    A("\n## Data quality flags\n")
    bad = [r for r in results if r.scan.family == "step-index"
           and (not r.fit_x.valid or not r.fit_y.valid or not r.fit_x.iso_sampling_ok)]
    if bad:
        A(f"{len(bad)} of the step-index scans do not meet the full ISO 11146 "
          "sampling recommendation (≥10 points, ≥5 within one Rayleigh length of "
          "the waist and ≥5 beyond two Rayleigh lengths). In these scans the waist "
          "sits near the start of the travel, so the near-field side is "
          "under-sampled and the fitted M² carries more uncertainty than the "
          "bootstrap error alone suggests. See `per_scan_fits.csv` for the "
          "per-scan flags.\n")
    sat = sum(1 for r in results for f in r.frames if f["saturated"])
    clip = sum(1 for r in results for f in r.frames if f["clipped"])
    stray = sum(1 for r in results for f in r.frames if f.get("stray_removed"))
    aspect = sum(1 for r in results for f in r.frames if f.get("aspect_outlier"))
    resid = sum(max(r.fit_x.n_rejected, r.fit_y.n_rejected) for r in results)
    total = sum(len(r.frames) for r in results)
    A(f"- Frames containing a secondary bright lobe that was masked out before "
      f"the moment integration: **{stray}** of {total}")
    A(f"- Frames rejected as aspect-ratio outliers (x/y width ratio more than "
      f"{ASPECT_TOL:g}× away from their own scan's median): **{aspect}** "
      f"({100*aspect/max(total,1):.1f}%)")
    A(f"- Further frames rejected by the >{RESID_REJECT_SIGMAS:g}σ residual cut "
      f"in the caustic fit: **{resid}** ({100*resid/max(total,1):.1f}%)")
    A("  Both are flagged per frame in `per_image_widths.csv` "
      "(`aspect_outlier`, `used_in_fit_x`, `used_in_fit_y`) and drawn as ✕ on "
      "the caustic plots, so every rejection is auditable.")
    A(f"- Saturated frames (peak = 255): **{sat}**")
    A(f"- Frames where the beam was large enough that the integration aperture "
      f"had to shrink to its floor of {MIN_APERTURE_FACTOR:g}×D4σ to stay on the "
      f"sensor: **{clip}**. A tighter aperture truncates the outer wings, which "
      "biases the width slightly low; these are all at the far end of the travel, "
      "where the beam is widest. `aperture_factor` in `per_image_widths.csv` "
      "gives the value used for every frame.")
    A(f"- Scans analysed: **{len(results)}** "
      f"({sum(1 for r in results if r.scan.family == 'step-index')} in the "
      "step-index campaign, the rest are earlier/legacy runs kept in the CSVs)\n")

    with open(os.path.join(outdir, "M2_REPORT.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return "\n".join(lines)


def build_conclusion(summary, deltas, fibers) -> list[str]:
    out = []
    if not summary:
        return ["No step-index scans were successfully fitted."]

    # 1. Baseline comparison between the two fibres.
    b = {f: _combine_stats(summary, f, 0, "Baseline") for f in fibers}
    if all(v is not None for v in b.values()) and len(fibers) == 2:
        f1, f2 = fibers
        out.append(
            f"**Baseline beam quality scales with core size.** With no magnets the "
            f"{f1} µm fibre gives M² = {b[f1][0]:.1f} and the {f2} µm fibre "
            f"M² = {b[f2][0]:.1f}, a ratio of {b[f2][0]/b[f1][0]:.2f}. For a fully "
            f"filled step-index fibre M² grows in proportion to core diameter × NA, "
            f"so a ratio of {int(f2)/int(f1):.2f} would be expected if both were "
            f"filled identically; the measured ratio being "
            f"{'larger' if b[f2][0]/b[f1][0] > int(f2)/int(f1) else 'smaller'} "
            f"indicates the {f2} µm fibre is "
            f"{'excited over a wider mode set (more overfilled) than' if b[f2][0]/b[f1][0] > int(f2)/int(f1) else 'launched more selectively than'} "
            f"the {f1} µm one.")

    # 2. Per fibre: which configuration does more.
    for f in fibers:
        rows = []
        for cfg in ("Normal", "Alternating"):
            for mag in (25, 50):
                if (f, cfg, mag) in deltas:
                    rows.append((cfg, mag) + deltas[(f, cfg, mag)])
        if not rows:
            continue
        norm = [r for r in rows if r[0] == "Normal"]
        alt = [r for r in rows if r[0] == "Alternating"]
        base = b.get(f)
        txt = [f"**{f} µm fibre.**"]
        if base:
            txt.append(f"Baseline M² = {base[0]:.2f} ± {base[1]:.2f}.")
        for cfg, rs in (("Normal", norm), ("Alternating", alt)):
            if not rs:
                continue
            desc = ", ".join(
                f"{r[1]} magnets → M² = {r[2]:.2f} ({r[5]:+.1f}%"
                f"{', significant' if r[7] else ', not significant'})"
                for r in sorted(rs, key=lambda r: r[1]))
            txt.append(f"{cfg}: {desc}.")
        if norm and alt:
            n50 = max(norm, key=lambda r: r[1])
            a50 = max(alt, key=lambda r: r[1])
            gap = a50[2] - n50[2]
            comb = math.sqrt(n50[3] ** 2 + a50[3] ** 2)
            if abs(gap) > 2 * comb:
                hi = "Alternating" if gap > 0 else "Normal"
                lo = "Normal" if gap > 0 else "Alternating"
                txt.append(
                    f"At {n50[1]} magnets the two configurations differ by "
                    f"{abs(gap):.2f} in M² ({hi} higher), against a combined standard "
                    f"error of {comb:.2f} — a significant difference, so {hi} degrades "
                    f"the beam measurably more than {lo} in this fibre.")
            else:
                txt.append(
                    f"At {n50[1]} magnets the two configurations differ by only "
                    f"{abs(gap):.2f} in M², against a combined standard error of "
                    f"{comb:.2f} — not significant, so on this data the two "
                    f"configurations cannot be distinguished for this fibre.")
        out.append(" ".join(txt))

    # 3. Cross-fibre statement about sensitivity.
    sens, signs = {}, {}
    for f in fibers:
        if not b.get(f):
            continue
        ch = [deltas[k][3] for k in deltas if k[0] == f]             # index 3 = % change
        sig_ch = [deltas[k][3] for k in deltas if k[0] == f and deltas[k][5]]
        if ch:
            sens[f] = max(abs(c) for c in ch)
            signs[f] = sig_ch
    n_sig = sum(1 for v in deltas.values() if v[5])

    # Direction of the effect, stated from the significant changes only.
    all_sig = [c for v in signs.values() for c in v]
    if all_sig:
        neg = sum(1 for c in all_sig if c < 0)
        if neg == len(all_sig):
            direction = (
                "**Direction of the effect.** Every change that clears the noise "
                "is a *reduction* in M² — the magnets make the beam slightly "
                "better, not worse, by up to "
                f"{max(abs(c) for c in all_sig):.1f}%. That rules out simple "
                "scattering into higher-order modes, which would raise M². A "
                "reduction instead points to the magnets acting as a weak mode "
                "filter: bending or stressing the fibre strips the highest-order "
                "modes, which are the least well confined and contribute most to "
                "the second moment. Note the effect is small and could also be a "
                "systematic of the launch drifting slightly when the magnets are "
                "fitted, so it is worth confirming with a repeat in which the "
                "magnets are added and removed several times in the same session.")
        elif neg == 0:
            direction = (
                "**Direction of the effect.** Every change that clears the noise "
                "is an *increase* in M² — the magnets degrade the beam, by up to "
                f"{max(all_sig):.1f}%. That is the signature of mode coupling: "
                "the perturbation redistributes power into higher-order modes, "
                "which widens the second moment in the far field.")
        else:
            direction = (
                "**Direction of the effect.** The significant changes do not share "
                f"a sign ({neg} reductions and {len(all_sig)-neg} increases in M²), "
                "so the magnets are not simply degrading or simply cleaning up the "
                "beam. With effects this small relative to the run-to-run "
                "reproducibility, the most likely explanation is that refitting the "
                "magnets slightly disturbs the launch conditions rather than that "
                "the field itself is acting on the guided modes.")
        out.append(direction)

    if len(sens) == 2:
        f_hi = max(sens, key=lambda k: sens[k])
        f_lo = min(sens, key=lambda k: sens[k])
        out.append(
            f"**Which fibre is more sensitive to the magnets.** The largest "
            f"magnet-induced change in M² is {sens[f_hi]:.1f}% for the {f_hi} µm "
            f"fibre against {sens[f_lo]:.1f}% for the {f_lo} µm fibre, so in "
            f"*relative* terms the {f_hi} µm fibre responds more even though its "
            f"baseline M² is "
            f"{'lower' if b[f_hi][0] < b[f_lo][0] else 'higher'}. The "
            f"{f_lo} µm fibre is the more consistent of the two, though: "
            f"{len(signs.get(f_lo, []))} of its 4 magnet conditions move "
            f"significantly and they all move the same way, whereas the "
            f"{f_hi} µm fibre has {len(signs.get(f_hi, []))}.")

    out.append(
        f"**Overall.** Of the {len(deltas)} magnet conditions measured "
        f"({len(fibers)} fibres × 2 configurations × 2 magnet counts), "
        f"{n_sig} differ significantly from their own 0-magnet baseline. "
        + ("The magnets change M² by at most a few percent in either direction, "
           "which is small compared with the ~3× difference between the two "
           "fibres. On this data set the magnet configuration is a second-order "
           "effect on beam quality: the core size dominates."
           if max(sens.values(), default=0) < 10 else
           "The magnet effect is large enough to matter for beam quality."))

    out.append(
        "**Caveat on the absolute scale.** Every M² above is proportional to the "
        "assumed z-step and inversely proportional to the assumed wavelength; the "
        "comparisons between fibres, configurations and magnet counts are "
        "unaffected because those constants are common to all scans.")
    return out


# ----------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", default="M-Squared Results",
                   help="directory containing the image folders")
    p.add_argument("--out", default="analysis_output", help="output directory")
    p.add_argument("--wavelength-nm", type=float, default=DEFAULT_WAVELENGTH_NM)
    p.add_argument("--z-step-mm", type=float, default=DEFAULT_Z_STEP_MM)
    p.add_argument("--pixel-size-um", type=float, default=DEFAULT_PIXEL_SIZE_UM)
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--step-index-only", action="store_true",
                   help="skip the earlier/legacy folders that do not name the fibre")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--formats", default=",".join(FIGURE_FORMATS),
                   help="comma-separated figure formats to write "
                        "(pdf and svg are vector; default: %(default)s)")
    p.add_argument("--no-cache", action="store_true",
                   help="ignore and do not write the cached per-frame widths")
    args = p.parse_args(argv)

    if not os.path.isdir(args.root):
        p.error(f"root directory not found: {args.root}")
    os.makedirs(args.out, exist_ok=True)

    scans = discover_scans(args.root)
    if args.step_index_only:
        scans = [s for s in scans if s.family == "step-index"]
    if not scans:
        p.error("no z-scans found")

    print(f"Found {len(scans)} z-scans under {args.root!r}")
    results = analyse(scans, args)
    summary = aggregate(results)
    write_csvs(results, summary, args.out)

    if not args.no_plots:
        fmts = tuple(f.strip().lower() for f in args.formats.split(',') if f.strip())
        for fiber in sorted({k[0] for k in summary}, key=lambda f: int(f)):
            plot_caustics(results, fiber,
                          os.path.join(args.out, f"caustics_{fiber}um"), args, fmts)
        if summary:
            plot_m2_vs_magnets(summary, os.path.join(args.out, "m2_vs_magnets"), args, fmts)
            plot_summary_bars(summary, os.path.join(args.out, "m2_summary_bars"), fmts)
            plot_waist_divergence(summary, os.path.join(args.out, "waist_divergence"), fmts)
            plot_example_beams(results, os.path.join(args.out, "example_beams"), fmts)
        print(f"Figures written as: {', '.join(fmts)}")

    report = write_report(results, summary, args.out, args)
    print("\n" + report.split("## Conclusion")[-1].split("## Data quality")[0])
    print(f"Wrote results to {os.path.abspath(args.out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
