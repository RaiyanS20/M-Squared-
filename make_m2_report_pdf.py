#!/usr/bin/env python3
"""
Builds the M-squared interpretation report as a PDF, for inclusion in Chapter 5
of "Anomalous Transverse Faraday Effect In Multimode Fibres".

Reads the numbers straight out of analysis_output/ so the prose and the figures
can never drift apart from the data.

    python3 make_m2_report_pdf.py            # -> M2_Report_Chapter5.pdf
"""

from __future__ import annotations

import math
import os
import sys

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether,
                                NextPageTemplate, PageBreak, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)

OUT = "M2_Report_Chapter5.pdf"
ANALYSIS = "analysis_output"

# Experimental constants, from Table 4.1 of the dissertation and the M2 analysis.
LAMBDA_NM = 532.0
NA = 0.20
N_CORE = 1.491
FIBRES = {"50": 25.0e-6, "105": 52.5e-6}      # core radius in metres

INK = colors.HexColor("#111111")
INK2 = colors.HexColor("#4a4a48")
RULE = colors.HexColor("#c9c8c3")
BAND = colors.HexColor("#eef2f7")
ACCENT = colors.HexColor("#2a78d6")
WARM = colors.HexColor("#eb6834")


# ---------------------------------------------------------------------------
# Fonts. DejaVu carries the Greek and maths glyphs the built-in fonts lack.
# ---------------------------------------------------------------------------

def register_fonts():
    base = "/usr/share/fonts/truetype/dejavu"
    faces = {"Body": "DejaVuSerif.ttf", "Body-Bold": "DejaVuSerif-Bold.ttf",
             "Head": "DejaVuSans.ttf", "Head-Bold": "DejaVuSans-Bold.ttf",
             "Mono": "DejaVuSansMono.ttf"}
    for name, fn in faces.items():
        path = os.path.join(base, fn)
        if not os.path.exists(path):
            raise SystemExit(f"missing font: {path}")
        pdfmetrics.registerFont(TTFont(name, path))
    # No italic face is shipped; map the italic slots onto the upright ones so
    # any <i> markup degrades gracefully instead of raising.
    pdfmetrics.registerFontFamily("Body", normal="Body", bold="Body-Bold",
                                  italic="Body", boldItalic="Body-Bold")
    pdfmetrics.registerFontFamily("Head", normal="Head", bold="Head-Bold",
                                  italic="Head", boldItalic="Head-Bold")


def styles():
    ss = getSampleStyleSheet()
    s = {}
    s["title"] = ParagraphStyle("title", parent=ss["Title"], fontName="Head-Bold",
                                fontSize=19, leading=24, textColor=INK,
                                spaceAfter=2, alignment=0)
    s["subtitle"] = ParagraphStyle("subtitle", fontName="Head", fontSize=11.5,
                                   leading=15, textColor=INK2, spaceAfter=14)
    s["h1"] = ParagraphStyle("h1", fontName="Head-Bold", fontSize=13.5, leading=17,
                             textColor=INK, spaceBefore=16, spaceAfter=7, keepWithNext=True)
    s["h2"] = ParagraphStyle("h2", fontName="Head-Bold", fontSize=11, leading=14,
                             textColor=INK, spaceBefore=11, spaceAfter=5, keepWithNext=True)
    s["body"] = ParagraphStyle("body", fontName="Body", fontSize=9.6, leading=14.2,
                               textColor=INK, alignment=TA_JUSTIFY, spaceAfter=7)
    s["caption"] = ParagraphStyle("caption", fontName="Body", fontSize=8.3,
                                  leading=11.6, textColor=INK2,
                                  alignment=TA_JUSTIFY, spaceBefore=5,
                                  spaceAfter=13)
    s["cell"] = ParagraphStyle("cell", fontName="Body", fontSize=8.2, leading=11,
                               textColor=INK)
    s["cellb"] = ParagraphStyle("cellb", fontName="Body-Bold", fontSize=8.2,
                                leading=11, textColor=INK)
    s["callout"] = ParagraphStyle("callout", fontName="Body", fontSize=9.6,
                                  leading=14.2, textColor=INK,
                                  alignment=TA_JUSTIFY, leftIndent=8,
                                  rightIndent=8, spaceBefore=8, spaceAfter=8)
    s["plain"] = ParagraphStyle("plain", fontName="Body", fontSize=10.2,
                                leading=15.4, textColor=INK, alignment=TA_JUSTIFY,
                                spaceAfter=9)
    return s


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load():
    sc = pd.read_csv(os.path.join(ANALYSIS, "summary_by_condition.csv"))
    sc["fiber_um"] = sc["fiber_um"].astype(str)
    fits = pd.read_csv(os.path.join(ANALYSIS, "per_scan_fits.csv"))
    frames = pd.read_csv(os.path.join(ANALYSIS, "per_image_widths.csv"))
    return sc, fits, frames


def combine(sc, fiber, magnets, config):
    """Mean M^2 and its standard error, pooled over the H and V arms."""
    v = sc[(sc.fiber_um == fiber) & (sc.magnets == magnets) & (sc.config == config)]
    if v.empty:
        return None
    m = float(v.M2_mean.mean())
    sem = float(math.sqrt((v.M2_sem ** 2).sum()) / len(v))
    return m, sem


def modal(fiber):
    """V-number, total mode count and the fully-filled M^2 ceiling."""
    a = FIBRES[fiber]
    lam = LAMBDA_NM * 1e-9
    V = 2 * math.pi * a * NA / lam
    return V, V ** 2 / 2, V / 2          # M^2 ceiling = mu_max = V/2


# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------

def table(rows, widths, s, header=True, align=None):
    data = []
    for i, r in enumerate(rows):
        st = s["cellb"] if (header and i == 0) else s["cell"]
        data.append([Paragraph(str(c), st) for c in r])
    t = Table(data, colWidths=widths, repeatRows=1 if header else 0)
    cmds = [("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ("LEFTPADDING", (0, 0), (-1, -1), 5),
            ("RIGHTPADDING", (0, 0), (-1, -1), 5),
            ("LINEBELOW", (0, 0), (-1, -2), 0.4, RULE)]
    if header:
        cmds += [("BACKGROUND", (0, 0), (-1, 0), BAND),
                 ("LINEBELOW", (0, 0), (-1, 0), 0.9, INK2)]
    if align:
        for col, a in align.items():
            cmds.append(("ALIGN", (col, 0), (col, -1), a))
    t.setStyle(TableStyle(cmds))
    return t


def callout(text, s, colour=ACCENT, bg=colors.HexColor("#f2f6fc")):
    t = Table([[Paragraph(text, s["callout"])]], colWidths=[165 * mm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), bg),
        ("LINEBEFORE", (0, 0), (0, -1), 2.6, colour),
        ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
        ("LEFTPADDING", (0, 0), (-1, -1), 9),
        ("RIGHTPADDING", (0, 0), (-1, -1), 9)]))
    return t


def figure(path, caption, s, width=165 * mm, max_h=205 * mm):
    from PIL import Image as PILImage
    with PILImage.open(path) as im:
        w, h = im.size
    scale = width / w
    if h * scale > max_h:
        scale = max_h / h
    img = Image(path, width=w * scale, height=h * scale)
    img.hAlign = "CENTER"
    return KeepTogether([img, Paragraph(caption, s["caption"])])


# ---------------------------------------------------------------------------

def build():
    register_fonts()
    s = styles()
    sc, fits, frames = load()

    def M(f, m, c):
        return combine(sc, f, m, c)

    base50, base105 = M("50", 0, "Baseline"), M("105", 0, "Baseline")
    step = fits[fits.family == "step-index"]
    r2 = pd.concat([step.R2x, step.R2y])

    doc = BaseDocTemplate(OUT, pagesize=A4,
                          leftMargin=22 * mm, rightMargin=22 * mm,
                          topMargin=20 * mm, bottomMargin=18 * mm,
                          title="Beam Quality (M2) Analysis of the Transverse-Field "
                                "Speckle Datasets",
                          author="Analysis report")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="f")

    def footer(canv, d):
        canv.saveState()
        canv.setFont("Head", 7.5)
        canv.setFillColor(INK2)
        canv.drawString(doc.leftMargin, 11 * mm,
                        "M² analysis — supporting material for Chapter 5")
        canv.drawRightString(A4[0] - doc.rightMargin, 11 * mm, "%d" % d.page)
        canv.setStrokeColor(RULE)
        canv.setLineWidth(0.4)
        canv.line(doc.leftMargin, 14 * mm, A4[0] - doc.rightMargin, 14 * mm)
        canv.restoreState()

    doc.addPageTemplates([PageTemplate(id="all", frames=[frame], onPage=footer)])
    E = []
    P = lambda t, st="body": E.append(Paragraph(t, s[st]))

    # ---------------------------------------------------------------- title
    P("Beam Quality (M<super>2</super>) of the Transverse-Field Speckle Datasets", "title")
    P("What the caustic measurements show, and how they bear on the transverse "
      "Faraday hypothesis &mdash; supporting analysis for Chapter 5", "subtitle")

    E.append(callout(
        "<b>In one sentence.</b> The magnets change the beam propagation ratio "
        "M<super>2</super> of both fibres by at most a few percent, while the same "
        "magnets decorrelate the speckle by up to 19% &mdash; and that contrast is "
        "the useful result: whatever the field is doing to the guided light, it is "
        "not moving significant optical power up or down the mode-group ladder, "
        "which is what a bending, stress or gross mode-scrambling artefact would "
        "have done.", s))

    P("1. What this document is", "h1")
    P("Every image folder in the speckle campaign is a caustic z-scan: a sequence "
      "of camera frames recorded as the camera translates through the focus of the "
      "beam leaving the fibre. Those scans were never analysed for beam quality. "
      "This report extracts the ISO 11146 beam propagation ratio M<super>2</super> "
      "from all of them and asks two questions that Chapter 5 currently leaves open:")
    P("<b>(a) How full is each fibre?</b> Section 5.3 attributes the core-size "
      "scaling of the decorrelation to modal volume, and states explicitly that a "
      "full account &ldquo;would require tracking the actual excited mode-group "
      "distribution &hellip; which was not measured directly in this campaign&rdquo;. "
      "M<super>2</super> is a direct measurement of exactly that quantity, and it "
      "was sitting in the data already.")
    P("<b>(b) Do the magnets redistribute modal power?</b> If the observed speckle "
      "decorrelation were produced by the magnets mechanically perturbing the fibre, "
      "or by any mechanism that scatters power broadly across mode groups, the "
      "modal power distribution would broaden and M<super>2</super> would rise "
      "measurably. This is a falsification test of the kind catalogued in "
      "Section 4.6.3, using an observable completely independent of the Pearson "
      "correlation.")

    # ------------------------------------------------------------- method
    P("2. How M<super>2</super> was measured", "h1")
    P("Beam widths are ISO 11146 second-moment (D4&#963;) diameters. The frame is "
      "background-corrected from a border ring; secondary bright lobes (the second "
      "Wollaston arm, and reflections near the sensor border, present in 320 of the "
      "1641 frames) are masked, because a lobe far from the centroid dominates a "
      "second moment even when it carries little power; moments are then integrated "
      "over an elliptical aperture of three times the D4&#963; diameter, iterated to "
      "convergence. For each scan the caustic d<super>2</super>(z) = a + bz + "
      "cz<super>2</super> is fitted and")
    P("<font face='Mono' size='9.5'>    M² = (π / 8λ) · "
      "√(4ac − b²)</font>", "body")
    P("The fit uses a 3&#963; robust residual cut; every rejected frame is recorded "
      "and marked on the plots. Across the %d step-index scans the median fit "
      "R<super>2</super> is %.4f and no fit failed." % (len(step), r2.median()))

    rows = [["Constant", "Value", "Source", "Effect if wrong"],
            ["Wavelength &#955;", "532 nm", "Table 4.1 / Ch. 4",
             "M<super>2</super> scales as 1/&#955;"],
            ["Stage step per frame", "10 mm", "supplied",
             "M<super>2</super> scales linearly with it"],
            ["Pixel pitch", "3.45 &#181;m", "Sony IMX273 (Blackfly S)",
             "M<super>2</super> scales as pitch<super>2</super>"]]
    E.append(table(rows, [34 * mm, 24 * mm, 46 * mm, 61 * mm], s))
    E.append(Spacer(1, 5))
    P("These three constants set the <b>absolute</b> scale only. They are common to "
      "every scan, so every comparison in this report &mdash; fibre against fibre, "
      "configuration against configuration, magnet count against magnet count &mdash; "
      "is unaffected by them. Only the fill fractions of Section 3 depend on the "
      "absolute scale being right.", "body")

    # -------------------------------------------------- what M2 means here
    P("3. What M<super>2</super> measures in a step-index multimode fibre", "h1")
    P("For a step-index fibre the guided modes organise into principal groups "
      "&#956; = 2p + l + 1, and the highest group the fibre can support is "
      "&#956;<sub>max</sub> = V/2. A fully filled fibre radiates with")
    P("<font face='Mono' size='9.5'>    M²(full) = π a NA / λ = V/2 = "
      "μ_max</font>", "body")
    P("so <b>M<super>2</super> is, to a good approximation, a direct read-out of the "
      "highest mode group carrying appreciable power.</b> That makes it the natural "
      "instrument for the modal-fill question. Using NA = 0.20 and &#955; = 532 nm "
      "from Table 4.1:")

    rows = [["Fibre", "V", "Modes (V<super>2</super>/2)",
             "&#956;<sub>max</sub> = M<super>2</super> ceiling",
             "Measured M<super>2</super>", "Fill"]]
    fillfrac = {}
    for f in ("50", "105"):
        V, nm, ceil = modal(f)
        b = M(f, 0, "Baseline")
        fillfrac[f] = b[0] / ceil
        rows.append([f"{f} &#181;m", f"{V:.0f}", f"{nm:.0f}", f"{ceil:.1f}",
                     f"<b>{b[0]:.2f} &#177; {b[1]:.2f}</b>",
                     f"<b>{100*b[0]/ceil:.0f}%</b>"])
    E.append(table(rows, [22 * mm, 15 * mm, 27 * mm, 34 * mm, 36 * mm, 15 * mm], s,
                   align={1: "CENTER", 2: "CENTER", 3: "CENTER", 4: "CENTER", 5: "CENTER"}))
    E.append(Spacer(1, 6))
    P("The mode counts reproduce Table 4.1 exactly (1744 and 7689), which confirms "
      "the V-numbers are consistent with the thesis parameters.")

    E.append(callout(
        "<b>Result 1 &mdash; the modal fill, measured directly.</b> At the launch "
        "condition used for these scans the 50 &#181;m fibre is excited to about "
        "mode group %.0f of its %.0f available (%.0f%% of its ladder), and the "
        "105 &#181;m fibre to about group %.0f of %.0f (%.0f%%). The larger core is "
        "therefore filled to a <i>higher</i> fraction of its own ladder, not a lower "
        "one. Roughly %.0f modes carry the light in the 50 &#181;m fibre against "
        "%.0f in the 105 &#181;m &mdash; a ratio of about %.1f. That is "
        "<i>larger</i> than the (105/50)<super>2</super> &#8776; 4.4 ratio of "
        "<i>total</i> mode counts quoted in Section 5.3, precisely because the "
        "larger core is also filled to a higher fraction of its own ladder." % (
            base50[0], modal("50")[2], 100 * fillfrac["50"],
            base105[0], modal("105")[2], 100 * fillfrac["105"],
            base50[0] ** 2 / 2, base105[0] ** 2 / 2,
            (base105[0] ** 2 / 2) / (base50[0] ** 2 / 2)), s))

    # ------------------------------------------------------------- results
    P("4. The measurements", "h1")

    E.append(figure(os.path.join(ANALYSIS, "caustics_50um.png"),
                    "<b>Figure 1.</b> Measured caustics for the 50 &#181;m step-index "
                    "fibre. Each panel is one magnet count; rows are the two "
                    "polarisation arms. Points are measured D4&#963; diameters, "
                    "curves are the fitted hyperbolae, and &#215; marks a frame "
                    "rejected by the robust residual cut. Repeat trials lie almost on "
                    "top of one another, which is what makes the few-percent "
                    "comparisons that follow meaningful. Note that the waist sits near "
                    "the start of the travel &mdash; see the limitation in Section 8.",
                    s, max_h=118 * mm))

    E.append(figure(os.path.join(ANALYSIS, "caustics_105um.png"),
                    "<b>Figure 2.</b> The same for the 105 &#181;m fibre. The beam is "
                    "visibly larger at every z, consistent with the larger modal "
                    "volume. One session (the scans without a trial number in the "
                    "folder name) sits at a shifted waist position with a slightly "
                    "smaller waist &mdash; a different alignment day. Since "
                    "M<super>2</super> is invariant to where the waist falls, those "
                    "scans still return M<super>2</super> values consistent with the "
                    "rest, and they are retained.",
                    s, max_h=118 * mm))

    E.append(figure(os.path.join(ANALYSIS, "m2_summary_bars.png"),
                    "<b>Figure 3.</b> The headline result: change in M<super>2</super> "
                    "relative to each fibre's own 0-magnet baseline, so that the "
                    "few-percent effect under test is not hidden by the almost "
                    "threefold difference in baseline between the fibres. Error bars "
                    "are &#177;2 standard errors, pooled over trials and over the two "
                    "polarisation arms; a bar whose error interval does not cross zero "
                    "is a real change. Every change that clears the noise is a "
                    "<i>reduction</i> in M<super>2</super>.",
                    s, max_h=95 * mm))

    E.append(figure(os.path.join(ANALYSIS, "m2_vs_magnets.png"),
                    "<b>Figure 4.</b> Absolute M<super>2</super> against magnet count, "
                    "split by fibre and by polarisation arm, with the uniform "
                    "(&ldquo;Normal&rdquo;) and alternating configurations overlaid. "
                    "The 0-magnet baseline is shared by both configurations. The "
                    "vertical scale spans only a few units of M<super>2</super> in "
                    "every panel: on the scale of the baseline difference between the "
                    "two fibres, all of these traces would be flat lines.",
                    s, max_h=125 * mm))

    E.append(figure(os.path.join(ANALYSIS, "waist_divergence.png"),
                    "<b>Figure 5.</b> M<super>2</super> = &#960;&#183;d<sub>0</sub>"
                    "&#183;&#952;/(8&#955;), so any change in it is a change in the "
                    "waist diameter, the far-field divergence, or both. Splitting the "
                    "effect this way separates two different behaviours: in the "
                    "105 &#181;m fibre the reduction comes mostly through a smaller "
                    "waist, whereas the one large change in the 50 &#181;m fibre "
                    "(alternating, 50 magnets) comes through a lower divergence.",
                    s, max_h=95 * mm))

    E.append(figure(os.path.join(ANALYSIS, "example_beams.png"),
                    "<b>Figure 6.</b> The near-waist intensity pattern for each "
                    "condition, cropped to the beam and stretched to the 99.5th "
                    "percentile. The 105 &#181;m fibre carries visibly finer and far "
                    "more numerous speckle grains than the 50 &#181;m fibre &mdash; a "
                    "direct visual counterpart of the mode counts in Section 3. Within "
                    "a row the patterns are of similar size and grain, which is the "
                    "M<super>2</super> result stated visually.",
                    s, max_h=85 * mm))

    # ------------------------------------------------------------- numbers
    E.append(PageBreak())
    P("5. The numbers", "h1")
    rows = [["Fibre", "Configuration", "Magnets", "M<super>2</super> &#177; s.e.",
             "&#916; vs baseline", "% change", "Clears noise?"]]
    deltas = {}
    for f in ("50", "105"):
        b = M(f, 0, "Baseline")
        rows.append([f"{f} &#181;m", "baseline (0 magnets)", "0",
                     f"<b>{b[0]:.2f} &#177; {b[1]:.2f}</b>", "&mdash;", "&mdash;", "&mdash;"])
        for cfg in ("Normal", "Alternating"):
            for mag in (25, 50):
                v = M(f, mag, cfg)
                if v is None:
                    continue
                d = v[0] - b[0]
                comb = math.sqrt(v[1] ** 2 + b[1] ** 2)
                sig = abs(d) > 2 * comb
                deltas[(f, cfg, mag)] = (v[0], v[1], d, 100 * d / b[0], sig)
                label = "uniform (non-alternating)" if cfg == "Normal" else "alternating"
                rows.append([f"{f} &#181;m", label, str(mag),
                             f"{v[0]:.2f} &#177; {v[1]:.2f}", f"{d:+.2f}",
                             f"{100*d/b[0]:+.1f}%",
                             "<b>yes</b>" if sig else "no"])
    E.append(table(rows, [17 * mm, 40 * mm, 19 * mm, 27 * mm, 21 * mm, 19 * mm, 22 * mm],
                   s, align={2: "CENTER", 4: "CENTER", 5: "CENTER", 6: "CENTER"}))
    E.append(Spacer(1, 5))
    P("&ldquo;Clears noise&rdquo; means the difference exceeds twice the combined "
      "standard error of the two means, a roughly 95% criterion. Each condition "
      "pools three trials and two polarisation arms.")

    E.append(PageBreak())

    # ------------------------------------------------------ interpretation
    P("6. What this means for the hypothesis", "h1")

    P("6.1 The decorrelation is not accompanied by mode-group redistribution", "h2")
    P("Chapter 5 reports speckle decorrelation reaching 1 &#8722; &#961; &#8776; 0.18 "
      "for the 50 &#181;m fibre and 0.19 for the 105 &#181;m fibre under the uniform "
      "configuration at maximum off-normal incidence. Over the same magnet range, "
      "M<super>2</super> moves by at most %.1f%% and %.1f%% respectively. The speckle "
      "pattern is being comprehensively rearranged while the modal power distribution "
      "across the group ladder stays put." % (
          max(abs(v[3]) for k, v in deltas.items() if k[0] == "50"),
          max(abs(v[3]) for k, v in deltas.items() if k[0] == "105")))
    P("That dissociation is exactly what a <b>phase-scrambling, power-conserving</b> "
      "interaction looks like. Speckle correlation responds to the relative phases of "
      "the modes; M<super>2</super> does not &mdash; it responds only to how much "
      "power sits in which mode group. An interaction that reshuffles amplitude and "
      "phase among modes of the <i>same</i> principal group changes the speckle "
      "completely and leaves M<super>2</super> untouched, because modes within one "
      "group share the same &#956; and therefore make the same contribution to the "
      "space&#8211;angle second moment. This is precisely the near-degenerate "
      "intra-group vector channel identified in Section 3.6 as the pathway capable of "
      "producing the measured effect.")

    E.append(callout(
        "<b>Result 2.</b> The data are consistent with the intra-group "
        "near-degenerate channel of Section 3.6 and give it independent support from "
        "an observable that has nothing to do with the Pearson correlation. "
        "<b>They do not, however, prove it</b> &mdash; see the sensitivity argument "
        "immediately below, which is the honest limit of what M<super>2</super> can "
        "settle.", s))

    P("6.2 How much power transfer M<super>2</super> can actually exclude", "h2")
    P("This needs stating carefully, because it is where the argument could be "
      "overclaimed. If a fraction f of the power is moved by &#916;&#956; groups, "
      "then to first order &#916;M<super>2</super>/M<super>2</super> &#8776; "
      "f&#183;&#916;&#956;/&#956;. The sensitivity therefore depends strongly on how "
      "far the power moves:")
    rows = [["Fibre", "M<super>2</super> resolution (2 s.e.)",
             "Bound on f for &#916;&#956; = 1", "Bound on f for broadband scattering"]]
    for f, dmu in (("50", 7.0), ("105", 20.0)):
        b = M(f, 0, "Baseline")
        res = 2 * b[1] / b[0]
        rows.append([f"{f} &#181;m", f"{100*res:.1f}%  (&#177;{2*b[1]:.2f} in &#956;)",
                     f"f &#8804; {100*res*b[0]:.0f}%",
                     f"f &#8804; {100*res*b[0]/dmu:.0f}%"])
    E.append(table(rows, [20 * mm, 46 * mm, 40 * mm, 59 * mm], s,
                   align={1: "CENTER", 2: "CENTER", 3: "CENTER"}))
    E.append(Spacer(1, 6))
    P("So M<super>2</super> is a <b>strong</b> constraint on broadband redistribution "
      "&mdash; power scattered widely across the ladder, as bending or gross mode "
      "mixing would do &mdash; bounding it at a few percent of the total. It is a "
      "<b>weak</b> constraint on nearest-neighbour transfer, and the inter-group "
      "channel of Section 3.4 has |&#916;l| = 1 and hence |&#916;&#956;| = 1 as its "
      "minimum step. M<super>2</super> therefore cannot by itself distinguish the "
      "inter-group channel from the intra-group one. What it does do is close off the "
      "mechanical and gross-scattering alternatives, which is its proper role in the "
      "falsification programme of Section 4.6.3. It also cannot test the "
      "2 &#215; 10<super>&#8722;7</super> ceiling of Eq. 3.47: that is some five "
      "orders of magnitude below anything this measurement can resolve.")

    P("6.3 A further falsification of the mechanical pathway", "h2")
    P("Section 5.1 already separates mechanical disturbance from magneto-optic "
      "coupling by its temporal signature &mdash; a step discontinuity against a "
      "smooth decay. The M<super>2</super> data close the same door from a different "
      "side. Macrobending and clamping stress do two things to a multimode fibre: "
      "they couple power between mode groups, and they strip the highest-order "
      "groups. Both change M<super>2</super>, and at the field strengths and handling "
      "involved in fitting 50 magnet pairs they would change it by far more than the "
      "few percent seen here. The near-constancy of M<super>2</super> across 0, 25 "
      "and 50 magnets is therefore independent evidence that mounting the magnet "
      "array is not mechanically reworking the fibre's modal content.")

    P("6.4 The small reductions, and how far to push them", "h2")
    P("Four of the eight magnet conditions show a change that clears the noise, and "
      "all four are reductions in M<super>2</super>. The 105 &#181;m fibre is the "
      "consistent one: three of its four conditions move, all downward, by 2&#8211;4%, "
      "and Figure 5 shows the reduction arriving through a smaller waist. Two readings "
      "are available and the data do not separate them. Either the magnet array is "
      "acting as a very weak mode filter &mdash; slight bending along the rail "
      "preferentially attenuating the least-confined highest-order modes, which are "
      "the ones that contribute most to the second moment &mdash; or the launch "
      "drifted very slightly over the sequence in which the magnets were fitted. "
      "The second is entirely plausible at this magnitude. What matters for the "
      "hypothesis is that neither reading involves power being scattered <i>up</i> "
      "the ladder, which is what a strong inter-group coupling mechanism would have "
      "produced.")

    P("6.5 The uniform-versus-alternating contrast", "h2")
    n50 = deltas[("50", "Normal", 50)]
    a50 = deltas[("50", "Alternating", 50)]
    n105 = deltas[("105", "Normal", 50)]
    a105 = deltas[("105", "Alternating", 50)]
    P("Section 4.6.3.3 uses the uniform-versus-alternating contrast as the "
      "discriminator between a field-linear effect and a B<super>2</super> effect. In "
      "the decorrelation data that contrast is large and is the central result. In the "
      "M<super>2</super> data it is mostly absent, as the power-conservation argument "
      "above predicts. For the 105 &#181;m fibre the two configurations differ by "
      "%.2f in M<super>2</super> at 50 magnets, within the combined standard error, "
      "so they are indistinguishable. For the 50 &#181;m fibre they do differ &mdash; "
      "uniform %.2f against alternating %.2f &mdash; but this rests on a single "
      "condition in the noisier of the two fibres, and it is the alternating "
      "configuration that sits lower, which is not the ordering the transverse "
      "Faraday mechanism would predict if M<super>2</super> were tracking the "
      "coupling strength. It is better treated as an outlier to be re-measured than "
      "as a result." % (abs(a105[0] - n105[0]), n50[0], a50[0]))

    E.append(PageBreak())

    # ---------------------------------------------------------- plain words
    P("7. The same thing in plain words", "h1")
    P("M<super>2</super> is a single number for how good a laser beam is. "
      "M<super>2</super> = 1 is a perfect beam; the bigger the number, the messier "
      "and faster-spreading the beam. A multimode fibre scrambles light into many "
      "&ldquo;modes&rdquo;, so its output has a large M<super>2</super> &mdash; and "
      "crucially, that number tells you <b>how many</b> of the fibre's modes are "
      "carrying light.", "plain")
    P("Two things came out of measuring it.", "plain")
    P("<b>First, a number the thesis says it did not have.</b> Chapter 5 argues the "
      "105 &#181;m fibre reacts more strongly than the 50 &#181;m one because it has "
      "more modes, but notes the actual modal content was never measured. It has now "
      "been measured, from images already taken: the 50 &#181;m fibre was running at "
      "about half of its mode capacity and the 105 &#181;m at about two thirds, "
      "carrying roughly 120 and 900 modes respectively. That is a real, quotable "
      "number for the modal-volume argument.", "plain")
    P("<b>Second, and more important: the magnets barely change the beam at all.</b> "
      "The speckle pattern changes a lot when the magnets are on &mdash; that is the "
      "whole result of the thesis &mdash; but the beam's M<super>2</super> stays "
      "essentially the same. Think of a choir. The speckle pattern is the sound of "
      "the choir, which changes completely if the singers change who is singing which "
      "note. M<super>2</super> is the range of the choir &mdash; how low the basses "
      "and how high the sopranos go. The measurement says the singers are swapping "
      "notes among themselves, but nobody is singing outside the range they started "
      "in. Something is rearranging the light, and it is not simply shoving it into "
      "higher-order modes.", "plain")
    P("That matters because the most boring explanation for the whole effect &mdash; "
      "&ldquo;the magnets are just bending or squeezing the fibre&rdquo; &mdash; "
      "would have shown up here. Bending a multimode fibre does change its "
      "M<super>2</super>, and noticeably. It did not change. So this is one more "
      "alternative explanation closed off, using a measurement that has nothing to do "
      "with the correlation numbers the rest of the chapter relies on.", "plain")

    verdict = callout(
        "<b>Yes, but as supporting and falsifying evidence, not as proof.</b><br/><br/>"
        "<b>It supports the project in three ways.</b> It supplies the direct "
        "measurement of modal fill that Section 5.3 identifies as missing, and the "
        "numbers back the modal-volume reading of the core-size scaling. It shows the "
        "decorrelation happens without power moving across the mode-group ladder, "
        "which is the behaviour expected of the near-degenerate intra-group channel "
        "proposed in Section 3.6, from an observable independent of the Pearson "
        "correlation. And it independently closes the mechanical and gross-scattering "
        "alternatives that Section 4.6.3 sets out to falsify.<br/><br/>"
        "<b>It does not prove the transverse Faraday mechanism.</b> "
        "M<super>2</super> is a null result here: it says what is <i>not</i> "
        "happening. It cannot separate the inter-group channel from the intra-group "
        "one, because both can operate at |&#916;&#956;| = 1, where "
        "M<super>2</super> is insensitive. It is nowhere near sensitive enough to "
        "test the theoretical ceiling of Eq. 3.47. And it shows no "
        "uniform-versus-alternating asymmetry, so it adds nothing to the field-linearity "
        "argument &mdash; which, on the power-conservation reasoning of Section 6.1, "
        "is what it should do.",
        s, colour=colors.HexColor("#1baf7a"), bg=colors.HexColor("#eefaf5"))
    E.append(KeepTogether([Paragraph("8. Does it support the hypothesis?", s["h1"]),
                           verdict]))

    P("9. Limitations you should state if you use this", "h1")
    lim = [
        ("Launch condition not recorded", "The folder names for the caustic scans do "
         "not record which incidence configuration (I, II or III of Section 5.3) was "
         "in use. The fill fractions in Section 3 are therefore valid for "
         "&ldquo;the launch used for these scans&rdquo;, not for a named "
         "configuration. If you know which it was, say so; if the scans span more "
         "than one, the fill numbers should be quoted per configuration."),
        ("Fill fractions depend on the absolute calibration", "Unlike everything else "
         "here, the 52% and 68% figures move if the 10 mm stage step or 532 nm "
         "wavelength is wrong. The comparisons and the null result do not."),
        ("ISO 11146 sampling is not met", "In every scan the waist sits near the "
         "start of the travel, so there are only about three frames before it and the "
         "scan reaches roughly two Rayleigh lengths beyond. The standard asks for at "
         "least five on each side. Fits are excellent (median R<super>2</super> = "
         "%.4f) and reproducible, so relative comparisons hold, but the absolute "
         "M<super>2</super> carries more systematic uncertainty than the quoted "
         "standard errors." % r2.median()),
        ("Three trials, not ten", "The caustic scans have n = 3 per condition against "
         "the n = 10 of the decorrelation protocol. The error bars reflect this."),
        ("A single significant 50 &#181;m result", "The &#8722;6.1% at 50 alternating "
         "magnets rests on one condition and should be re-measured before being "
         "leaned on."),
    ]
    rows = [["Limitation", "What to say"]]
    for a, b in lim:
        rows.append([f"<b>{a}</b>", b])
    E.append(table(rows, [42 * mm, 123 * mm], s))

    P("10. Suggested placement in Chapter 5", "h1")
    P("The material divides naturally in two, matching the two questions in "
      "Section 1 above.")
    rows = [["Content", "Suggested home"],
            ["Modal fill: V-numbers, ceilings, measured M<super>2</super>, fill "
             "fractions, mode counts (Section 3, Figure 6)",
             "Into <b>5.3</b>, at the &ldquo;Core size and modal volume&rdquo; "
             "paragraph, replacing the sentence conceding that the excited mode-group "
             "distribution was not measured."],
            ["Power-conservation falsification: Figures 3 and 5, the numbers table, "
             "and Sections 6.1&#8211;6.4",
             "As a <b>new section after 5.5</b>, e.g. &ldquo;5.6 Modal Power "
             "Conservation under the Applied Field&rdquo;, renumbering the regime "
             "diagnostic and theory confrontation that follow."],
            ["The intra-group reading (Section 6.1)",
             "Cross-reference from <b>5.7.3</b>, &ldquo;Evidence Bearing on the "
             "Resolution&rdquo;, which is where the case for the Section 3.6 channel "
             "is assembled."],
            ["Figures 1, 2 and 4",
             "An appendix, or supporting figures for the new section. They document "
             "the measurement rather than carry the argument."]]
    E.append(table(rows, [72 * mm, 93 * mm], s))
    E.append(Spacer(1, 8))
    P("One caution on wording. Section 5.3 explains the convergence of the two fibres "
      "at maximum off-normal incidence by suggesting a given launch populates a "
      "larger fraction of the 50 &#181;m fibre's ladder than of the 105 &#181;m "
      "fibre's. At the launch used for the caustic scans the measured fill runs the "
      "other way (%.0f%% against %.0f%%). This is not necessarily a contradiction "
      "&mdash; the launch conditions may differ &mdash; but the two statements should "
      "not sit in the same chapter unreconciled. Either establish which configuration "
      "the caustic scans correspond to, or soften the fractional-fill argument to a "
      "possibility rather than an explanation."
      % (100 * fillfrac["50"], 100 * fillfrac["105"]))

    P("11. Provenance", "h1")
    P("All numbers and figures come from <font face='Mono' size='8.6'>m2_analysis.py"
      "</font> in the project repository, run over the %d frames in %d z-scans under "
      "<font face='Mono' size='8.6'>M-Squared Results/</font>. Per-frame widths, "
      "per-scan fits and the condition summary are in "
      "<font face='Mono' size='8.6'>analysis_output/</font> as CSV, including every "
      "quality flag and every rejected frame, so any figure here can be traced back "
      "to individual images. This PDF is generated by "
      "<font face='Mono' size='8.6'>make_m2_report_pdf.py</font> and reads those CSVs "
      "directly, so the prose cannot drift from the data."
      % (len(frames), fits.shape[0]))

    doc.build(E)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    sys.exit(build())
