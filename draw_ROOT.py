import ROOT
import os
from pathlib import Path

ROOT.gROOT.SetBatch(True)  # 🔥 kluczowe: bez GUI

# --- output directory setup ---
SCRIPT_DIR = Path(__file__).resolve().parent
OUTDIR = SCRIPT_DIR / "results"
OUTDIR.mkdir(parents=True, exist_ok=True)


def save_single(hist, outname, title=None):
    c = ROOT.TCanvas("c", "c", 800, 600)
    hist.SetLineWidth(2)
    if title:
        hist.SetTitle(title)
    hist.Draw("hist")
    c.SaveAs(str(OUTDIR / outname))
    c.Close()


def save_overlay(hists, labels, colors, outname, title=None):
    c = ROOT.TCanvas("c", "c", 800, 600)

    leg = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)

    for i, h in enumerate(hists):
        h.SetLineColor(colors[i])
        h.SetLineWidth(2)
        if i == 0:
            if title:
                h.SetTitle(title)
            h.Draw("hist")
        else:
            h.Draw("hist same")
        leg.AddEntry(h, labels[i], "l")

    leg.Draw()
    c.SaveAs(str(OUTDIR / outname))
    c.Close()


def main():
    f = ROOT.TFile.Open("kinfit_3pr_tt.root")

    # --- get histograms ---
    mvis = f.Get("mvis")
    mtt = f.Get("mtt")

    dpt1_kinfit = f.Get("dpt1_kinfit")
    dpt1_fast = f.Get("dpt1_FastMTT")
    dpt1_const = f.Get("dpt1_FastMTT_const")

    dpt2_kinfit = f.Get("dpt2_kinfit")
    dpt2_fast = f.Get("dpt2_FastMTT")
    dpt2_const = f.Get("dpt2_FastMTT_const")

    # --- single plots ---
    save_single(mvis, "mvis.png", "Visible mass")
    save_single(mtt, "mtt.png", "KinFit mass")

    # --- overlays ---
    save_overlay(
        [mvis, mtt],
        ["mvis", "kinfit"],
        [ROOT.kBlack, ROOT.kRed],
        "mass_comparison.png",
        "Mass comparison",
    )

    save_overlay(
        [dpt1_kinfit, dpt1_fast, dpt1_const],
        ["kinfit", "FastMTT", "FastMTT const"],
        [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2],
        "dpt1_comparison.png",
        "Tau1 pT resolution",
    )

    save_overlay(
        [dpt2_kinfit, dpt2_fast, dpt2_const],
        ["kinfit", "FastMTT", "FastMTT const"],
        [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2],
        "dpt2_comparison.png",
        "Tau2 pT resolution",
    )

    print(f"✅ PNG files saved in: {OUTDIR}")


if __name__ == "__main__":
    main()