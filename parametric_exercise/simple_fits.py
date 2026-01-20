import ROOT
from config import plot_dir

def pdf_to_hist(
    pdf,
    obs,
    nbins,
    xmin,
    xmax,
    total_yield,
    hist_name="h_exp"
):
    """
    Convert a RooAbsPdf into a binned TH1 scaled to a given total expected yield.

    Parameters
    ----------
    pdf : ROOT.RooAbsPdf
        PDF stored in the workspace
    obs : ROOT.RooRealVar
        Observable
    nbins : int
        Number of bins
    xmin, xmax : float
        Histogram range
    total_yield : float
        Total expected yield in [xmin, xmax]
    hist_name : str
        Name of output histogram

    Returns
    -------
    ROOT.TH1
    """

    # # Define range (important!)
    # obs.setRange("histRange", xmin, xmax)

    # Create histogram from the PDF
    h = pdf.createHistogram(
        hist_name,
        obs,
        ROOT.RooFit.Binning(nbins, xmin, xmax),
        # ROOT.RooFit.Range("histRange")
    )

    # h.SetDirectory(0)

    # Scale histogram to the expected yield
    # integral = h.Integral(1, nbins)
    integral = h.Integral()
    if integral > 0:
        h.Scale(total_yield / integral)

    return h


ROOT.gROOT.SetBatch(True)
print("Plotting directory: %s" % plot_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Open the best-fit file and load the postfit workspace
f = ROOT.TFile("higgsCombine.bestfit.MultiDimFit.mH125.root")
w = f.Get("w")
w.Print("v")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lets plot the post-fit and prefit model to the data
n_bins = 80
binning = ROOT.RooFit.Binning(n_bins,100,180)

can = ROOT.TCanvas()
plot = w.var("CMS_hgg_mass").frame()
w.data("data_obs").plotOn( plot, binning )
print(f"data obs sum: {w.data("data_obs").sumEntries()}")
# Load the S+B model
sb_model = w.pdf("model_s").getPdf("Tag0")

# Prefit
sb_model.plotOn(plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("prefit") )

# Postfit
w.loadSnapshot("MultiDimFit")
sb_model.plotOn(plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("postfit") )
r_bestfit = w.var("r").getVal()

plot.Draw()

leg = ROOT.TLegend(0.55,0.6,0.85,0.85)
leg.AddEntry("prefit", "Prefit S+B model (r=1.00)", "L")
leg.AddEntry("postfit", "Postfit S+B model (r=%.2f)"%r_bestfit, "L")
leg.Draw("Same")

can.Update()
can.SaveAs("%s/part2_sb_model.png"%plot_dir)

b_model = w.pdf("shapeBkg_bkg_mass_Tag0")
print("---"*50)
# b_model.Print()
# sb_model.Print("v")
b_model_norm = w.obj(b_model.GetName()+"__norm")
print(f"b_model_norm: {b_model_norm.getVal()}")
s_model_norm = w.obj("shapeSig_ggH_Tag0__norm")
print(f"s_model_norm: {s_model_norm.getVal()}")

# xmin=100
# xmax=180
# n_bins = 80
# h_bkg = pdf_to_hist(
#     pdf=bkg_pdf,
#     obs=obs,
#     nbins=n_bins,
#     xmin=xmin,
#     xmax=xmax,
#     total_yield=nbkg,
#     hist_name="h_bkg"
# )
