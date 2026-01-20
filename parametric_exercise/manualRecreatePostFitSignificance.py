import ROOT
import math


ROOT.gROOT.SetBatch(True)
def poisson_nll(h_data, h_exp):
    """
    Compute the Poisson log-likelihood for binned data.

    Parameters
    ----------
    h_data : ROOT.TH1
        Observed data histogram
    h_exp : ROOT.TH1
        Expected histogram (signal + background)

    Returns
    -------
    logL : float
        Log-likelihood value
    """

    if h_data.GetNbinsX() != h_exp.GetNbinsX():
        raise ValueError("Histograms must have the same number of bins")

    logL = 0.0

    for i in range(1, h_data.GetNbinsX() + 1):
        n = h_data.GetBinContent(i)
        lam = h_exp.GetBinContent(i)

        if lam <= 0:
            continue  # avoid log(0)

        logL += n * math.log(lam) - lam - math.lgamma(n + 1)

    return logL



def plot_two_hists(
    h1,
    h2,
    out_png="two_histograms.png",
    title="Two histograms",
    x_title="x",
    y_title="Events",
    label1="Histogram 1",
    label2="Histogram 2"
):
    """
    Plot two TH1 histograms on the same canvas and save to a PNG.

    Parameters
    ----------
    h1, h2 : ROOT.TH1
        Histograms to be plotted
    out_png : str
        Output PNG filename
    title : str
        Canvas title
    x_title, y_title : str
        Axis labels
    label1, label2 : str
        Legend labels
    """

    ROOT.gROOT.SetBatch(True)

    # Style
    h1.SetLineColor(ROOT.kBlue)
    h1.SetLineWidth(2)

    h2.SetLineColor(ROOT.kRed)
    h2.SetLineWidth(2)

    # Canvas
    c = ROOT.TCanvas("c", title, 800, 600)

    # Axis titles
    h1.SetTitle(f"{title};{x_title};{y_title}")

    # Draw
    h1.Draw("HIST")
    h2.Draw("HIST SAME")

    # Legend
    leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
    leg.AddEntry(h1, label1, "l")
    leg.AddEntry(h2, label2, "l")
    leg.Draw()

    # Save
    c.SaveAs(out_png)

    return c


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


def plot_pdf_vs_hist(
    pdf,
    obs,
    th1,
    xmin,
    xmax,
    data,
    fit_range="loSB,hiSB",
    out_png="pdf_vs_hist.png",
):
    """
    Overlay a RooAbsPdf and a binned expected histogram and save to a PNG.
    """
    # Define range (important!)
    # plot_range = "histRange"
    # obs.setRange(plot_range, xmin, xmax)
    
    # RooPlot frame
    frame = obs.frame(ROOT.RooFit.Range(xmin, xmax))
    


    binning = ROOT.RooFit.Binning(n_bins,xmin,xmax)
    data.plotOn(frame, binning, Invisible=True)
    # data.plotOn(frame, binning)
    # pdf.plotOn(frame, ROOT.RooFit.NormRange(fit_range), ROOT.RooFit.Range(plot_range)) 
    pdf.plotOn( frame, ROOT.RooFit.NormRange(fit_range), ROOT.RooFit.Range("full"),
                ROOT.RooFit.LineColor(2), Name=pdf.GetName()
    )
    
    # Canvas
    c = ROOT.TCanvas("c", "PDF vs Histogram", 800, 600)

    frame.Draw()

    # Histogram styling
    th1.SetLineColor(ROOT.kBlack)
    th1.SetLineWidth(2)
    th1.SetMarkerStyle(20)
    th1.SetMarkerSize(0.9)

    th1.Draw("HIST SAME")

    # Legend
    leg = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg.AddEntry(th1, "Binned expectation", "lep")
    leg.AddEntry(frame.findObject(pdf.GetName()), "PDF (scaled)", "l")
    leg.Draw()

    c.Update()
    c.SaveAs(out_png)

    return c

def roodataset_to_th1(
    data,
    obs,
    nbins,
    xmin,
    xmax,
    hist_name="h_data"
):
    """
    Convert an unbinned RooDataSet into a TH1F.
    """

    # Create histogram
    h = ROOT.TH1F(hist_name, "", nbins, xmin, xmax)
    h.Sumw2()

    # Loop over entries
    for i in range(data.numEntries()):
        row = data.get(i)
        x = row.getRealValue(obs.GetName())
        w = data.weight() if data.isWeighted() else 1.0
        h.Fill(x, w)

    return h

def fq0(q0):
    """
    Source: equation 35 from https://inspirehep.net/literature/1723231
    Asymptotic PDF f(q0 | 0) for the test statistic q0.

    Parameters
    ----------
    q0 : float
        Test statistic value (q0 >= 0)

    Returns
    -------
    float
        PDF value (delta component handled explicitly)
    """

    if q0 < 0:
        return 0.0

    # Dirac delta contribution at q0 = 0
    if q0 == 0:
        return float("inf")  # formal representation of delta(q0)

    # Continuous part for q0 > 0
    return (
        0.5
        * (1.0 / math.sqrt(2.0 * math.pi))
        * (1.0 / math.sqrt(q0))
        * math.exp(-0.5 * q0)
    )

def scale_roodataset_weights(data, scale, name="data_scaled"):
    """
    Return a new RooDataSet with all weights scaled by `scale`.
    """

    obs = data.get()  # RooArgSet of observables

    data_scaled = ROOT.RooDataSet(
        name,
        name,
        obs,
        ROOT.RooFit.WeightVar(data.weightVar().GetName())
    )

    for i in range(data.numEntries()):
        obs_vals = data.get(i)
        w = data.weight()
        data_scaled.add(obs_vals, scale * w)

    return data_scaled

if __name__ == "__main__":
    root_file = "workspace_bkg.root"
    f = ROOT.TFile.Open(root_file)
    f.ls("v")
    ws_name= "workspace_bkg"
    ws = f.Get(ws_name)
    # ws.Print()
    bkg_pdf_name="model_bkg_Tag0"
    bkg_pdf_norm_name= bkg_pdf_name + "_norm"
    bkg_pdf = ws.pdf(bkg_pdf_name)
    bkg_pdf.Print()
    nbkg = ws.obj(bkg_pdf_norm_name).getVal()
    
    obs = ws.var("CMS_hgg_mass")
    xmin=100
    xmax=180
    n_bins = 80
    h_bkg = pdf_to_hist(
        pdf=bkg_pdf,
        obs=obs,
        nbins=n_bins,
        xmin=xmin,
        xmax=xmax,
        total_yield=nbkg,
        hist_name="h_bkg"
    )
    data_name = "data_Tag0"
    data = ws.obj(data_name)

    # ------------------------------------------------------
    # Sanity check: plot background RooAbsPdf and histogram
    # that we got out of it
    # ------------------------------------------------------
    plot_pdf_vs_hist(
        pdf=bkg_pdf,
        obs=obs,
        th1=h_bkg,
        xmin=xmin,
        xmax=xmax,
        data=data,
        out_png="bkg_pdf_sanity.png"
    )
    """
    NOTE: the h_bkg has slightly higher SumEntries compared to the fit function
    this is expected since nbkg also accounts for a "bump" in the H-peak region
    that the fit function doesn't know bc it fit only over the H-sidebands region
    """

    # same on signal model
    root_file = "workspace_sig_with_norm.root"
    f_sig = ROOT.TFile.Open(root_file)
    f_sig.ls("v")
    ws_name= "workspace_sig"
    ws_sig = f_sig.Get(ws_name)
    # ws.Print()
    sig_pdf_name="model_ggH_Tag0"
    sig_pdf_norm_name= sig_pdf_name + "_norm"
    sig_pdf = ws_sig.pdf(sig_pdf_name)
    sig_pdf.Print()
    Run2Lumi = 138000 # this information is added to "rate" section of the datacard
    nsig = ws_sig.obj(sig_pdf_norm_name).getVal()*Run2Lumi
    
    # nsig
    h_sig = pdf_to_hist(
        pdf=sig_pdf,
        obs=obs,
        nbins=n_bins,
        xmin=xmin,
        xmax=xmax,
        total_yield=nsig,
        hist_name="h_sig"
    )
    sig_mc = ws_sig.obj("ggH_Tag0")

    # ------------------------------------------------------
    # Sanity check: plot signal RooAbsPdf and histogram
    # that we got out of it
    # ------------------------------------------------------
    scale = Run2Lumi
    sig_mc4plotting = scale_roodataset_weights(sig_mc, scale)
    plot_pdf_vs_hist(
        pdf=sig_pdf,
        obs=obs,
        th1=h_sig,
        xmin=xmin,
        xmax=xmax,
        data=sig_mc4plotting,
        fit_range="full", #signal MC, so fit over full
        out_png="sig_pdf_sanity.png"
    )

    # now calculate logliklihood of null Hypothesis and MLE hypothesis,
    # which we freeze to r=1, since we're re-creating post-fit expected
    # significance
    h_sigPlusBkg = h_bkg.Clone("h_sigPlusBkg")
    h_sigPlusBkg.Add(h_sig)
    h_simulation = h_sigPlusBkg # this is the simulated truth
    h_expected = h_bkg # this is the bkg only expectation
    
    ll_bkgOnly = poisson_nll(h_simulation, h_expected)
    print(f"ll_bkgOnly: {ll_bkgOnly}")

    h_simulation = h_sigPlusBkg # this is the simulated truth
    h_expected = h_sigPlusBkg # this is the bkg only expectation
    
    ll_MLE = poisson_nll(h_simulation, h_expected)
    print(f"ll_MLE: {ll_MLE}")

    test_statistic = 2*(ll_MLE- ll_bkgOnly) # Source: eqn 34 from https://inspirehep.net/literature/1723231
    print(f"test_statistic: {test_statistic}")
    significance = math.sqrt(test_statistic) # Source: eqn 36 from https://inspirehep.net/literature/1723231
    print(f"significance: {significance}")
    print("End of program: check if the significance value above matches with combine's significance")

    # sanity check on the histograms
    plot_two_hists(h_sigPlusBkg, h_bkg)
