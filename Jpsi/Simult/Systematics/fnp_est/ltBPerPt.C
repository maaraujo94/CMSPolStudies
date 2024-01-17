// macro to plot and fit the psi(2S) lifetime distribution with extra exponential

TF1 *fres;
TF1 *fNP;
TF1 *fbkL, *fbkR, *fbkD;

// define negative exponential only for positive x
double pos_exp(double x, double ld)
{
  if(x > 0) return exp(-x/ld);
  else return 0;
}

// define final fit function summing the PR and NP contributions
double func_sum(double *xx, double *pp)
{
  double lt = xx[0];
  double N_PR = pp[0], N_NP = pp[1], f2 = pp[2], mu = pp[3], sig1 = pp[4], sig2 = pp[5], ld = pp[6];

  fres->SetParameters(N_PR, f2, mu, sig1, sig2);
  fNP->SetParameters(N_NP, f2, mu, sig1, sig2, ld);

  return fres->Eval(lt) + fNP->Eval(lt);
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

// MAIN
double ltPerPt(double binLow, double binHigh, int binN)
{
  // PART 1 : GETTING THE HISTOS
  TH2D *hist2d = new TH2D();
  TH1D *ltHist = new TH1D();
  
  // open and read the histo store
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Simult/PR_fit/files/ltStore.root");
  fin->GetObject("ltH_SR", hist2d);
  ltHist = hist2d->ProjectionX(Form("ltH%.0f", binLow), binN+1, binN+1);
  ltHist->SetDirectory(0);
  fin->Close();

  int tbins = ltHist->GetNbinsX();
  double lowt = ltHist->GetXaxis()->GetBinLowEdge(1);
  double hit = ltHist->GetXaxis()->GetBinUpEdge(tbins);
  double wbin = (hit-lowt)/(double)tbins;

  double pr_lim = 0.05;
  double np_lim = 0.1;
  double lowPlot = -0.1;
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.13);
  c->SetLogy();
  
  // PART 2 : FITTING THE HISTOS
  
  // define the resolution (=PR) function
  fres = new TF1("fres", "[0]*([1]*TMath::Gaus(x, [2],[3]) + (1.-[1])*TMath::Gaus(x, [2], [4]))", 5*lowt, 5*hit);

  // define the NP function by convolution
  TF1 *fexp = new TF1("fexp", "pos_exp(x,[0])", 5*lowt, 5*hit);
  TF1Convolution *fcNP = new TF1Convolution(fres, fexp, 5*lowt, 5*hit);
  fcNP->SetRange(5*lowt, 5*hit);
  fcNP->SetNofPointsFFT(1000);
  fNP = new TF1("fNP", *fcNP, lowt, hit, fcNP->GetNpar());

  // define the fit function as the sum of all contributions
  TF1 *fitS = new TF1("fitS", func_sum, lowt, hit, 7);
  fitS->SetParNames("N_PR", "N_NP", "f", "mu", "sigma1", "sigma2", "t_NP");
  fitS->SetParameters(ltHist->GetMaximum(), ltHist->GetMaximum()*5., 0.75, 0, 1e-2, 2.5e-2, 0.3);
  fitS->SetLineColor(kBlue);
  TFitResultPtr fits = ltHist->Fit(fitS, "RS");
  
  // get the NP fraction in the signal region (+- 100 mum)
  fNP->SetParameters(fitS->GetParameter(1), fitS->GetParameter(2), fitS->GetParameter(3), fitS->GetParameter(4), fitS->GetParameter(5), fitS->GetParameter(6));
  double evt_NP = fNP->Integral(-pr_lim, pr_lim);
    
  double min_bin = ltHist->GetXaxis()->FindBin(-(pr_lim-1e-6));
  double max_bin = ltHist->GetXaxis()->FindBin(pr_lim-1e-6);
  double evt_all = ltHist->Integral(min_bin, max_bin, "width");
  double fracNP = evt_NP / evt_all;

  // PART 3 : PLOTTING AND STORING RESULTS
  
  // draw results
  ltHist->SetMaximum(ltHist->GetMaximum()*1.2);
  ltHist->SetMinimum(ltHist->GetMaximum()*5e-3);

  TH1F *fh = c->DrawFrame(lowPlot, ltHist->GetMinimum(), hit, ltHist->GetMaximum());
  fh->SetXTitle("c#tau (mm)");
  fh->SetYTitle(Form("Events per %.0f #mum", wbin*1000.));
  fh->GetYaxis()->SetTitleOffset(1.8);
  fh->GetYaxis()->SetLabelOffset(0.01);
  fh->GetXaxis()->SetLabelOffset(0.015);
  fh->GetXaxis()->SetTitleOffset(1.3);
  //  fh->SetTitle(ltHist->GetTitle());

  ltHist->SetMarkerStyle(20);
  ltHist->SetLineColor(kBlack);
  ltHist->SetMarkerColor(kBlack);
  ltHist->SetMarkerSize(.75);
  ltHist->Draw("error same");

  // draw fit contributions
  fres->SetParameters(fitS->GetParameter(0), fitS->GetParameter(2), fitS->GetParameter(3), fitS->GetParameter(4), fitS->GetParameter(5));
  fres->SetLineStyle(kDashDotted);
  fres->SetLineColor(kGreen+3);
  fres->Draw("lsame");
  fNP->SetLineStyle(kDashDotted);
  fNP->SetLineColor(kViolet);
  fNP->Draw("lsame");

  // individual Gaussians
  TF1 **fGs = new TF1*[2];
  fGs[0] = new TF1("fGs_0", "[0]*[1]*TMath::Gaus(x, [2],[3])", 5*lowt, 5*hit);
  fGs[0]->SetParameters(fitS->GetParameter(0), fitS->GetParameter(2), fitS->GetParameter(3), fitS->GetParameter(4));
  fGs[0]->SetLineStyle(kDashDotted);
  fGs[0]->SetLineColor(kRed+2);
  fGs[0]->Draw("lsame");

  fGs[1] = new TF1("fGs_1", "[0]*(1.-[1])*TMath::Gaus(x, [2],[3])", 5*lowt, 5*hit);
  fGs[1]->SetParameters(fitS->GetParameter(0), fitS->GetParameter(2), fitS->GetParameter(3), fitS->GetParameter(5));
  fGs[1]->SetLineStyle(kDashDotted);
  fGs[1]->SetLineColor(kOrange+2);
  fGs[1]->Draw("lsame");

  // aux lines for the 2.5 sigma and 4 sigma limits
  TLine *lsig1 = new TLine(-pr_lim, ltHist->GetMinimum(), -pr_lim, ltHist->GetMaximum());
  lsig1->SetLineStyle(kDashed);
  lsig1->Draw("lsame");
  TLine *lsig2 = new TLine(pr_lim, ltHist->GetMinimum(), pr_lim, ltHist->GetMaximum());
  lsig2->SetLineStyle(kDashed);
  lsig2->Draw("lsame");
  TLine *lsig3 = new TLine(np_lim, ltHist->GetMinimum(), np_lim, ltHist->GetMaximum());
  lsig3->SetLineStyle(kDashed);
  lsig3->Draw("lsame");

  TLatex lc;
  lc.SetTextSize(0.03);

  // draw CMS text    
  double xp = getPos(lowt, hit, 0.35, 0);
  double yp = getPos(ltHist->GetMinimum(), ltHist->GetMaximum(), 0.95, 1);
  lc.DrawLatex(xp, yp, "CMS");
  // draw L
  yp = getPos(ltHist->GetMinimum(), ltHist->GetMaximum(), 0.9, 1);
  lc.DrawLatex(xp, yp, "#bf{L = 103.3 fb^{-1}}");
  // draw sqrt(s)
  yp = getPos(ltHist->GetMinimum(), ltHist->GetMaximum(), 0.85, 1);
  lc.DrawLatex(xp, yp, "#bf{#sqrt{s} = 13 TeV}");
  // draw pT
  yp = getPos(ltHist->GetMinimum(), ltHist->GetMaximum(), 0.75, 1);
  lc.DrawLatex(xp, yp, Form("#bf{%.0f < #it{p}_{T} < %.0f GeV}", binLow, binHigh));
  // draw y
  yp = getPos(ltHist->GetMinimum(), ltHist->GetMaximum(), 0.7, 1);
  lc.DrawLatex(xp, yp, "#bf{|#it{y}| < 1.2}");
  // draw the chi^2/ndf
  xp = getPos(lowt, hit, 0.7, 0);
  yp = getPos(ltHist->GetMinimum(), ltHist->GetMaximum(), 0.9, 1);
  lc.DrawLatex(xp, yp, Form("#bf{#chi^{2}/ndf = %.0f / %d}", fitS->GetChisquare(), fitS->GetNDF()));
  // draw the state
  lc.SetTextSize(0.04);
  xp = getPos(lowt, hit, 0.6, 0);
  yp = getPos(ltHist->GetMinimum(), ltHist->GetMaximum(), 0.5, 1);
  lc.DrawLatex(xp, yp, "#bf{#psi(2S)}");

  TLegend *leg = new TLegend(0.424, 0.2, 0.724, 0.5);
  //TLegend *leg = new TLegend(0.67, 0.785, 0.974, 0.985);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  leg->AddEntry(ltHist, "Data", "pl");
  leg->AddEntry(fitS, "Total fit", "l");
  leg->AddEntry(fres, "Prompt cont.", "l");
  for(int i = 0; i < 2; i++){
    leg->AddEntry(fGs[i], Form("G_{%d}", i+1), "l");
  }
  leg->AddEntry(fNP, "Non-prompt cont.", "l");
  leg->Draw();

  c->SaveAs(Form("plots/lifetimeB/fit_pt%.0f.pdf", binLow));
  c->Clear();
  
  // get the pulls and rel dev distribution
  double xv[tbins], pv[tbins], dv[tbins];
  for(int i_x = 0 ; i_x < tbins; i_x++) {
    xv[i_x] = ltHist->GetBinCenter(i_x+1);
    double fitv = fitS->Eval(xv[i_x]);
    double datav = ltHist->GetBinContent(i_x+1);
    double datau = ltHist->GetBinError(i_x+1);
    if(datau > 0 ) pv[i_x] = (datav-fitv)/datau;
    else pv[i_x] = 0;
    if(fitv > 0 ) dv[i_x] = (datav-fitv)/fitv * 100.;
    else dv[i_x] = 0;
  }

  // plotting the pulls
  c->SetTopMargin(0.1);
  c->SetLogy(0);
  
  TH1F *fp = c->DrawFrame(lowPlot, -9, hit, 9);
  fp->SetXTitle("c#tau (mm)");
  fp->SetYTitle("pulls");
  fp->GetYaxis()->SetTitleOffset(1.3);
  fp->GetYaxis()->SetLabelOffset(0.01);
  fp->SetTitle(Form("Lifetime fit pulls (%.0f < p_{T} < %.0f GeV)", binLow, binHigh));
  
  TGraph *g_pull = new TGraph(tbins, xv, pv);
  g_pull->SetLineColor(kBlack);
  g_pull->SetMarkerColor(kBlack);
  g_pull->SetMarkerStyle(20);
  //g_pull->SetMarkerSize(.75);
  g_pull->Draw("psame");

  // aux lines - pull=0 and sigma limits
  TF1 *fcons = new TF1("fcons", "[0]", lowPlot, hit);
  fcons->SetParameter(0, 0);
  fcons->SetLineStyle(kDashed);
  fcons->SetLineColor(kBlack);
  fcons->Draw("lsame");

  TLine *psig1 = new TLine(-pr_lim, -9, -pr_lim, 9);
  psig1->SetLineStyle(kDashed);
  psig1->Draw("lsame");
  TLine *psig2 = new TLine(pr_lim, -9, pr_lim, 9);
  psig2->SetLineStyle(kDashed);
  psig2->Draw("lsame");
  TLine *psig3 = new TLine(np_lim, -9, np_lim, 9);
  psig3->SetLineStyle(kDashed);
  psig3->Draw("lsame");

  TLine *plim1 = new TLine(lowPlot, -5, hit, -5);
  plim1->SetLineStyle(kDotted);
  plim1->Draw("lsame");
  TLine *plim2 = new TLine(lowPlot, -3, hit, -3);
  plim2->SetLineStyle(kDotted);
  plim2->Draw("lsame");
  TLine *plim3 = new TLine(lowPlot, 3, hit, 3);
  plim3->SetLineStyle(kDotted);
  plim3->Draw("lsame");
  TLine *plim4 = new TLine(lowPlot, 5, hit, 5);
  plim4->SetLineStyle(kDotted);
  plim4->Draw("lsame");
  
  c->SaveAs(Form("plots/lifetimeB/pulls_%.0f.pdf", binLow));
  c->Clear();

  // plotting the devs
  TH1F *fd = c->DrawFrame(lowPlot, -15, hit, 15);
  fd->SetXTitle("c#tau (mm)");
  fd->SetYTitle("relative difference (%)");
  fd->GetYaxis()->SetTitleOffset(1.3);
  fd->GetYaxis()->SetLabelOffset(0.01);
  fd->SetTitle(Form("Lifetime rel. difference (%.0f < p_{T} < %.0f GeV)", binLow, binHigh));
  
  TGraph *g_dev = new TGraph(tbins, xv, dv);
  g_dev->SetLineColor(kBlack);
  g_dev->SetMarkerColor(kBlack);
  g_dev->SetMarkerStyle(20);
  //g_dev->SetMarkerSize(.75);
  g_dev->Draw("psame");

  // aux lines - pull = 0 and sigma limits
  fcons->Draw("lsame");

  TLine *dsig1 = new TLine(-pr_lim, -15, -pr_lim, 15);
  dsig1->SetLineStyle(kDashed);
  dsig1->Draw("lsame");
  TLine *dsig2 = new TLine(pr_lim, -15, pr_lim, 15);
  dsig2->SetLineStyle(kDashed);
  dsig2->Draw("lsame");
  TLine *dsig3 = new TLine(np_lim, -15, np_lim, 15);
  dsig3->SetLineStyle(kDashed);
  dsig3->Draw("lsame");
  
  c->SaveAs(Form("plots/lifetimeB/devs_%.0f.pdf", binLow));
  c->Clear();

  TFile *fout = new TFile("files/ltfitresB.root", "update");
  fits->SetName(Form("fitres_%.0f", binLow));
  fits->Write(0, TObject::kOverwrite);
  fout->Close();

  c->Destructor();

  return fracNP;
}

