// macro to plot and fit the 2018 J/psi lifetime distribution

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

// MAIN
void ltPerPt(double binLow, double binHigh, int binN)
{
  // PART 1 : GETTING THE HISTOS
  TH2D *hist2d = new TH2D();
  TH1D *ltHist = new TH1D();
  
  // open and read the histo store
  TFile *fin = new TFile("files/ltStore.root");
  fin->GetObject("ltH", hist2d);
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
  fitS->SetParNames("N_PR", "N_NP", "f", "mu", "sigma1", "sigma2", "lambda");
  fitS->SetParameters(ltHist->GetMaximum(), ltHist->GetMaximum()*5., 0.8, 0, 1e-2, 3e-2, 0.35);
  fitS->SetLineColor(kBlue);
  int fits = ltHist->Fit(fitS, "RQ");
  
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
  ltHist->SetMinimum(ltHist->GetMaximum()*5e-4);

  TH1F *fh = c->DrawFrame(lowPlot, ltHist->GetMinimum(), hit, ltHist->GetMaximum());
  fh->SetXTitle("c#tau (mm)");
  fh->SetYTitle(Form("Events per %.0f #mum", wbin*1000.));
  fh->GetYaxis()->SetTitleOffset(1.3);
  fh->GetYaxis()->SetLabelOffset(0.01);
  fh->SetTitle(ltHist->GetTitle());

  ltHist->SetMarkerStyle(20);
  ltHist->SetLineColor(kBlack);
  ltHist->SetMarkerColor(kBlack);
  ltHist->SetMarkerSize(.75);
  ltHist->Draw("error same");

  // draw fit contributions
  fres->SetParameters(fitS->GetParameter(0), fitS->GetParameter(2), fitS->GetParameter(3), fitS->GetParameter(4), fitS->GetParameter(5));
  fres->SetLineStyle(kDashDotted);
  fres->SetLineColor(kGreen);
  fres->Draw("lsame");
  fNP->SetParameters(fitS->GetParameter(1), fitS->GetParameter(2), fitS->GetParameter(3), fitS->GetParameter(4), fitS->GetParameter(5), fitS->GetParameter(6));
  fNP->SetLineStyle(kDashDotted);
  fNP->SetLineColor(kViolet);
  fNP->Draw("lsame");

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

  c->SaveAs(Form("plots/lifetime/fit_pt%.0f.pdf", binLow));
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
    if(fitv > 0 ) dv[i_x] = (datav-fitv)/fitv;
    else dv[i_x] = 0;
  }

  // plotting the pulls
  c->SetLogy(0);
  
  TH1F *fp = c->DrawFrame(lowPlot, -15, hit, 15);
  fp->SetXTitle("c#tau (mm)");
  fp->SetYTitle("pulls");
  fp->GetYaxis()->SetTitleOffset(1.3);
  fp->GetYaxis()->SetLabelOffset(0.01);
  fp->SetTitle(Form("Lifetime fit pulls (%.0f < p_{T} < %.0f GeV)", binLow, binHigh));
  
  TGraph *g_pull = new TGraph(tbins, xv, pv);
  g_pull->SetLineColor(kBlack);
  g_pull->SetMarkerColor(kBlack);
  g_pull->SetMarkerStyle(20);
  g_pull->Draw("psame");

  // aux lines - pull=0 and sigma limits
  TF1 *fcons = new TF1("fcons", "[0]", lowPlot, hit);
  fcons->SetParameter(0, 0);
  fcons->SetLineStyle(kDashed);
  fcons->SetLineColor(kBlack);
  fcons->Draw("lsame");

  TLine *psig1 = new TLine(-pr_lim, -15, -pr_lim, 15);
  psig1->SetLineStyle(kDashed);
  psig1->Draw("lsame");
  TLine *psig2 = new TLine(pr_lim, -15, pr_lim, 15);
  psig2->SetLineStyle(kDashed);
  psig2->Draw("lsame");
  TLine *psig3 = new TLine(np_lim, -15, np_lim, 15);
  psig3->SetLineStyle(kDashed);
  psig3->Draw("lsame");
  
  c->SaveAs(Form("plots/lifetime/pulls_%.0f.pdf", binLow));
  c->Clear();

  // plotting the devs
  TH1F *fd = c->DrawFrame(lowPlot, -1, hit, 1);
  fd->SetXTitle("c#tau (mm)");
  fd->SetYTitle("deviation");
  fd->GetYaxis()->SetTitleOffset(1.3);
  fd->GetYaxis()->SetLabelOffset(0.01);
  fd->SetTitle(Form("Lifetime fit deviations (%.0f < p_{T} < %.0f GeV)", binLow, binHigh));
  
  TGraph *g_dev = new TGraph(tbins, xv, dv);
  g_dev->SetLineColor(kBlack);
  g_dev->SetMarkerColor(kBlack);
  g_dev->SetMarkerStyle(20);
  g_dev->Draw("psame");

  // aux lines - pull = 0 and sigma limits
  fcons->Draw("lsame");

  TLine *dsig1 = new TLine(-pr_lim, -1, -pr_lim, 1);
  dsig1->SetLineStyle(kDashed);
  dsig1->Draw("lsame");
  TLine *dsig2 = new TLine(pr_lim, -1, pr_lim, 1);
  dsig2->SetLineStyle(kDashed);
  dsig2->Draw("lsame");
  TLine *dsig3 = new TLine(np_lim, -1, np_lim, 1);
  dsig3->SetLineStyle(kDashed);
  dsig3->Draw("lsame");
  
  c->SaveAs(Form("plots/lifetime/devs_%.0f.pdf", binLow));
  c->Clear();
  
  ofstream ftable;
  ftable.open("text_output/lt_fit.txt", std::ios::app);
  ftable << binLow << "\t " << binHigh << "\t ";
  for(int i = 0; i < 7; i++) {
    if(i < 2)
      ftable << fitS->GetParameter(i)/(binHigh-binLow) << "\t " << fitS->GetParError(i)/(binHigh-binLow) << "\t ";
    else
      ftable << fitS->GetParameter(i) << "\t " << fitS->GetParError(i) << "\t ";
  }
  ftable << fitS->GetChisquare() << "\t " << fitS->GetNDF() << "\t " << fracNP << "\n";
  ftable.close();

  c->Destructor();

}

