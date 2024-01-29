// macro to plot and fit the psi(2S) lifetime distribution with extra exponential

TF1 *fres;
TF1 *fNP;
TF1 *fbkL, *fbkR, *fbkD;

// define negative exponential only for positive x
double pos_exp(double x, double f2, double ld1, double ld2)
{
  if(x > 0) return f2*exp(-x/ld1) + (1.-f2)*exp(-x/ld2);
  else return 0;
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

// MAIN
void ltPerPt(string lbl, int iBin)
{
  // PART 1 : GETTING THE HISTOS
  TH1D *ltHist = new TH1D();
  
  // open and read the histo store
  TFile *fin = new TFile("files/ltStore.root");
  fin->GetObject(Form("ltH_%s", lbl.c_str()), ltHist);
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
  
  // define the fit function for the NP region
  TF1 *fitS = new TF1("fitS", "[0]*pos_exp(x,[1],[2],[3])", np_lim, hit);
  fitS->SetParNames("N_NP", "f_2", "t_NP1", "t_NP2");
  fitS->SetParameters(ltHist->GetMaximum()/3., 0.5, 0.4, 0.05);
  fitS->FixParameter(1,1);
  fitS->FixParameter(3, 1);
  fitS->SetLineColor(kBlue);
  TFitResultPtr fits = ltHist->Fit(fitS, "RS");
  
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

  ltHist->SetMarkerStyle(20);
  ltHist->SetLineColor(kBlack);
  ltHist->SetMarkerColor(kBlack);
  ltHist->SetMarkerSize(.75);
  ltHist->Draw("error same");

  // plot each exponential
  TF1 *fexp1 = new TF1("fexp1", "[0]*[1]*exp(-x/[2])", np_lim, hit);
  fexp1->SetParameters(fitS->GetParameter(0), fitS->GetParameter(1), fitS->GetParameter(2));
  fexp1->SetLineColor(kGreen);
  fexp1->SetLineStyle(kDashDotted);
  fexp1->Draw("lsame");
  TF1 *fexp2 = new TF1("fexp2", "[0]*(1.-[1])*exp(-x/[2])", np_lim, hit);
  fexp2->SetParameters(fitS->GetParameter(0), fitS->GetParameter(1), fitS->GetParameter(3));
  fexp2->SetLineColor(kRed);
  fexp2->SetLineStyle(kDashDotted);
  fexp2->Draw("lsame");

 
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

  // draw mass bin
  double xp = getPos(lowt, hit, 0.35, 0);
  double yp = getPos(ltHist->GetMinimum(), ltHist->GetMaximum(), 0.9, 1);
  lc.DrawLatex(xp, yp, Form("#bf{%s}", lbl.c_str()));
  yp = getPos(ltHist->GetMinimum(), ltHist->GetMaximum(), 0.83, 1);
  lc.DrawLatex(xp, yp, "#bf{1d fit}");
  // draw the chi^2/ndf
  yp = getPos(ltHist->GetMinimum(), ltHist->GetMaximum(), 0.9, 1);
  xp = getPos(lowt, hit, 0.7, 0);
  lc.DrawLatex(xp, yp, Form("#bf{#chi^{2}/ndf = %.0f / %d}", fitS->GetChisquare(), fitS->GetNDF()));

  c->SaveAs(Form("plots/lifetime/fit_%s.pdf", lbl.c_str()));
  c->Clear();
  
  // get the pulls and rel dev distribution
  double xv[tbins], pv[tbins], dv[tbins];
  for(int i_x = 0 ; i_x < tbins; i_x++) {
    xv[i_x] = ltHist->GetBinCenter(i_x+1);
    double fitv = fitS->Eval(xv[i_x]);
    double datav = ltHist->GetBinContent(i_x+1);
    double datau = ltHist->GetBinError(i_x+1);
    if(datau > 0 && xv[i_x] > np_lim) pv[i_x] = (datav-fitv)/datau;
    else pv[i_x] = 0;
    if(fitv > 0 && xv[i_x] > np_lim) dv[i_x] = (datav-fitv)/fitv * 100.;
    else dv[i_x] = 0;
  }

  // plotting the pulls
  c->SetLogy(0);
  c->SetTopMargin(0.1);
  
  TH1F *fp = c->DrawFrame(lowPlot, -9, hit, 9);
  fp->SetXTitle("c#tau (mm)");
  fp->SetYTitle("pulls");
  fp->GetYaxis()->SetTitleOffset(1.3);
  fp->GetYaxis()->SetLabelOffset(0.01);
  fp->SetTitle(Form("Lifetime fit pulls (%s)", lbl.c_str()));
  
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

  xp = getPos(lowt, hit, 0.4, 0);
  yp = getPos(-9, 9, 0.9, 0);
  lc.DrawLatex(xp, yp, "#bf{1d fit}");
  
  c->SaveAs(Form("plots/lifetime/pulls_%s.pdf", lbl.c_str()));
  c->Clear();

  // plotting the devs
  TH1F *fd = c->DrawFrame(lowPlot, -15, hit, 15);
  fd->SetXTitle("c#tau (mm)");
  fd->SetYTitle("relative difference (%)");
  fd->GetYaxis()->SetTitleOffset(1.3);
  fd->GetYaxis()->SetLabelOffset(0.01);
  fd->SetTitle(Form("Lifetime rel. difference (%s)", lbl.c_str()));
  
  TGraph *g_dev = new TGraph(tbins, xv, dv);
  g_dev->SetLineColor(kBlack);
  g_dev->SetMarkerColor(kBlack);
  g_dev->SetMarkerStyle(20);
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

  xp = getPos(lowt, hit, 0.4, 0);
  yp = getPos(-15, 15, 0.9, 0);
  lc.DrawLatex(xp, yp, "#bf{1d fit}");

  c->SaveAs(Form("plots/lifetime/devs_%s.pdf", lbl.c_str()));
  c->Clear();

  TFile *fout = new TFile("files/ltfitres.root", "update");
  fits->SetName(Form("fitres_%s", lbl.c_str()));
  fits->Write(0, TObject::kOverwrite);
  fout->Close();

  c->Destructor();

  
}

