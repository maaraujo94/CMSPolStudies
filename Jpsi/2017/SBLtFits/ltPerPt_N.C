#import "mbins.C"
#import "../ptbins.C"

// define negative exponential only for positive x
double pos_exp(double x, double ld1)
{
  if(x > 0) return exp(-x/ld1);
  else return 0;
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

// define final fit function (2d)
// parameters: N_NP (exponential), tnp2 (constant)
double func_sum(double *xx, double *pp)
{
  // get m bin
  double lt = xx[0], m = xx[1];
  int m_bin = 100;
  for(int i = 0; i < nmBins; i++)
    if(m_min[i] < m && m_max[i] > m)
      m_bin = i;
  if (m_bin > nmBins){
    TF2::RejectPoint();
    return 0;
  }
  double d_m = m_max[m_bin]-m_min[m_bin];
  
  double N_NP = pp[0]*d_m*exp(-m/pp[1]); // exp
  double tnp2 = pp[1*nmBins]; // constant

  return N_NP * pos_exp(lt, tnp2);
}

// define final fit function (1d)
double sum_1d(double *xx, double *pp)
{
  double lt = xx[0];
  double N_NP = pp[0], tnp2 = pp[1];

  return N_NP * pos_exp(lt, tnp2);
}

void ltPerPt_N(double binLo, double binHi, int nBin)
{
  // prepare binning and histograms for plots
  TH2D **h_d1d_pt = new TH2D*[nmBins];
  TH1D **h_d1d = new TH1D*[nmBins];
  TFile *fin = new TFile("files/ltStore_N.root");
  for(int i = 0; i < nmBins; i++) {
    fin->GetObject(Form("ltH_%s", lbl[i].c_str()), h_d1d_pt[i]);
    h_d1d_pt[i]->SetDirectory(0);
    h_d1d[i] = h_d1d_pt[i]->ProjectionX(Form("ltH1_%s", lbl[i].c_str()), nBin+1, nBin+1);
    h_d1d[i]->SetDirectory(0);
  }
  fin->Close();

  int tbins = h_d1d[0]->GetNbinsX();
  double lowt = h_d1d[0]->GetXaxis()->GetBinLowEdge(1);
  double hit = h_d1d[0]->GetXaxis()->GetBinUpEdge(tbins);
  double wbin = (hit-lowt)/(double)tbins;

  //need to define m bin edge array for the 2d histo
  double mBins[nmBins+2];
  mBins[0] = m_min[0];
  mBins[1] = m_max[0];
  mBins[2] = m_min[1];
  mBins[3] = m_max[1];

  // get 2d histo from the sb 1d histos
  TH2D *h_d2d = new TH2D(Form("ltH_%d", nBin), "2017 data c#tau (sidebands)", tbins, lowt, hit, nmBins+1, mBins);
  for(int i_x = 0; i_x < tbins; i_x++) {
    for(int i_y = 0; i_y < 1; i_y++) {
      h_d2d->SetBinContent(i_x+1, i_y+1, h_d1d[i_y]->GetBinContent(i_x+1));
      h_d2d->SetBinError(i_x+1, i_y+1, h_d1d[i_y]->GetBinError(i_x+1));
    }
    for(int i_y = 1; i_y < nmBins; i_y++) {
      h_d2d->SetBinContent(i_x+1, i_y+2, h_d1d[i_y]->GetBinContent(i_x+1));
      h_d2d->SetBinError(i_x+1, i_y+2, h_d1d[i_y]->GetBinError(i_x+1));
    }
  }
  
  // define aux vals for plotting
  double pr_lim = 0.05;
  double np_lim = 0.1;
  double lowPlot = -0.1;

  TFile *fin_SB = new TFile("files/ltfitres2d.root");
  TFitResult* fitres_s = (TFitResult*)fin_SB->Get("fitres");
  double tnp1 = fitres_s->Parameter(nmBins*1);
  fin_SB->Close();

  // the fit function
  TF2 *fitS = new TF2("fitS", func_sum, np_lim, hit, m_min[0], m_max[nmBins-1], 2*nmBins, 2);
  string par_n[] =  {"N_NP", "t_NP"};
  double par_v[] =  {1.e6, tnp1};
  double par2_v[] = {1.,   1.};

  // define the parameters
  for(int i = 0; i < nmBins; i++) {
    for(int j = 0; j < 2; j++) {
      // name and initialize all parameters
      fitS->SetParName(j*nmBins+i, Form("%s_%d", par_n[j].c_str(), i));
      fitS->SetParameter(j*nmBins+i, par_v[j]);

      // setting the scale parameter
      if(j==0 && i == 1) fitS->SetParameter(j*nmBins+i, par2_v[j]);
      else if(j==0 && i >1) fitS->FixParameter(j*nmBins+i, par_v[j]);
      // setting the constant parameters t_NP
      else if(j==1) fitS->FixParameter(j*nmBins+i, par_v[j]);
    }
  }
  // fit the 2d function to the lifetime:pT map
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.13);
  TFitResultPtr fitres = h_d2d->Fit("fitS", "SR");

  // tf1 for plotting in the 1D bins
  TF1 *f_1d = new TF1("f_1d", sum_1d, np_lim, hit, 2);
  f_1d->SetParNames("N_NP", "tnp2");

  double m_val[nmBins], m_err[nmBins];
  double pars[2][nmBins], epars[2][nmBins];
 
  // cycle over all m bins
  for(int i_m = 0; i_m < nmBins; i_m++) {
    m_val[i_m] = 0.5*(m_max[i_m]+m_min[i_m]);
    m_err[i_m] = 0.5*(m_max[i_m]-m_min[i_m]);

    // storing all parameters - N_NP, f1, tnp_1, tnp_2
    for(int j = 0; j < 2; j++) {
      if(j ==0) { // N is exponential
	pars[j][i_m] = 2.*m_err[i_m]*fitS->GetParameter(0)*exp(-m_val[i_m]/fitS->GetParameter(1));
	epars[j][i_m] = pars[0][i_m]*sqrt(pow(fitS->GetParError(0)/fitS->GetParameter(0), 2) + pow(fitS->GetParError(1)*m_val[i_m]/(fitS->GetParameter(1)*fitS->GetParameter(1)), 2));
      }
      else{
	pars[j][i_m] = fitS->GetParameter(j*nmBins); // constant
	epars[j][i_m] = fitS->GetParError(j*nmBins); // constant
      }
    }
    
    // initializing f_1d and plotting
    f_1d->SetParameters(pars[0][i_m],
			pars[1][i_m]);
    //f_1d->Eval(0.2);

    c->SetLogy();
    c->SetTopMargin(0.015);
     
    h_d1d[i_m]->SetMaximum(h_d1d[i_m]->GetMaximum()*1.2);
    h_d1d[i_m]->SetMinimum(h_d1d[i_m]->GetMaximum()*5e-3);

    TH1F *fh = c->DrawFrame(lowPlot, h_d1d[i_m]->GetMinimum(), hit, h_d1d[i_m]->GetMaximum());
    fh->SetXTitle("c#tau (mm)");
    fh->SetYTitle(Form("Events per %.0f #mum", wbin*1000.));
    fh->GetYaxis()->SetTitleOffset(1.8);
    fh->GetYaxis()->SetLabelOffset(0.01);
    fh->GetXaxis()->SetLabelOffset(0.015);
    fh->GetXaxis()->SetTitleOffset(1.3);

    h_d1d[i_m]->SetMarkerStyle(20);
    h_d1d[i_m]->SetMarkerColor(kBlack);
    h_d1d[i_m]->SetLineColor(kBlack);
    h_d1d[i_m]->SetMarkerSize(0.75);
    h_d1d[i_m]->Draw("error same");
    
    f_1d->SetLineColor(kBlue);
    f_1d->Draw("lsame");

    // aux lines for the 2.5 sigma and 4 sigma limits
    TLine *lsig1 = new TLine(-pr_lim, h_d1d[i_m]->GetMinimum(), -pr_lim, h_d1d[i_m]->GetMaximum());
    lsig1->SetLineStyle(kDashed);
    lsig1->Draw("lsame");
    TLine *lsig2 = new TLine(pr_lim, h_d1d[i_m]->GetMinimum(), pr_lim, h_d1d[i_m]->GetMaximum());
    lsig2->SetLineStyle(kDashed);
    lsig2->Draw("lsame");
    TLine *lsig3 = new TLine(np_lim, h_d1d[i_m]->GetMinimum(), np_lim, h_d1d[i_m]->GetMaximum());
    lsig3->SetLineStyle(kDashed);
    lsig3->Draw("lsame");

    TLatex lc;
    lc.SetTextSize(0.03);

    // draw mass bin
    double xp = getPos(lowt, hit, 0.35, 0);
    double yp = getPos(h_d1d[i_m]->GetMinimum(), h_d1d[i_m]->GetMaximum(), 0.9, 1);
    lc.DrawLatex(xp, yp, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV}", binLo, binHi));
    yp = getPos(h_d1d[i_m]->GetMinimum(), h_d1d[i_m]->GetMaximum(), 0.83, 1);
    lc.DrawLatex(xp, yp, Form("#bf{%.3f < M(#mu#mu) < %.3f GeV}", m_min[i_m], m_max[i_m]));
    // draw the chi^2/ndf
    yp = getPos(h_d1d[i_m]->GetMinimum(), h_d1d[i_m]->GetMaximum(), 0.75, 1);
    lc.DrawLatex(xp, yp, Form("#bf{global #chi^{2}/ndf = %.0f / %d}", fitS->GetChisquare(), fitS->GetNDF()));

    c->SaveAs(Form("plots/lifetime_N/fit/fit_pt%.0f_%s.pdf", binLo, lbl[i_m].c_str()));
    c->Clear();
    
    // calculating pulls
    double tv[tbins], pv[tbins], dv[tbins];
    for(int i_t = 0 ; i_t < tbins; i_t++) {
      tv[i_t] = h_d1d[i_m]->GetBinCenter(i_t+1);
      double fitv = f_1d->Eval(tv[i_t]);
      double datav = h_d1d[i_m]->GetBinContent(i_t+1);
      double datau = h_d1d[i_m]->GetBinError(i_t+1);
      if(datau > 0 && tv[i_t] > np_lim) pv[i_t] = (datav-fitv)/datau;
      else pv[i_t] = 0;
      if(fitv > 0 && tv[i_t] > np_lim ) dv[i_t] = (datav-fitv)/fitv * 100.;
      else dv[i_t] = 0;
    }
    
    // plotting the pulls
    c->SetLogy(0);
    c->SetTopMargin(0.1);

    TH1F *fp = c->DrawFrame(lowPlot, -9, hit, 9);
    fp->SetXTitle("c#tau (mm)");
    fp->SetYTitle("pulls");
    fp->GetYaxis()->SetTitleOffset(1.3);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("Lifetime fit pulls"));
  
    TGraph *g_pull = new TGraph(tbins, tv, pv);
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

    xp = getPos(lowPlot, hit, 0.4, 0);
    yp = getPos(-9, 9, 0.92, 0);
    lc.DrawLatex(xp, yp, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV}", binLo, binHi));
    yp = getPos(-9, 9, 0.85, 0);
    lc.DrawLatex(xp, yp, Form("#bf{%.3f < M(#mu#mu) < %.3f GeV}", m_min[i_m], m_max[i_m]));

    c->SaveAs(Form("plots/lifetime_N/fit/pulls_pt%.0f_%s.pdf", binLo, lbl[i_m].c_str()));
    c->Clear();

    // plotting the devs
    TH1F *fd = c->DrawFrame(lowPlot, -15, hit, 15);
    fd->SetXTitle("c#tau (mm)");
    fd->SetYTitle("relative difference (%)");
    fd->GetYaxis()->SetTitleOffset(1.3);
    fd->GetYaxis()->SetLabelOffset(0.01);
    fd->SetTitle(Form("Lifetime rel. difference (%s)", lbl[i_m].c_str()));

    TGraph *g_dev = new TGraph(tbins, tv, dv);
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

    xp = getPos(lowPlot, hit, 0.4, 0);
    yp = getPos(-15, 15, 0.92, 0);
    lc.DrawLatex(xp, yp, "#bf{free N_{NP}}");
    yp = getPos(-15, 15, 0.85, 0);
    lc.DrawLatex(xp, yp, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV}", binLo, binHi));

    c->SaveAs(Form("plots/lifetime_N/fit/devs_pt%.0f_%s.pdf", binLo, lbl[i_m].c_str()));
    c->Clear();
  }
  
  // storing the free parameters
  TFile *fout = new TFile("files/ltfitres_N.root", "update");
  fitres->SetName(Form("fitres_%.0f", binLo));
  fitres->Write(0, TObject::kOverwrite);
  fout->Close();

  cout << fitS->GetChisquare() << "/" << fitS->GetNDF() << endl;
  
  c->Destructor();

}
