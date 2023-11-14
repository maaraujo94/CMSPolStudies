// macro to fit the data mass background
// f, mu constant in pT
// sigma_1,2 linear in pT
// n, alpha fixed from the MC results

#import "ptbins.C"

double m_min[] = {2.94, 3.0, 3.24};
double m_max[] = {2.95, 3.2, 3.26};

double gPI = TMath::Pi();

int do_round(double val)
{
  int valR = (int)val;
  if (val-valR > 0.5) return valR+1;
  else return valR;
}

double bkg_exp(double m, double p1, double p2)
{
  return p1 * exp( - m / p2 );
}
double bkg_exp_plot(double m, double p1, double p2)
{
  double func = p1 * exp( - m / p2 );

  if(m > m_min[0] && m < m_max[0])
    return func;
  else if(m > m_min[2] && m < m_max[2])
    return func;
  else {
    TF1::RejectPoint();
    return 0;
  }

}


// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


// MAIN
void massFit_1d()
{
  // PART 1 : FILLING THE MASS HISTO
  // prepare binning and histograms for plots
  TH2D *h_d2d = new TH2D();
  TFile *fin = new TFile("files/mStore.root");
  fin->GetObject("mH_NP", h_d2d);
  h_d2d->SetDirectory(0);
  fin->Close();

  int mbins = h_d2d->GetNbinsX();
  double lowm = h_d2d->GetXaxis()->GetBinLowEdge(1);
  double him = h_d2d->GetXaxis()->GetBinUpEdge(mbins);
  for(int i = 0; i <= nPtBins; i++) {
    ptBins[i] = h_d2d->GetYaxis()->GetXbins()->GetArray()[i];
    cout << ptBins[i] << ",";
  }
  cout << endl;
  
  // Make 1d histos
  TH1D **h_d1d = new TH1D*[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    h_d1d[i] = h_d2d->ProjectionX(Form("mH%.0f", ptBins[i]), i+1, i+1);
    h_d1d[i]->SetTitle(Form("Run 2 data M(#mu#mu) (%.1f < p_{T} < %.1f GeV)", ptBins[i], ptBins[i+1]));
  }
  
  // setup the plotting
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.03);
  
  // tf1 for fitting
  TF1 **fp3 = new TF1*[nPtBins];
 
  double pt_val[nPtBins], pt_err[nPtBins];
  double pars[2][nPtBins], epars[2][nPtBins];
  double fBkg[nPtBins], efz[nPtBins];
  double chis[nPtBins], ndf[nPtBins], stat_d[nPtBins];
  Int_t stat[nPtBins];

  TFile *foutF = new TFile("files/mfit_bkg.root", "recreate");

  // cycle over all pT bins
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    pt_val[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
    pt_err[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);

    // run fit
    //    fp3->SetParameter(0, h_d1d[i_pt]->Integral()/(1.5*i_pt+1));
    fp3[i_pt] = new TF1(Form("fp3_%d", i_pt), "bkg_exp_plot(x,[0],[1])", m_min[0], m_max[2]);
    fp3[i_pt]->SetParNames("NB", "lambda");
    fp3[i_pt]->SetParameter(0, h_d1d[i_pt]->Integral());
    fp3[i_pt]->SetParameter(1, 0.5);

    h_d1d[i_pt]->Fit(Form("fp3_%d", i_pt), "R0");
    TFitResultPtr fitres = h_d1d[i_pt]->Fit(Form("fp3_%d", i_pt), "RS0");
    fitres->SetName(Form("fitres_%d", i_pt));
    fitres->Write();
    stat[i_pt] = fitres;
    stat_d[i_pt] = (double)stat[i_pt];
    
    // storing parameters
    for(int j = 0; j < 2; j++) {
      pars[j][i_pt] = fp3[i_pt]->GetParameter(j);
      epars[j][i_pt] = fp3[i_pt]->GetParError(j);
    }
    cout << endl << endl << fp3[i_pt]->GetChisquare() << endl << endl;
    chis[i_pt] = fp3[i_pt]->GetChisquare();
    ndf[i_pt] = (double)fp3[i_pt]->GetNDF();
    
    h_d1d[i_pt]->SetMaximum(h_d1d[i_pt]->GetMaximum()*1.1);
    h_d1d[i_pt]->SetMinimum(0);
    h_d1d[i_pt]->SetStats(0);
    h_d1d[i_pt]->GetYaxis()->SetTitle(Form("Events per %.0f MeV", (him-lowm)/mbins*1000));
    h_d1d[i_pt]->GetYaxis()->SetTitleOffset(1.7);
    h_d1d[i_pt]->GetXaxis()->SetTitle(Form("M(#mu#mu) (GeV)"));
    h_d1d[i_pt]->GetXaxis()->SetRangeUser(m_min[0], m_max[2]);
    h_d1d[i_pt]->SetMarkerStyle(20);
    h_d1d[i_pt]->SetMarkerSize(0.75);
    h_d1d[i_pt]->SetMarkerColor(kBlack);
    h_d1d[i_pt]->Draw("error");

    // tf1 for plotting in the 1D bins
    fp3[i_pt]->SetLineColor(kBlue);
    fp3[i_pt]->Draw("lsame");

    TLine *l1 = new TLine(m_max[0], 0, m_max[0], h_d1d[i_pt]->GetMaximum());
    l1->SetLineColor(kBlack);
    l1->SetLineStyle(kDashed);
    l1->Draw();

    TLine *l2 = new TLine(m_min[2], 0, m_min[2], h_d1d[i_pt]->GetMaximum());
    l2->SetLineColor(kBlack);
    l2->SetLineStyle(kDashed);
    l2->Draw();

    TLatex lc;
    lc.SetTextSize(0.03);

    // draw chi2/ndf    
    double xp = getPos(m_min[0], m_max[2], 0.05, 0);
    double yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.95, 0);
    lc.DrawLatex(xp, yp, Form("#chi^2/ndf = %.1f/%.0f", chis[i_pt], ndf[i_pt]));

    
    c->SaveAs(Form("plots/mass/fit_pt%d.pdf", i_pt));
    c->Clear();

    // get the bkg fraction in the signal region (3.0 - 3.2 GeV)
    // tf1 for plotting in the 1D bins
    TF1 *f_fb = new TF1("fp3", "bkg_exp(x,[0],[1])", m_min[0], m_max[2]);
    f_fb->SetParNames("NB", "lambda");
    f_fb->SetParameters(pars[0][i_pt], pars[1][i_pt]);
 
    double evt_bkg = f_fb->Integral(m_min[1], m_max[1]);
    
    double min_bin = h_d1d[i_pt]->GetXaxis()->FindBin(m_min[1]+1e-6);
    double max_bin = h_d1d[i_pt]->GetXaxis()->FindBin(m_max[1]-1e-6);
    double evt_all = h_d1d[i_pt]->Integral(min_bin, max_bin, "width");
    fBkg[i_pt] = evt_bkg / evt_all;
    efz[i_pt] = 0;

    // calculating pulls
    double mv[mbins], pv[mbins], dv[mbins];
    for(int i_m = 0 ; i_m < mbins; i_m++) {
      mv[i_m] = h_d1d[i_pt]->GetBinCenter(i_m+1);
      double fitv = fp3[i_pt]->Eval(mv[i_m]);
      double datav = h_d1d[i_pt]->GetBinContent(i_m+1);
      double datau = h_d1d[i_pt]->GetBinError(i_m+1);
      if(datau > 0 && ((mv[i_m] > m_min[0] && mv[i_m] < m_max[0]) || ( mv[i_m] < m_max[2] && mv[i_m] > m_min[2]))) pv[i_m] = (datav-fitv)/datau;
      else pv[i_m] = 100;
      if(fitv > 0 && ((mv[i_m] > m_min[0] && mv[i_m] < m_max[0]) || ( mv[i_m] < m_max[2] && mv[i_m] > m_min[2]))) dv[i_m] = (datav-fitv)/fitv * 100.;
      else dv[i_m] = 1000;

      //      cout << i_pt << " " << pt_val[i_pt] << " " << i_m << " " << mv[i_m] << " " << pv[i_m] << " " << dv[i_m] << endl;
    }
  
    c->SetLogy(0);
    
    // plotting the pulls
    TH1F *fl = c->DrawFrame(m_min[0], -9, m_max[2], 9);
    fl->SetXTitle("M(#mu#mu) (GeV)");
    fl->SetYTitle("pulls");
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("Data mass fit pulls (%.1f < p_{T} < %.1f GeV)", ptBins[i_pt], ptBins[i_pt+1]));

    TGraph *g_pull = new TGraph(mbins, mv, pv);
    g_pull->SetLineColor(kBlack);
    g_pull->SetMarkerColor(kBlack);
    g_pull->SetMarkerStyle(20);
    g_pull->Draw("p");
    
    TLine *zero = new TLine(m_min[0], 0, m_max[2], 0);
    zero->SetLineStyle(kDashed);
    zero->Draw();

    TLine *plim1 = new TLine(m_min[0], -5, m_max[2], -5);
    plim1->SetLineStyle(kDotted);
    plim1->Draw("lsame");
    TLine *plim2 = new TLine(m_min[0], -3, m_max[2], -3);
    plim2->SetLineStyle(kDotted);
    plim2->Draw("lsame");
    TLine *plim3 = new TLine(m_min[0], 3, m_max[2], 3);
    plim3->SetLineStyle(kDotted);
    plim3->Draw("lsame");
    TLine *plim4 = new TLine(m_min[0], 5, m_max[2], 5);
    plim4->SetLineStyle(kDotted);
    plim4->Draw("lsame");

    TLine *lp1 = new TLine(m_max[0], -9, m_max[0], 9);
    lp1->SetLineColor(kBlack);
    lp1->SetLineStyle(kDashed);
    lp1->Draw();

    TLine *lp2 = new TLine(m_min[2], -9, m_min[2], 9);
    lp2->SetLineColor(kBlack);
    lp2->SetLineStyle(kDashed);
    lp2->Draw();

    c->SaveAs(Form("plots/mass/pulls_pt%d.pdf", i_pt));
    c->Clear();

    // plotting the devs
    TH1F *fd = c->DrawFrame(m_min[0], -15, m_max[2], 15);
    fd->SetXTitle("M(#mu#mu) (GeV)");
    fd->SetYTitle("relative difference (%)");
    fd->GetYaxis()->SetTitleOffset(1.3);
    fd->GetYaxis()->SetLabelOffset(0.01);
    fd->SetTitle(Form("Data mass rel. difference (%.1f < p_{T} < %.1f GeV)",  ptBins[i_pt], ptBins[i_pt+1]));
  
    TGraph *g_dev = new TGraph(mbins, mv, dv);
    g_dev->SetLineColor(kBlack);
    g_dev->SetMarkerColor(kBlack);
    g_dev->SetMarkerStyle(20);
    g_dev->Draw("psame");
    
    // aux lines - pull = 0 and sigma limits
    zero->Draw("lsame");

    TLine *ld1 = new TLine(m_max[0], -15, m_max[0], 15);
    ld1->SetLineColor(kBlack);
    ld1->SetLineStyle(kDashed);
    ld1->Draw();

    TLine *ld2 = new TLine(m_min[2], -15, m_min[2], 15);
    ld2->SetLineColor(kBlack);
    ld2->SetLineStyle(kDashed);
    ld2->Draw();


    c->SaveAs(Form("plots/mass/devs_pt%d.pdf", i_pt));
    c->Clear();
  }
  
  // storing the free parameters
  string parlab[] = {"NB", "lambda"};

  for(int i_p = 0; i_p < 2; i_p++) {
    TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[i_p], pt_err, epars[i_p]);
    g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
  }

  TGraphErrors *g_fBG = new TGraphErrors(nPtBins, pt_val, fBkg, pt_err, efz);
  g_fBG->Write("fit_fBG");

  TGraph *g_stat = new TGraph(nPtBins, pt_val, stat_d);
  g_stat->Write("fit_status");
  TGraph *g_chi = new TGraph(nPtBins, pt_val, chis);
  g_chi->Write("fit_chisquare");
  TGraph *g_ndf = new TGraph(nPtBins, pt_val, ndf);
  g_ndf->Write("fit_ndf");
  
  foutF->Close();

  c->Destructor();
}
