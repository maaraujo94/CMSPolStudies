// macro to fit the data mass background
// f, mu, alpha constant in pT
// sigma_1,2,G linear in pT
// n, fG fixed from the MC results

#import "../../ptbins.C"

double gPI = TMath::Pi();

// crystal ball function
double cb_exp(double m, double N, double sig, double m0, double n, double alpha)
{
  double delta_m = (m-m0)/sig;
  double f_val = 0;
  double norm = N/(sqrt(2*gPI)*sig);

  if(delta_m > -alpha) {
    f_val = exp(-0.5*delta_m*delta_m);
  }
  else {
    double a = abs(alpha);
    f_val = pow(n/a, n) * exp(-a*a/2.) * pow(n/a-a-delta_m, -n);
  }

  return norm * f_val;
}
// gaussian function
double g_exp(double m, double N, double sig, double m0)
{
  double delta_m = (m-m0)/sig;
  double norm = N/(sqrt(2*gPI)*sig);
  double f_val = exp(-0.5*delta_m*delta_m);

  return norm * f_val;
}
// bkg expression - negative exponential
double bkg_exp(double m, double p1, double p2)
{
  return p1 * exp( - m / p2 );
}

// crystal ball function parser - called by TF2
// parameters: NS, f, mu, sig1, sig2, n, alpha, NB, lambda, fG, sigG
double mmod_func(double *x, double *par)
{
  // get m, pt and corresp pt bin
  double m = x[0], pt = x[1];
  int pt_bin;
  for(int i = 0; i < nPtBins; i++)
    if(ptBins[i] < pt && ptBins[i+1] > pt)
      pt_bin = i;

  double NS = par[pt_bin];
  double f = par[nPtBins]; // f is constant in pt, only take the first value and use in all cases

  double mu = par[2*nPtBins]; // mu is constant in pt
  double sig1 = par[3*nPtBins] * pt + par[3*nPtBins+1]; 
  double sig2 = par[3*nPtBins] * pt + par[4*nPtBins+1]; // sigmas linear in pt
  
  double n = par[5*nPtBins]; // n is constant in pt
  double alpha = par[6*nPtBins]; // alpha is constant in pt

  double NB = par[7*nPtBins+pt_bin];
  double ld = par[8*nPtBins+pt_bin];

  double fG = par[9*nPtBins]; // fG constant in pT
  double sigG = par[3*nPtBins] * pt + par[10*nPtBins+1]; // sigma_G linear in pt - always sigma1 slope

  double func = f * cb_exp(m, NS, sig1, mu, n, alpha) + (1.-f-fG) * cb_exp(m, NS, sig2, mu, n, alpha) + fG * g_exp(m, NS, sigG, mu) + bkg_exp(m, NB, ld);

  return func;
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


// MAIN
void newDatamass_2()
{
  // PART 1 : FILLING THE MASS HISTO
  // prepare binning and histograms for plots
  TH2D *h_d2d = new TH2D();
  TFile *fin = new TFile("../../bkgFits/files/mStore.root");
  fin->GetObject("mH", h_d2d);
  h_d2d->SetDirectory(0);
  fin->Close();

  int mbins = h_d2d->GetNbinsX();
  double lowm = h_d2d->GetXaxis()->GetBinLowEdge(1);
  double him = h_d2d->GetXaxis()->GetBinUpEdge(mbins);
  
  // Make 1d histos
  TH1D **h_d1d = new TH1D*[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    h_d1d[i] = h_d2d->ProjectionX(Form("mH%.0f", ptBins[i]), i+1, i+1);
    h_d1d[i]->SetTitle(Form("Run 2 data M(#mu#mu) (%.1f < p_{T} < %.1f GeV)", ptBins[i], ptBins[i+1]));
  }

  // define aux vals for plotting
  double m_min[] = {2.94, 3.0, 3.21};
  double m_max[] = {2.95, 3.2, 3.26};

  // get alpha + unc from the prev fit
  TFile *finF = new TFile("../../bkgFits/files/mfit_2.root");
  double alpha_v = ((TGraphErrors*)finF->Get("fit_alpha"))->GetY()[0];
  double alpha_e = ((TGraphErrors*)finF->Get("fit_alpha"))->GetEY()[0];
  finF->Close();
  alpha_v+=alpha_e;

  // fix n_v to a given value, give initial alpha
  double n_v = 2.5, fG_v = 0.035;

  // define 2d function for fitting
  TF2 *f_cb = new TF2("f_cb", mmod_func, m_min[0], m_max[2], ptBins[0], ptBins[nPtBins], 11*nPtBins, 2);
  string par_n[] =  {"NS", "f",  "mu",  "sig1", "sig2", "n", "alpha", "NB", "lambda", "fG", "sigG"};
  double par_v[] =  {1.,   0.55, 3.095, 1e-4,   1e-4,   n_v, alpha_v, 1.,   0.5,      fG_v, 1.};
  double par2_v[] = {1.,   1.,   1.,    2e-2,   3e-2,   1.,  1.,      1.,   1.,       1.,   6.5e-2};
  
  // define parameters
  for(int i = 0; i < nPtBins; i++) {
    // normalizations
    f_cb->SetParName(i, Form("NS_%d", i));
    f_cb->SetParameter(i, h_d1d[i]->Integral()/100.);
    f_cb->SetParName(7*nPtBins+i, Form("NB_%d", i));
    f_cb->SetParameter(7*nPtBins+i, h_d1d[i]->Integral()/(1.5*i+1));
    
    for(int j = 1; j < 11; j++) { // between NS, NB
      if(j != 7) { // removing NB
	f_cb->SetParName(j*nPtBins+i, Form("%s_%d", par_n[j].c_str(), i));
	f_cb->SetParameter(j*nPtBins+i, par_v[j]);
	// setting the constant parameters f, mu, alpha
	if((j < 3 || j == 6)  && i > 0) f_cb->FixParameter(j*nPtBins+i, par_v[j]);
	// fixing the linear parameters sigma_1,2 - shared slope, different intercept
	else if(j == 3 && i > 1) f_cb->FixParameter(j*nPtBins+i, par_v[j]);
	else if(j == 3 && i == 1) f_cb->SetParameter(j*nPtBins+i, par2_v[j]);
	else if((j == 4 || j == 10) && i != 1) f_cb->FixParameter(j*nPtBins+i, par_v[j]);
	else if((j == 4 || j == 10) && i == 1) f_cb->SetParameter(j*nPtBins+i, par2_v[j]);
	// fixing n, fG, alpha
	else if(j == 5 || j == 6 || j == 9) f_cb->FixParameter(j*nPtBins+i, par_v[j]);
	// lambda left free
	else if(j == 8) f_cb->SetParameter(j*nPtBins+i, 0.012*(ptBins[i+1]+ptBins[i])/2.+0.18);
      }
    }
  }
  
  // fit the 2d function to the mass:pT map
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.03);
  
  f_cb->SetNpx(1000);
  h_d2d->Fit("f_cb", "RS");
  TFitResultPtr fitres = h_d2d->Fit("f_cb", "RS");

  // tf1 for plotting in the 1D bins
  TF1 *f_1d = new TF1("f_1d", "[1]*cb_exp(x,[0],[3],[2],[5],[6]) + (1.-[1]-[9]) * cb_exp(x,[0],[4],[2],[5],[6]) + [9]*g_exp(x, [0], [10], [2])+bkg_exp(x,[7],[8])", m_min[0], m_max[2]);
  f_1d->SetParNames("NS", "f", "mu", "sigma1", "sigma2", "n", "alpha", "p1", "lambda", "fG", "sigG");

  // separate parts of the fit function
  TF1 *fp1 = new TF1("fp1", "[1]*cb_exp(x,[0],[3],[2],[4],[5])", m_min[0], m_max[2]);
  fp1->SetParNames("NS", "f", "mu", "sigma1", "n", "alpha");
  TF1 *fp2 = new TF1("fp2", "(1.-[1]-[6]) * cb_exp(x,[0],[3],[2],[4],[5])", m_min[0], m_max[2]);
  fp2->SetParNames("NS", "f", "mu", "sigma2", "n", "alpha", "fG");
  TF1 *fp3 = new TF1("fp3", "bkg_exp(x,[0],[1])", m_min[0], m_max[2]);
  fp3->SetParNames("NB", "lambda");
  TF1 *fp4 = new TF1("fp4", "[1]*g_exp(x,[0],[3],[2])", m_min[0], m_max[2]);
  fp4->SetParNames("NS", "fG", "mu", "sigmaG");

  double pt_val[nPtBins], pt_err[nPtBins];
  double pars[11][nPtBins], epars[11][nPtBins];
  double fBkg[nPtBins], efz[nPtBins];

  // cycle over all pT bins
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    pt_val[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
    pt_err[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);

     // storing parameters
    for(int j = 0; j < 11; j++) {
      if(j == 0 || j == 7 || j == 8) { // free parameters NS, NB, lambda
	pars[j][i_pt] = f_cb->GetParameter(j*nPtBins+i_pt);
	epars[j][i_pt] = f_cb->GetParError(j*nPtBins+i_pt);
      }
      else if ( j == 1 || j == 2 || j == 5 || j == 6 || j == 9) { // constant parameters mu, f, n, alpha, fG
	pars[j][i_pt] = f_cb->GetParameter(j*nPtBins);
	epars[j][i_pt] = f_cb->GetParError(j*nPtBins);
      }
      else if ( j == 3 || j == 4 || j == 10) { // linear parameters sig1, sig2, sigG
	pars[j][i_pt] = f_cb->GetParameter(3*nPtBins) * pt_val[i_pt] + f_cb->GetParameter(j*nPtBins+1);
	epars[j][i_pt] = sqrt(pow(f_cb->GetParError(3*nPtBins) * pt_val[i_pt], 2) + pow(f_cb->GetParError(j*nPtBins+1), 2));
      }
    }

    c->SetTopMargin(0.04);
    c->SetLogy(0); 	
    
    // initializing f_1d and plotting
    f_1d->SetParameters(pars[0][i_pt],
			pars[1][i_pt],
			pars[2][i_pt],
			pars[3][i_pt],
			pars[4][i_pt],
			pars[5][i_pt],
			pars[6][i_pt],
			pars[7][i_pt],
			pars[8][i_pt],
			pars[9][i_pt],
			pars[10][i_pt]);
    fp1->SetParameters(pars[0][i_pt], pars[1][i_pt], pars[2][i_pt], pars[3][i_pt], pars[5][i_pt], pars[6][i_pt]);
    fp2->SetParameters(pars[0][i_pt], pars[1][i_pt], pars[2][i_pt], pars[4][i_pt], pars[5][i_pt], pars[6][i_pt], pars[9][i_pt]);
    fp3->SetParameters(pars[7][i_pt], pars[8][i_pt]);
    fp4->SetParameters(pars[0][i_pt], pars[9][i_pt], pars[2][i_pt], pars[10][i_pt]);

    h_d1d[i_pt]->SetMaximum(h_d1d[i_pt]->GetMaximum()*1.1);
    h_d1d[i_pt]->SetMinimum(0);
    h_d1d[i_pt]->SetStats(0);
    h_d1d[i_pt]->SetTitle("");
    h_d1d[i_pt]->GetYaxis()->SetTitle(Form("Events per %.0f MeV", (him-lowm)/mbins*1000));
    h_d1d[i_pt]->GetYaxis()->SetTitleOffset(1.7);
    h_d1d[i_pt]->GetXaxis()->SetTitle(Form("M(#mu#mu) (GeV)"));
    h_d1d[i_pt]->GetXaxis()->SetRangeUser(m_min[0], m_max[2]);
    h_d1d[i_pt]->SetMarkerStyle(20);
    h_d1d[i_pt]->SetMarkerSize(0.75);
    h_d1d[i_pt]->SetMarkerColor(kBlack);
    h_d1d[i_pt]->Draw("error");
    f_1d->SetLineColor(kBlue);
    f_1d->Draw("lsame");

    // tf1 for plotting in the 1D bins
    fp1->SetLineColor(kRed);
    fp1->SetLineStyle(kDashed);
    fp1->Draw("lsame");
    fp2->SetLineColor(kGreen);
    fp2->SetLineStyle(kDashed);
    fp2->Draw("lsame");
    fp3->SetLineColor(kOrange+4);
    fp3->SetLineStyle(kDashDotted);
    fp3->Draw("lsame");
    fp4->SetLineColor(kViolet);
    fp4->SetLineStyle(kDashed);
    fp4->Draw("lsame");

    TLatex lc;
    lc.SetTextSize(0.03);

    // draw CMS text    
    double xp = getPos(m_min[0], m_max[2], 0.05, 0);
    double yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.95, 0);
    lc.DrawLatex(xp, yp, "CMS");
    // draw L
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.9, 0);
    lc.DrawLatex(xp, yp, "#bf{L = 103.3 fb^{-1}}");
    // draw sqrt(s)
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.85, 0);
    lc.DrawLatex(xp, yp, "#bf{#sqrt{s} = 13 TeV}");
    // draw pT
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.75, 0);
    lc.DrawLatex(xp, yp, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV}", ptBins[i_pt], ptBins[i_pt+1]));
    // draw y
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.7, 0);
    lc.DrawLatex(xp, yp, "#bf{|#it{y}| < 1.2}");
    // draw chi2
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.6, 0);
    lc.DrawLatex(xp, yp, Form("#bf{#chi^{2}/ndf = %.0f/%d}", f_cb->GetChisquare(), f_cb->GetNDF()));
    // draw the state
    lc.SetTextSize(0.04);
    xp = getPos(m_min[0], m_max[2], 0.75, 0);
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.95, 0);
    lc.DrawLatex(xp, yp, "#bf{J/#psi}");

    TLegend *leg = new TLegend(0.7, 0.65, 1., 0.9);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(kWhite,0);
    leg->AddEntry(h_d1d[i_pt], "Data", "pl");
    leg->AddEntry(f_1d, "Total fit", "l");
    leg->AddEntry(fp1, "CB1", "l");
    leg->AddEntry(fp2, "CB2", "l");
    leg->AddEntry(fp4, "G", "l");
    leg->AddEntry(fp3, "Comb. bkg.", "l");
    leg->Draw();

    c->SaveAs(Form("plots/mass/fit_2/fit_pt%d.pdf", i_pt));
    c->Clear();

    // get the bkg fraction in the signal region (3.0 - 3.2 GeV)
    double evt_bkg = fp3->Integral(m_min[1], m_max[1]);
    
    double min_bin = h_d1d[i_pt]->GetXaxis()->FindBin(m_min[1]+1e-6);
    double max_bin = h_d1d[i_pt]->GetXaxis()->FindBin(m_max[1]-1e-6);
    double evt_all = h_d1d[i_pt]->Integral(min_bin, max_bin, "width");
    fBkg[i_pt] = evt_bkg / evt_all;
    efz[i_pt] = 0;

    // calculating pulls
    double mv[mbins], pv[mbins], dv[mbins];
    for(int i_m = 0 ; i_m < mbins; i_m++) {
      mv[i_m] = h_d1d[i_pt]->GetBinCenter(i_m+1);
      double fitv = f_1d->Eval(mv[i_m]);
      double datav = h_d1d[i_pt]->GetBinContent(i_m+1);
      double datau = h_d1d[i_pt]->GetBinError(i_m+1);
      if(datau > 0 && mv[i_m] > m_min[0] && mv[i_m] < m_max[2]) pv[i_m] = (datav-fitv)/datau;
      else pv[i_m] = 0;
      if(fitv > 0 && mv[i_m] > m_min[0] && mv[i_m] < m_max[2]) dv[i_m] = (datav-fitv)/fitv * 100.;
      else dv[i_m] = 0;

      if(abs(pv[i_m]) > 7) cout << i_pt << " " << pt_val[i_pt] << " " << i_m << " " << mv[i_m] << " " << pv[i_m] << endl;
    }
  
    c->SetTopMargin(0.1);
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

    
    c->SaveAs(Form("plots/mass/fit_2/pulls_pt%d.pdf", i_pt));
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

    c->SaveAs(Form("plots/mass/fit_2/devs_pt%d.pdf", i_pt));
    c->Clear();
  }
  
  // storing the free parameters
  TFile *foutF = new TFile("files/mfit_2.root", "recreate");
  string parlab[] = {"NS", "f", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda", "fG", "sigG"};

  for(int i_p = 0; i_p < 11; i_p++) {
    TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[i_p], pt_err, epars[i_p]);
    g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
  }

  TGraphErrors *g_fBG = new TGraphErrors(nPtBins, pt_val, fBkg, pt_err, efz);
  g_fBG->Write("fit_fBG");
  
  fitres->SetName("fitres");
  fitres->Write();
  
  foutF->Close();

  cout << f_cb->GetChisquare() << "/" << f_cb->GetNDF() << endl;

  c->Destructor();
}
