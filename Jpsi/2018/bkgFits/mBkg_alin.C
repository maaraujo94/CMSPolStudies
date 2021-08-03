double gPI = TMath::Pi();
//pt bins defined globally for access from functions
const int nPtBins = 7;
double ptBins[nPtBins+1];

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
double bkg_exp(double m, double p1, double p2)
{
  return p1 * exp( - m / p2 );
}

// crystal ball function parser - called by TF2
// parameters: NS, f, mu, sig1, sig2, n, alpha, NB, lambda
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
  double sig2 = par[4*nPtBins] * pt + par[4*nPtBins+1]; // sigmas linear in pt
  
  double n = par[5*nPtBins]; // n is constant in pt
  double alpha = par[6*nPtBins] * pt + par[6*nPtBins+1]; // alpha is linear in pt

  double NB = par[7*nPtBins+pt_bin];
  double ld = par[8*nPtBins+pt_bin];
  
  double func = f * cb_exp(m, NS, sig1, mu, n, alpha) + (1.-f) * cb_exp(m, NS, sig2, mu, n, alpha) + bkg_exp(m, NB, ld);
  return func;
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


// MAIN
void mBkg_alin()
{
  // PART 1 : FILLING THE MASS HISTO
  // prepare binning and histograms for plots
  for(int i=0; i<3; i++) ptBins[i] = 7.*i+25.;
  for(int i=0; i<4; i++) ptBins[i+3] = 46.+10.*i;
  ptBins[7] = 120;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // prepare mass histograms
  TH1D **h_d1d = new TH1D*[nPtBins];
  TFile *fin = new TFile("../PR_fit/files/mStore.root");
  for(int ip = 0; ip < nPtBins; ip++) {
    fin->GetObject(Form("mH%.0f", ptBins[ip]), h_d1d[ip]);
    h_d1d[ip]->SetDirectory(0);
  }
  fin->Close();

  int mbins = h_d1d[0]->GetNbinsX();
  double lowm = h_d1d[0]->GetXaxis()->GetBinLowEdge(1);
  double him = h_d1d[0]->GetXaxis()->GetBinUpEdge(mbins);

  // define aux vals for plotting
  double m_min[] = {2.94, 3.0, 3.21};
  double m_max[] = {2.95, 3.2, 3.26};

  // Fill 2d histo
  TH2D *h_d2d = new TH2D("h_d2d", "2018 data M(#mu#mu)", mbins, lowm, him, nPtBins, ptBins);

  // scale 1d histos and fill 2d histo
  for(int i = 0; i < nPtBins; i++) {
    for(int j = 0; j < mbins; j++) {
      h_d2d->SetBinContent(j+1, i+1, h_d1d[i]->GetBinContent(j+1));
      h_d2d->SetBinError(j+1, i+1, h_d1d[i]->GetBinError(j+1));
    }
  } 

  // define 2d function for fitting
  TF2 *f_cb = new TF2("f_cb", mmod_func, m_min[0], m_max[2], ptBins[0], ptBins[nPtBins], 9*nPtBins, 2);
  string par_n[] = {"NS", "f", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda"};
  double par_v[] = {1., 0.7, 3.1, 1e-4, 1e-4, 1., 1e-3, 1., 1.5};
  double par2_v[] = {1., 1., 1., 2e-2, 3e-2, 1., 2.1, 1., 1};
  // define parameters
  for(int i = 0; i < nPtBins; i++) {
    // normalizations - N_SR, N_BG
    f_cb->SetParName(i, Form("NS_%d", i));
    f_cb->SetParameter(i, h_d1d[i]->Integral()/100.);
    f_cb->SetParName(7*nPtBins+i, Form("NB_%d", i));
    f_cb->SetParameter(7*nPtBins+i, h_d1d[i]->Integral()/100.);

    for(int j = 1; j < 7; j++) { // between NS, NB
      f_cb->SetParName(j*nPtBins+i, Form("%s_%d", par_n[j].c_str(), i));
      f_cb->SetParameter(j*nPtBins+i, par_v[j]);
      // fixing mu, n, f so only one value matters
      if((j < 3 || j == 5 ) && i > 0) f_cb->FixParameter(j*nPtBins+i, par_v[j]);
      // setting the linear parameters sigma, alpha
      else if((j == 3 || j == 4 || j == 6) && i > 1) f_cb->FixParameter(j*nPtBins+i, par_v[j]);
      else if((j == 3 || j == 4 || j == 6) && i == 1) f_cb->SetParameter(j*nPtBins+i, par2_v[j]);
      // fixing n to MC value -- SETUP READING MC (1.18, 2.159)
      if(j == 5 && i == 0) f_cb->FixParameter(j*nPtBins+i, 1.18);
    }

    // lambda parameter
    f_cb->SetParName(8*nPtBins+i, Form("lambda_%d", i));
    f_cb->SetParameter(8*nPtBins+i, par_v[8]);
  }
  
  // fit the 2d function to the mass:pT map
  TCanvas *c = new TCanvas("", "", 700, 700);
  h_d2d->Fit("f_cb", "R");

  // tf1 for plotting in the 1D bins
  TF1 *f_1d = new TF1("f_1d", "[1]*cb_exp(x,[0],[3],[2],[5],[6]) + (1.-[1]) * cb_exp(x,[0],[4],[2],[5],[6])+bkg_exp(x,[7],[8])", m_min[0], m_max[2]);
  f_1d->SetParNames("NS", "f", "mu", "sigma1", "sigma2", "n", "alpha", "p1", "lambda");

  // separate parts of the fit function
  TF1 *fp1 = new TF1("fp1", "[1]*cb_exp(x,[0],[3],[2],[4],[5])", m_min[0], m_max[2]);
  fp1->SetParNames("NS", "f", "mu", "sigma1", "n", "alpha");
  TF1 *fp2 = new TF1("fp2", "(1.-[1]) * cb_exp(x,[0],[3],[2],[4],[5])", m_min[0], m_max[2]);
  fp2->SetParNames("NS", "f", "mu", "sigma2", "n", "alpha");
  TF1 *fp3 = new TF1("fp3", "bkg_exp(x,[0],[1])", m_min[0], m_max[2]);
  fp3->SetParNames("NB", "lambda");

  double pt_val[nPtBins], pt_err[nPtBins];
  double pars[9][nPtBins], epars[9][nPtBins];
  double fBkg[nPtBins], efz[nPtBins];

  // cycle over all pT bins
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    pt_val[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
    pt_err[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);

    // storing parameters
    for(int j = 0; j < 9; j++) {
      if(j == 0 || j == 7 || j == 8) { // free parameters NS, NB, lambda
	pars[j][i_pt] = f_cb->GetParameter(j*nPtBins+i_pt);
	epars[j][i_pt] = f_cb->GetParError(j*nPtBins+i_pt);
      }
      else if ( j == 1 || j == 2 || j == 5) { // constant parameters mu, f, n
	pars[j][i_pt] = f_cb->GetParameter(j*nPtBins);
	epars[j][i_pt] = f_cb->GetParError(j*nPtBins);
      }
      else if ( j == 3 || j == 4 || j == 6) { // linear parameters sig1, sig2, alpha
	pars[j][i_pt] = f_cb->GetParameter(j*nPtBins) * pt_val[i_pt] + f_cb->GetParameter(j*nPtBins+1);
	epars[j][i_pt] = sqrt(pow(f_cb->GetParError(j*nPtBins) * pt_val[i_pt], 2) + pow(f_cb->GetParError(j*nPtBins+1), 2));
      }
    }

    // initializing f_1d and plotting
    f_1d->SetParameters(pars[0][i_pt],
			pars[1][i_pt],
			pars[2][i_pt],
			pars[3][i_pt],
			pars[4][i_pt],
			pars[5][i_pt],
			pars[6][i_pt],
			pars[7][i_pt],
			pars[8][i_pt]);
  
    h_d1d[i_pt]->SetMaximum(h_d1d[i_pt]->GetMaximum()*1.1);
    h_d1d[i_pt]->SetMinimum(0);
    h_d1d[i_pt]->SetStats(0);
    h_d1d[i_pt]->GetYaxis()->SetTitle(Form("Events per %.0f MeV", (him-lowm)/mbins*1000));
    h_d1d[i_pt]->GetYaxis()->SetTitleOffset(1.4);
    h_d1d[i_pt]->GetXaxis()->SetTitle(Form("M(#mu#mu) (GeV)"));
    h_d1d[i_pt]->GetXaxis()->SetRangeUser(m_min[0], m_max[2]);
    h_d1d[i_pt]->SetMarkerStyle(20);
    h_d1d[i_pt]->SetMarkerSize(0.75);
    h_d1d[i_pt]->SetMarkerColor(kBlack);
    h_d1d[i_pt]->Draw("error");
    f_1d->SetLineColor(kBlue);
    f_1d->Draw("lsame");

    // tf1 for plotting in the 1D bins
    fp1->SetParameters(pars[0][i_pt], pars[1][i_pt], pars[2][i_pt], pars[3][i_pt], pars[5][i_pt], pars[6][i_pt]);
    fp1->SetLineColor(kRed);
    fp1->SetLineStyle(kDashed);
    fp1->Draw("lsame");
    fp2->SetParameters(pars[0][i_pt], pars[1][i_pt], pars[2][i_pt], pars[4][i_pt], pars[5][i_pt], pars[6][i_pt]);
    fp2->SetLineColor(kGreen);
    fp2->SetLineStyle(kDashed);
    fp2->Draw("lsame");
    fp3->SetParameters(pars[7][i_pt], pars[8][i_pt]);
    fp3->SetLineColor(kViolet);
    fp3->SetLineStyle(kDashed);
    fp3->Draw("lsame");

    TLine **lims = new TLine*[2];
    lims[0] = new TLine(m_min[1], 0, m_min[1], getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.5, 0));
    lims[1] = new TLine(m_max[1], 0, m_max[1], getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.5, 0));
    for(int j = 0; j < 2; j++) {
      lims[j]->SetLineColor(kRed);
      lims[j]->SetLineStyle(kDashed);
      lims[j]->Draw();
    }
    
    c->SaveAs(Form("plots/mass/fitl_pt%d.pdf", i_pt));
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
    }
  
    c->SetLogy(0);
    
    // plotting the pulls
    TH1F *fl = c->DrawFrame(m_min[0], -20, m_max[2], 20);
    fl->SetXTitle("M(#mu#mu) (GeV)");
    fl->SetYTitle("pulls");
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("Data J/#psi Pulls (%.0f < p_{T} < %.0f GeV)", ptBins[i_pt], ptBins[i_pt+1]));

    TGraph *g_pull = new TGraph(mbins, mv, pv);
    g_pull->SetLineColor(kBlack);
    g_pull->SetMarkerColor(kBlack);
    g_pull->SetMarkerStyle(20);
    g_pull->Draw("p");
    
    TLine *zero = new TLine(m_min[0], 0, m_max[2], 0);
    zero->SetLineStyle(kDashed);
    zero->Draw();

    TLine **limp = new TLine*[2];
    limp[0] = new TLine(m_min[1], -20, m_min[1], 20);
    limp[1] = new TLine(m_max[1], -20, m_max[1], 20);
    for(int j = 0; j < 2; j++) {
      limp[j]->SetLineColor(kRed);
      limp[j]->SetLineStyle(kDashed);
      limp[j]->Draw();
    }

    
    c->SaveAs(Form("plots/mass/pullsl_pt%d.pdf", i_pt));
    c->Clear();

    // plotting the devs
    TH1F *fd = c->DrawFrame(m_min[0], -15, m_max[2], 15);
    fd->SetXTitle("M(#mu#mu) (GeV)");
    fd->SetYTitle("relative difference (%)");
    fd->GetYaxis()->SetTitleOffset(1.3);
    fd->GetYaxis()->SetLabelOffset(0.01);
    fd->SetTitle(Form("Data J/#psi rel. difference (%.0f < p_{T} < %.0f GeV)",  ptBins[i_pt], ptBins[i_pt+1]));
  
    TGraph *g_dev = new TGraph(mbins, mv, dv);
    g_dev->SetLineColor(kBlack);
    g_dev->SetMarkerColor(kBlack);
    g_dev->SetMarkerStyle(20);
    g_dev->Draw("psame");
    
    // aux lines - pull = 0 and sigma limits
    zero->Draw("lsame");

    TLine **limd = new TLine*[2];
    limd[0] = new TLine(m_min[1], -15, m_min[1], 15);
    limd[1] = new TLine(m_max[1], -15, m_max[1], 15);
    for(int j = 0; j < 2; j++) {
      limd[j]->SetLineColor(kRed);
      limd[j]->SetLineStyle(kDashed);
      limd[j]->Draw();
    }
  
    c->SaveAs(Form("plots/mass/devsl_pt%d.pdf", i_pt));
    c->Clear();
  }
  
  // storing the free parameters
  TFile *fout = new TFile("files/mfitl.root", "recreate");
  string parlab[] = {"NS", "f", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda"};

  for(int i_p = 0; i_p < 9; i_p++) {
    TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[i_p], pt_err, epars[i_p]);
    g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
  }

  TGraphErrors *g_fBG = new TGraphErrors(nPtBins, pt_val, fBkg, pt_err, efz);
  g_fBG->Write("fit_fBG");
  
  TLine *l_chi = new TLine(ptBins[0], f_cb->GetChisquare()/f_cb->GetNDF(), ptBins[nPtBins], f_cb->GetChisquare()/f_cb->GetNDF());
  l_chi->Write("fit_chiN");
  fout->Close();

  ofstream ftex;
  ftex.open(Form("text_output/mfitl_res.tex"));
  ftex << "\\begin{tabular}{c||c|c|c||c}\n";
  ftex << "$\\pt$ (GeV) & $N_{SR}$ & $N_{BG}$ & $\\lambda$ (GeV) & $f_{bkg}$ (\\%) \\\\\n";
  ftex << "\\hline\n";

  for(int i = 0; i < nPtBins; i++) {
    // pT bin
    ftex << Form("$[%.0f, %.0f]$", ptBins[i], ptBins[i+1]);
    for(int i_p = 0; i_p < 9; i_p++) {
      if(i_p == 0 || i_p > 6) {
	double mult = 1.;
	if(i_p == 0 || i_p == 7) mult = 1./(ptBins[i+1]-ptBins[i]);
	else if(i_p != 8)  mult = 1e3;
	
	double val = pars[i_p][i]*mult, unc = epars[i_p][i]*mult;
	if (unc > 0) {
	  int p_norm = 1.; 
	  if(unc < 1 ) 
	    p_norm = ceil(-log10(unc))+1;	
	  ftex << " & " <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc;
	}
	else {
	  int p_norm = 3.;
	  ftex << " & " <<  setprecision(p_norm) << fixed << val ;
	}
      }
    }
    ftex << " & " << setprecision(2) << fixed << fBkg[i]*100.;
    ftex <<  "\\\\\n";
  }
  ftex << "\\end{tabular}\n";
  ftex.close();

  ofstream fout2;
  fout2.open(Form("text_output/mfitl_resA.tex"));
  fout2 << "\\begin{tabular}{c|c|cc|cc|c|cc||c}\n";
  fout2 << " \\multirow{2}{*}{$f$ $(\\%)$} & \\multirow{2}{*}{$\\mu$ $(MeV)$} & \\multicolumn{2}{|c|}{$\\sigma_1$} & \\multicolumn{2}{|c|}{$\\sigma_2$}  & \\multirow{2}{*}{$n$} & \\multicolumn{2}{|c||}{$\\alpha$} & \\multirow{2}{*}{$\\chi^2/$ndf} \\\\\n";
  fout2 << " & & $m$ ($\\times1e5$) & $b$ (MeV) & $m$ ($\\times1e5$) & $b$ (MeV) & & $m$ ($\\times1e4$) & $b$ & \\\\\n";
  fout2 << "\\hline\n";

  // f, mu
  for(int i = 1; i < 3; i++) {
    double mult = 1e3;
    if(i == 1) mult = 100;
    double val = pars[i][0]*mult;
    double unc = epars[i][0]*mult;
    int p_norm = 1.;
    if(unc < 1) p_norm = ceil(-log10(unc))+1;
    fout2 <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
  }
  // sigma_1, 2
  for(int i = 3; i < 5; i++) {
    double val = f_cb->GetParameter(i*nPtBins)*1e5;
    double unc = f_cb->GetParError(i*nPtBins)*1e5;
    int p_norm = 1.;
    if(unc < 1) p_norm = ceil(-log10(unc))+1;	
    fout2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
    val = f_cb->GetParameter(i*nPtBins+1)*1e3;
    unc = f_cb->GetParError(i*nPtBins+1)*1e3;
    p_norm = 1.;
    if(unc < 1) p_norm = ceil(-log10(unc))+1;	
    fout2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
  }
  // n
  for(int i = 5; i < 6; i++) {
    double val = pars[i][0];
    double unc = epars[i][0];
    int p_norm = 1.;
    if(unc > 0 && unc < 1){
      p_norm = ceil(-log10(unc))+1;
      fout2 <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
    }
    else if(unc == 0){
      p_norm = 3;
      fout2 <<  setprecision(p_norm) << fixed << val  << " & ";
    }
  }
  // alpha
  for(int i = 6; i < 7; i++) {
    double val = f_cb->GetParameter(i*nPtBins)*1e4;
    double unc = f_cb->GetParError(i*nPtBins)*1e4;
    int p_norm = 1.;
    if(unc < 1) p_norm = ceil(-log10(unc))+1;	
    fout2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
    val = f_cb->GetParameter(i*nPtBins+1);
    unc = f_cb->GetParError(i*nPtBins+1);
    p_norm = 1.;
    if(unc < 1) p_norm = ceil(-log10(unc))+1;	
    fout2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
  }
  // chi^2
  fout2 << setprecision(0) << f_cb->GetChisquare() << "/" << f_cb->GetNDF() << "\\\\\n";
  fout2 << "\\end{tabular}\n";
  fout2.close();

  cout << f_cb->GetChisquare() << "/" << f_cb->GetNDF() << endl;

  c->Destructor();
}
