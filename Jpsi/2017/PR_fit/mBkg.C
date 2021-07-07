#import "plotDMPars.C"

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
// parameters: f, N(per bin), mu, sig1(linear), sig2(linear), n(per bin), alpha (per bin), p1 (per bin), p2 (per bin)
double mmod_func(double *x, double *par)
{
  // get m, pt and corresp pt bin
  double m = x[0], pt = x[1];
  int pt_bin;
  for(int i = 0; i < nPtBins; i++)
    if(ptBins[i] < pt && ptBins[i+1] > pt)
      pt_bin = i;

  double f = par[0];
  double mu = par[nPtBins+1];
  double sig1 = par[nPtBins+2] * pt + par[nPtBins+3];
  double sig2 = par[nPtBins+4] * pt + par[nPtBins+5];
  
  double N = par[1+pt_bin];
  double n = par[6+nPtBins+pt_bin];
  double alpha = par[6+2*nPtBins+pt_bin];
  double p1 = par[6+3*nPtBins+pt_bin];
  double p2 = par[6+4*nPtBins+pt_bin];
  
  double func = f * cb_exp(m, N, sig1, mu, n, alpha) + (1.-f) * cb_exp(m, N, sig2, mu, n, alpha) + bkg_exp(m, p1, p2);
  return func;
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


// MAIN
void mBkg()
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
  TFile *fin = new TFile("files/mStore.root");
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
  TH2D *h_d2d = new TH2D("h_d2d", "2017 data M(#mu#mu)", mbins, lowm, him, nPtBins, ptBins);

  // scale 1d histos and fill 2d histo
  for(int i = 0; i < nPtBins; i++) {
    for(int j = 0; j < mbins; j++) {
      h_d2d->SetBinContent(j+1, i+1, h_d1d[i]->GetBinContent(j+1));
      h_d2d->SetBinError(j+1, i+1, h_d1d[i]->GetBinError(j+1));
    }
  } 

  // define 2d function for fitting
  TF2 *f_cb = new TF2("f_cb", mmod_func, m_min[0], m_max[2], ptBins[0], ptBins[nPtBins], 5*nPtBins+6, 2);
  // define constant parameters - f, mu
  f_cb->SetParName(0, "f");
  f_cb->SetParameter(0, 0.7);
  f_cb->SetParName(nPtBins+1, "mu");
  f_cb->SetParameter(nPtBins+1, 3.1);
  // define linear parameters - sigma1, sigma2
  f_cb->SetParName(nPtBins+2, "m_sig1");
  f_cb->SetParameter(nPtBins+2, 1.e-4);
  f_cb->SetParName(nPtBins+3, "b_sig1");
  f_cb->SetParameter(nPtBins+3, 2.e-2);
  f_cb->SetParName(nPtBins+4, "m_sig2");
  f_cb->SetParameter(nPtBins+4, 1.e-4);
  f_cb->SetParName(nPtBins+5, "b_sig2");
  f_cb->SetParameter(nPtBins+5, 3.e-2);
  // define free parameters - N, alpha, n
  for(int i = 0; i < nPtBins; i++) {
    f_cb->SetParName(i+1, Form("NS_%d", i));
    f_cb->SetParameter(i+1, h_d1d[i]->Integral()/100.);  

    f_cb->SetParName(nPtBins+6+i, Form("n_%d", i));
    f_cb->FixParameter(nPtBins+6+i, 1.2);

    f_cb->SetParName(i+6+2*nPtBins, Form("alpha_%d", i));
    f_cb->FixParameter(i+6+2*nPtBins, 2.15);

    f_cb->SetParName(i+6+3*nPtBins, Form("NB_%d", i));
    f_cb->SetParameter(i+6+3*nPtBins, h_d1d[i]->Integral()/2.);

    f_cb->SetParName(i+6+4*nPtBins, Form("lambda_%d", i));
    f_cb->SetParameter(i+6+4*nPtBins, 0.7);
  }
  
  // fit the 2d function to the mass:pT map
  TCanvas *c = new TCanvas("", "", 700, 700);
  h_d2d->Fit("f_cb", "R");

  // tf1 for plotting in the 1D bins
  TF1 *f_1d = new TF1("f_1d", "[0]*cb_exp(x,[1],[3],[2],[5],[6]) + (1.-[0]) * cb_exp(x,[1],[4],[2],[5],[6])+bkg_exp(x,[7],[8])", m_min[0], m_max[2]);
  f_1d->SetParNames("f", "NS", "mu", "sigma1", "sigma2", "n", "alpha", "p1", "lambda");

  // separate parts of the fit function
  TF1 *fp1 = new TF1("fp1", "[0]*cb_exp(x,[1],[3],[2],[4],[5])", m_min[0], m_max[2]);
  fp1->SetParNames("f", "NS", "mu", "sigma1", "n", "alpha");
  TF1 *fp2 = new TF1("fp2", "(1.-[0]) * cb_exp(x,[1],[3],[2],[4],[5])", m_min[0], m_max[2]);
  fp2->SetParNames("f", "NS", "mu", "sigma2", "n", "alpha");
  TF1 *fp3 = new TF1("fp3", "bkg_exp(x,[0],[1])", m_min[0], m_max[2]);
  fp3->SetParNames("NB", "lambda");

  double pt_val[nPtBins], pt_err[nPtBins];
  double pars[5][nPtBins], epars[5][nPtBins];
  double fBkg[nPtBins], efz[nPtBins];

  // cycle over all pT bins
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    pt_val[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
    pt_err[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);

    // storing free parameters
    pars[0][i_pt] = f_cb->GetParameter(i_pt+1);
    pars[1][i_pt] = f_cb->GetParameter(i_pt+6+nPtBins);
    pars[2][i_pt] = f_cb->GetParameter(i_pt+6+2*nPtBins);
    pars[3][i_pt] = f_cb->GetParameter(i_pt+6+3*nPtBins);
    pars[4][i_pt] = f_cb->GetParameter(i_pt+6+4*nPtBins);
    epars[0][i_pt] = f_cb->GetParError(i_pt+1);
    epars[1][i_pt] = f_cb->GetParError(i_pt+6+nPtBins);
    epars[2][i_pt] = f_cb->GetParError(i_pt+6+2*nPtBins);
    epars[3][i_pt] = f_cb->GetParError(i_pt+6+3*nPtBins);
    epars[4][i_pt] = f_cb->GetParError(i_pt+6+4*nPtBins);

    // initializing f_1d and plotting
    f_1d->SetParameters(f_cb->GetParameter(0),
			pars[0][i_pt],
			f_cb->GetParameter(nPtBins+1),
			f_cb->GetParameter(nPtBins+2) * pt_val[i_pt] + f_cb->GetParameter(nPtBins+3),
			f_cb->GetParameter(nPtBins+4) * pt_val[i_pt] + f_cb->GetParameter(nPtBins+5),
			pars[1][i_pt],
			pars[2][i_pt],
			pars[3][i_pt],
			pars[4][i_pt]);
  
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
    fp1->SetParameters(f_1d->GetParameter(0), f_1d->GetParameter(1), f_1d->GetParameter(2), f_1d->GetParameter(3), f_1d->GetParameter(5), f_1d->GetParameter(6));
    fp1->SetLineColor(kRed);
    fp1->SetLineStyle(kDashed);
    fp1->Draw("lsame");
    fp2->SetParameters(f_1d->GetParameter(0), f_1d->GetParameter(1), f_1d->GetParameter(2), f_1d->GetParameter(4), f_1d->GetParameter(5), f_1d->GetParameter(6));
    fp2->SetLineColor(kGreen);
    fp2->SetLineStyle(kDashed);
    fp2->Draw("lsame");
    fp3->SetParameters(pars[3][i_pt], pars[4][i_pt]);
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
    
    c->SaveAs(Form("plots/mass/fit_pt%d.pdf", i_pt));
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

    
    c->SaveAs(Form("plots/mass/pulls_pt%d.pdf", i_pt));
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
  
    c->SaveAs(Form("plots/mass/devs_pt%d.pdf", i_pt));
    c->Clear();
  }
  
  // storing the free parameters
  TFile *fout = new TFile("files/mfit.root", "recreate");
  string parlab[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda"};
  int ct_p = 0, ct_pos = 0;
  int p_pos[] = {0, nPtBins+1};  

  for(int i_p = 0; i_p < 9; i_p++) {
    // free parameters: N, n, alpha
    if(i_p == 1 || i_p >= 5) {
      TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[ct_p], pt_err, epars[ct_p]);
      g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
      ct_p++;
    }
    // linear parameters: sig1, sig2
    else if(i_p == 3 || i_p == 4) {
      double par_val[nPtBins], par_err[nPtBins];
      int n_v = i_p - (i_p%2);
      for(int i = 0; i < nPtBins; i++) {
	par_val[i] = f_cb->GetParameter(nPtBins+n_v) * pt_val[i] + f_cb->GetParameter(nPtBins+n_v+1);
	par_err[i] = sqrt(pow(f_cb->GetParError(nPtBins+n_v) * pt_val[i], 2) + pow(f_cb->GetParError(nPtBins+n_v+1), 2));
      }
      TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, par_val, pt_err, par_err);
      g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
    }
    // constant parameters: f, mu
    else {
      double par_val[nPtBins], par_err[nPtBins];
      for(int i = 0; i < nPtBins; i++) {
	par_val[i] = f_cb->GetParameter(p_pos[ct_pos]);
	par_err[i] = f_cb->GetParError(p_pos[ct_pos]);
      }
      TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, par_val, pt_err, par_err);
      g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
      ct_pos++;
    }
  }
  TGraphErrors *g_fBG = new TGraphErrors(nPtBins, pt_val, fBkg, pt_err, efz);
  g_fBG->Write("fit_fBG");
  
  TLine *l_chi = new TLine(ptBins[0], f_cb->GetChisquare()/f_cb->GetNDF(), ptBins[nPtBins], f_cb->GetChisquare()/f_cb->GetNDF());
  l_chi->Write("fit_chiN");
  fout->Close();

  ofstream ftex;
  ftex.open(Form("text_output/mfit_res.tex"));
  ftex << "\\begin{tabular}{c||c|c|c|c|c||c}\n";
  ftex << "$\\pt$ (GeV) & $N_{SR}$ & $n$ & $\\alpha$ & $N_{BG}$  & $\\lambda$ (MeV) & $f_{bkg}$ (\\%) \\\\\n";
  ftex << "\\hline\n";

  for(int i = 0; i < nPtBins; i++) {
    // pT bin
    ftex << Form("$[%.0f, %.0f]$", ptBins[i], ptBins[i+1]);
    for(int i_p = 0; i_p < 5; i_p++) {
      double mult = 1.;
      if(i_p == 0 || i_p == 3) mult = 1./(ptBins[i+1]-ptBins[i]);
      else if(i_p == 4) mult = 1e3;
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
    ftex << " & " << setprecision(2) << fixed << fBkg[i]*100.;
    ftex <<  "\\\\\n";
  }
  ftex << "\\end{tabular}\n";
  ftex.close();

  ofstream fout2;
  fout2.open(Form("text_output/mfit_resA.tex"));
  fout2 << "\\begin{tabular}{c|c|cc|cc||c}\n";
  fout2 << " \\multirow{2}{*}{$f$ $(\\%)$} & \\multirow{2}{*}{$\\mu$ $(MeV)$} & \\multicolumn{2}{|c|}{$\\sigma_1$} & \\multicolumn{2}{|c||}{$\\sigma_2$} & \\multirow{2}{*}{$\\chi^2/$ndf} \\\\\n";
  fout2 << " & & $m$ ($\\times1e5$) & $b$ (MeV) & $m$ ($\\times1e5$) & $b$ (MeV) & \\\\\n";
  fout2 << "\\hline\n";

  // f
  double val = f_cb->GetParameter(p_pos[0])*100.;
  double unc = f_cb->GetParError(p_pos[0])*100.;
  int p_norm = 1.;
  if(unc < 1) p_norm = ceil(-log10(unc))+1;	
  fout2 <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
  // mu
  val = f_cb->GetParameter(p_pos[1])*1000.;
  unc = f_cb->GetParError(p_pos[1])*1000.;
  p_norm = 1.;
  if(unc < 1) p_norm = ceil(-log10(unc))+1;	
  fout2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
  // sigma_1
  val = f_cb->GetParameter(nPtBins+2)*1e5;
  unc = f_cb->GetParError(nPtBins+2)*1e5;
  p_norm = 1.;
  if(unc < 1) p_norm = ceil(-log10(unc))+1;	
  fout2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
  val = f_cb->GetParameter(nPtBins+3)*1000.;
  unc = f_cb->GetParError(nPtBins+3)*1000.;
  p_norm = 1.;
  if(unc < 1) p_norm = ceil(-log10(unc))+1;	
  fout2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
  // sigma_2 
  val = f_cb->GetParameter(nPtBins+4)*1e5;
  unc = f_cb->GetParError(nPtBins+4)*1e5;
  p_norm = 1.;
  if(unc < 1) p_norm = ceil(-log10(unc))+1;	
  fout2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
  val = f_cb->GetParameter(nPtBins+5)*1000.;
  unc = f_cb->GetParError(nPtBins+5)*1000.;
  p_norm = 1.;
  if(unc < 1) p_norm = ceil(-log10(unc))+1;	
  fout2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
  // chi^2
  fout2 << setprecision(0) << f_cb->GetChisquare() << "/" << f_cb->GetNDF() << "\\\\\n";
  fout2 << "\\end{tabular}\n";
  fout2.close();

  cout << f_cb->GetChisquare() << "/" << f_cb->GetNDF() << endl;

  c->Destructor();

  plotDMPars();
}
