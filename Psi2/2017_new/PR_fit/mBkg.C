// macro to fit the data mass background
// f, mu constant in pT
// sigma_1,2 linear in pT
// n, alpha fixed from the MC results

#import "plotDMPars.C"

int do_round(double val)
{
  int valR = (int)val;
  if (val-valR > 0.5) return valR+1;
  else return valR;
}

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
// gaussian function
double g_exp(double m, double N, double sig, double m0)
{
  double delta_m = (m-m0)/sig;
  double norm = N/(sqrt(2*gPI)*sig);
  double f_val = exp(-0.5*delta_m*delta_m);

  return norm * f_val;
}
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
  double sig2 = par[4*nPtBins] * pt + par[4*nPtBins+1]; // sigmas linear in pt
  
  double n = par[5*nPtBins]; // n is constant in pt
  double alpha = par[6*nPtBins]; // alpha is constant in pt

  double NB = par[7*nPtBins+pt_bin];
  double ld = par[8*nPtBins+pt_bin];

  double fG = par[9*nPtBins]; // fG constant in pT
  double sigG = par[10*nPtBins] * pt + par[10*nPtBins+1]; // sigG linear in pT
  
  double func = f * cb_exp(m, NS, sig1, mu, n, alpha) + (1.-f-fG) * cb_exp(m, NS, sig2, mu, n, alpha) + fG * g_exp(m, NS, sigG, mu) + bkg_exp(m, NB, ld);
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
  double m_min[] = {3.4, 3.57, 3.82};
  double m_max[] = {3.52, 3.81, 4.0};

  // Fill 2d histo
  TH2D *h_d2d = new TH2D("h_d2d", "2017 data M(#mu#mu)", mbins, lowm, him, nPtBins, ptBins);

  // scale 1d histos and fill 2d histo
  for(int i = 0; i < nPtBins; i++) {
    for(int j = 0; j < mbins; j++) {
      h_d2d->SetBinContent(j+1, i+1, h_d1d[i]->GetBinContent(j+1));
      h_d2d->SetBinError(j+1, i+1, h_d1d[i]->GetBinError(j+1));
    }
  } 

  // get the MC n and alpha values for fixing
  TFile *inMC = new TFile("../bkgFits/files/MCfit_G.root");
  double n_v = ((TGraphErrors*)inMC->Get("fit_n"))->GetY()[0];
  double alpha_v = ((TGraphErrors*)inMC->Get("fit_alpha"))->GetY()[0]; 
  double fG_v = ((TGraphErrors*)inMC->Get("fit_fG"))->GetY()[0]/100.;
  double sigG_v1 = ((TGraphErrors*)inMC->Get("sigG_lin"))->GetY()[0];
  double sigG_v2 = ((TGraphErrors*)inMC->Get("sigG_lin"))->GetY()[1];
  inMC->Close();

  // fix with 3 decimal cases
  int nm = ceil(-log10(n_v))+3;	
  n_v = do_round(n_v*pow(10, nm))/pow(10, nm);
  nm = ceil(-log10(alpha_v))+3;	
  alpha_v = do_round(alpha_v*pow(10, nm))/pow(10, nm);
  nm = ceil(-log10(fG_v))+3;	
  fG_v = do_round(fG_v*pow(10, nm))/pow(10, nm);
  nm = ceil(-log10(sigG_v1))+3;	
  sigG_v1 = do_round(sigG_v1*pow(10, nm))/pow(10, nm);
  nm = ceil(-log10(sigG_v2))+3;	
  sigG_v2 = do_round(sigG_v2*pow(10, nm))/pow(10, nm);
  
  // define 2d function for fitting
  TF2 *f_cb = new TF2("f_cb", mmod_func, m_min[0], m_max[2], ptBins[0], ptBins[nPtBins], 11*nPtBins, 2);
  string par_n[] = {"NS", "f", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda", "fG", "sigG"};
  double par_v[] =  {1., 0.5, 3.7, 1e-4, 1e-4, n_v, alpha_v, 1., 0.7, fG_v, sigG_v1};
  double par2_v[] = {1., 1.,  1.,  2e-2, 3e-2, 1.,  1.,      1., 1.,  1.,   sigG_v2};
  // define parameters
  for(int i = 0; i < nPtBins; i++) {
    // normalizations
    f_cb->SetParName(i, Form("NS_%d", i));
    f_cb->SetParameter(i, h_d1d[i]->Integral()/100.);
    f_cb->SetParName(7*nPtBins+i, Form("NB_%d", i));
    f_cb->SetParameter(7*nPtBins+i, h_d1d[i]->Integral()/2.);

    for(int j = 1; j < 11; j++) { // between NS, NB
      if(j != 7) {
	f_cb->SetParName(j*nPtBins+i, Form("%s_%d", par_n[j].c_str(), i));
	f_cb->SetParameter(j*nPtBins+i, par_v[j]);
	// fixing mu, f so only one value matters
	if(j < 3  && i > 0) f_cb->FixParameter(j*nPtBins+i, par_v[j]);
	// fixing n, alpha, fG to MC value
	else if((j == 5 || j == 6 || j == 9)) f_cb->FixParameter(j*nPtBins+i, par_v[j]);

	// setting the linear parameters sigma_1,2
	else if((j == 3 || j == 4) && i > 1) f_cb->FixParameter(j*nPtBins+i, par_v[j]);
	else if((j == 3 || j == 4) && i == 1) f_cb->SetParameter(j*nPtBins+i, par2_v[j]);
	// setting sigma_G from MC
	else if((j == 10) && i != 1) f_cb->FixParameter(j*nPtBins+i, par_v[j]);
	else if((j == 10) && i == 1) f_cb->FixParameter(j*nPtBins+i, par2_v[j]);
      }
    }
  }
  
  // fit the 2d function to the mass:pT map
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLeftMargin(0.12);
  f_cb->SetNpx(1000);
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
			pars[8][i_pt],
			pars[9][i_pt],
			pars[10][i_pt]);
  
    h_d1d[i_pt]->SetMaximum(h_d1d[i_pt]->GetMaximum()*1.1);
    h_d1d[i_pt]->SetMinimum(0);
    h_d1d[i_pt]->SetStats(0);
    h_d1d[i_pt]->GetYaxis()->SetTitle(Form("Events per %.0f MeV", (him-lowm)/mbins*1000));
    h_d1d[i_pt]->GetYaxis()->SetTitleOffset(1.8);
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
    fp2->SetParameters(pars[0][i_pt], pars[1][i_pt], pars[2][i_pt], pars[4][i_pt], pars[5][i_pt], pars[6][i_pt], pars[9][i_pt]);
    fp2->SetLineColor(kGreen);
    fp2->SetLineStyle(kDashed);
    fp2->Draw("lsame");
    fp3->SetParameters(pars[7][i_pt], pars[8][i_pt]);
    fp3->SetLineColor(kOrange+4);
    fp3->SetLineStyle(kDashDotted);
    fp3->Draw("lsame");
    fp4->SetParameters(pars[0][i_pt], pars[9][i_pt], pars[2][i_pt], pars[10][i_pt]);
    fp4->SetLineColor(kViolet);
    fp4->SetLineStyle(kDashed);
    fp4->Draw("lsame");

    TLine **lims = new TLine*[2];
    lims[0] = new TLine(m_min[1], 0, m_min[1], getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.5, 0));
    lims[1] = new TLine(m_max[1], 0, m_max[1], getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.5, 0));
    for(int j = 0; j < 2; j++) {
      lims[j]->SetLineColor(kRed);
      lims[j]->SetLineStyle(kDashed);
      //lims[j]->Draw();
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

      if(abs(pv[i_m]) > 7) cout << i_pt << " " << pt_val[i_pt] << " " << i_m << " " << mv[i_m] << " " << pv[i_m] << endl;
    }
  
    c->SetLogy(0);
    
    // plotting the pulls
    TH1F *fl = c->DrawFrame(m_min[0], -9, m_max[2], 9);
    fl->SetXTitle("M(#mu#mu) (GeV)");
    fl->SetYTitle("pulls");
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("Data mass fit pulls (%.0f < p_{T} < %.0f GeV)", ptBins[i_pt], ptBins[i_pt+1]));

    TGraph *g_pull = new TGraph(mbins, mv, pv);
    g_pull->SetLineColor(kBlack);
    g_pull->SetMarkerColor(kBlack);
    g_pull->SetMarkerStyle(20);
    g_pull->Draw("p");
    
    TLine *zero = new TLine(m_min[0], 0, m_max[2], 0);
    zero->SetLineStyle(kDashed);
    zero->Draw();

    /*    TLine **limp = new TLine*[2];
    limp[0] = new TLine(m_min[1], -20, m_min[1], 20);
    limp[1] = new TLine(m_max[1], -20, m_max[1], 20);
    for(int j = 0; j < 2; j++) {
      limp[j]->SetLineColor(kRed);
      limp[j]->SetLineStyle(kDashed);
      limp[j]->Draw();
      }*/

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
  string parlab[] = {"NS", "f", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda", "fG", "sigG"};

  for(int i_p = 0; i_p < 11; i_p++) {
    TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[i_p], pt_err, epars[i_p]);
    g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
  }

  TGraphErrors *g_fBG = new TGraphErrors(nPtBins, pt_val, fBkg, pt_err, efz);
  g_fBG->Write("fit_fBG");
  
  TLine *l_chi = new TLine(ptBins[0], f_cb->GetChisquare()/f_cb->GetNDF(), ptBins[nPtBins], f_cb->GetChisquare()/f_cb->GetNDF());
  l_chi->Write("fit_chiN");

  fitres->SetName("fitres");
  fitres->Write();
  
  fout->Close();

  double mult[] = {1., 1e2, 1e3, 1e3, 1e3, 1., 1., 1., 1., 1e2, 1e3};

  // tex table with values per pt bin
  ofstream ftex;
  ftex.open(Form("text_output/mfit_res.tex"));
  ftex << "\\begin{tabular}{c||c|c|c|c|c|c|c|c|c|c|c||c}\n";
  ftex << "$\\pt$ (GeV) & $N_{SR}$ & $f_{CB1}$ (\\%) & $\\mu$ (MeV) & $\\sigma_1$ (MeV) & $\\sigma_2$ (MeV) & $n$ & $\\alpha$ & $N_{BG}$ & $\\lambda$ (GeV) & $f_G$ (\\%) & $\\sigma_G$ (MeV) & $f_{bkg}$ (\\%) \\\\\n";
  ftex << "\\hline\n";

  for(int i = 0; i < nPtBins; i++) {
    // pT bin
    ftex << Form("$[%.0f, %.0f]$", ptBins[i], ptBins[i+1]);
    for(int i_p = 0; i_p < 11; i_p++) {
      // plot all pT values - N (0), sig1,2 (3,4), N_BG (7), lambda (8), sigG (10)
      if(i_p == 0 || i_p == 3 || i_p == 4 || i_p == 7 || i_p == 8 || i_p == 10) {
	double val = pars[i_p][i]*mult[i_p], unc = epars[i_p][i]*mult[i_p];
	if(i_p == 0 || i_p == 7) {
	  val /= (ptBins[i+1]-ptBins[i]);
	  unc /= (ptBins[i+1]-ptBins[i]);
	}
	if(unc > 0) {
	  int p_norm = 1.; 
	  if(unc < 1) p_norm = ceil(-log10(unc))+1;	
	  ftex << " & " << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc;
	}
	else {
	  int p_norm = 2.;
	  ftex << " & " <<  setprecision(p_norm) << fixed << val ;
	}
      }
      // plot single value: f (1), mu (2), n, alpha (5,6), fG (9)
      else if((i_p == 1 || i_p == 2 || i_p == 5 || i_p == 6 || i_p == 9) && i == 0) {
	double val = pars[i_p][i]*mult[i_p], unc = epars[i_p][i]*mult[i_p];
	if(unc > 0) {
	  int p_norm = 1.; 
	  if(unc < 1) p_norm = ceil(-log10(unc))+1;	
	  ftex << " & \\multirow{" << nPtBins << "}{*}{" <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << "}" ;
	}
	else {
	  int p_norm = 3.;
	  ftex << " & \\multirow{" << nPtBins << "}{*}{" <<  setprecision(p_norm) << fixed << val << "}" ;
	}
      }
      else
	ftex << " & ";
    }
    ftex << " & " << setprecision(2) << fixed << fBkg[i]*100.;
    ftex <<  "\\\\\n";
  }
  ftex << "\\end{tabular}\n";
  ftex.close();
  
  // sigma parameters
  ofstream ftex2;
  ftex2.open("text_output/mfit_resA.tex");
  ftex2 << "\\begin{tabular}{cc|cc|cc||c}\n";
  ftex2 << "\\multicolumn{2}{c|}{$\\sigma_1$} & \\multicolumn{2}{|c}{$\\sigma_2$} & \\multicolumn{2}{|c}{$\\sigma_G$}  & \\multirow{2}{*}{$\\chi^2/$ndf}\\\\\n";
  ftex2 << "$m$ ($\\times1e5$) & $b$ (MeV) & $m$ ($\\times1e5$) & $b$ (MeV) & $m$ ($\\times1e5$) & $b$ (MeV) & \\\\\n";
  ftex2 << "\\hline\n";

  for(int j = 3; j < 5; j++) {
    double val = f_cb->GetParameter(j*nPtBins)*1e5;
    double unc = f_cb->GetParError(j*nPtBins)*1e5;
    int p_norm = ceil(-log10(unc))+1;	
    ftex2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
    val = f_cb->GetParameter(j*nPtBins+1)*1e3;
    unc = f_cb->GetParError(j*nPtBins+1)*1e3;
    p_norm = ceil(-log10(unc))+1;	
    ftex2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc;
    ftex2 << " & ";
  }
  for(int j = 10; j < 11; j++) {
    double val = f_cb->GetParameter(j*nPtBins)*1e5;
    int p_norm = ceil(-log10(val))+3;	
    ftex2 << setprecision(p_norm) << fixed << val << " & ";
    val = f_cb->GetParameter(j*nPtBins+1)*1e3;
    p_norm = ceil(-log10(val))+3;	
    ftex2 << setprecision(p_norm) << fixed << val << " & ";
  }
  // chi^2
  ftex2 << setprecision(0) << f_cb->GetChisquare() << "/" << f_cb->GetNDF() << "\\\\\n";
  ftex2 << "\\end{tabular}\n";
  ftex2.close();

  cout << f_cb->GetChisquare() << "/" << f_cb->GetNDF() << endl;

  c->Destructor();

  plotDMPars();
}
