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
// background function
double bkg_exp(double m, double p1, double p2)
{
  return p1 * exp( - m / p2 );
}


void fSBUnc()
{
  int ngen = 1e5; // nr of iterations for generation
  gRandom = new TRandom3(0);
  
  // PART 1 : the f_bg unc

  // define aux vals for plotting
  double m_min[] = {2.94, 3.0, 3.21};
  double m_max[] = {2.95, 3.2, 3.26};
  
  // fit function for each pt bin
  TF1 *f_1d = new TF1("f_1d", "[1]*cb_exp(x,[0],[3],[2],[5],[6]) + (1.-[1]-[9]) * cb_exp(x,[0],[4],[2],[5],[6]) + [9]*g_exp(x, [0], [10], [2])+bkg_exp(x,[7],[8])", m_min[0], m_max[2]);
  f_1d->SetParNames("NS", "f", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda", "fG", "sigG");
  TF1 *fp3 = new TF1("fp3", "bkg_exp(x,[0],[1])", m_min[0], m_max[2]);
  fp3->SetParNames("NB", "lambda");
  

  // get fit parameters
  TFile *inBG = new TFile("../PR_fit/files/mfit.root");
  int n_par = f_1d->GetNpar();
  TGraphErrors** g_par = new TGraphErrors*[n_par];
  for(int i = 0; i < n_par; i++) {
    inBG->GetObject(Form("fit_%s", f_1d->GetParName(i)), g_par[i]);
  }
  inBG->Close();

  int n_pt = g_par[0]->GetN();
  double ptBins[n_pt+1];
  for(int i = 0; i < n_pt; i++) {
    ptBins[i] = g_par[0]->GetX()[i]-g_par[0]->GetEX()[i];
  }
  ptBins[n_pt] = g_par[0]->GetX()[n_pt-1]+g_par[0]->GetEX()[n_pt-1];
  
  // prepare mass histograms
  TH1D **h_d1d = new TH1D*[n_pt];
  TFile *fin = new TFile("../PR_fit/files/mStore.root");
  for(int ip = 0; ip < n_pt; ip++) {
    fin->GetObject(Form("mH%.0f", ptBins[ip]), h_d1d[ip]);
    h_d1d[ip]->SetDirectory(0);
  }
  fin->Close();

  // cycle over all pt bins
  TH1F **h_fbg = new TH1F*[n_pt];
  double a_par[n_par];
  TFile *fout = new TFile("fbg_gen.root", "recreate");
  for(int i_pt = 0; i_pt < n_pt; i_pt++) {
    cout << "running pt bin " << ptBins[i_pt] << " - " << ptBins[i_pt+1] << endl;
    h_fbg[i_pt] = new TH1F(Form("h_fbg_%d", i_pt), "f_BG", 100, 0, 0.1);
    for(int i_gen = 0; i_gen < ngen; i_gen++) {
      for(int i_par = 0; i_par < n_par; i_par++) {
	// generate acc'd to gaussian
	double mean = g_par[i_par]->GetY()[i_pt];
	double sigma = g_par[i_par]->GetEY()[i_pt];
	a_par[i_par] = gRandom->Gaus(mean, sigma);
	// fill tf1 value with generated par
	f_1d->SetParameter(i_par, a_par[i_par]);      
      }
      
      // get the bkg fraction in the signal region (3.0 - 3.2 GeV)
      fp3->SetParameters(a_par[7], a_par[8]);
      double evt_bkg = fp3->Integral(m_min[1], m_max[1]);
      
      double min_bin = h_d1d[i_pt]->GetXaxis()->FindBin(m_min[1]+1e-6);
      double max_bin = h_d1d[i_pt]->GetXaxis()->FindBin(m_max[1]-1e-6);
      double evt_all = h_d1d[i_pt]->Integral(min_bin, max_bin, "width");
      h_fbg[i_pt]->Fill(evt_bkg/evt_all);
    }
    h_fbg[i_pt]->Write();
  }
  fout->Close();
  
}
