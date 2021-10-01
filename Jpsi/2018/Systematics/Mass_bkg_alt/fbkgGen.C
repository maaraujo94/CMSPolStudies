// macro to generate distribution plots for f_bkg, to then get the uncertainty

// background function
double bkg_exp(double m, double p1, double p2)
{
  return p2 * ( - p1 * m + 1.);
}

void fbkgGen()
{
  int ngen = 1e5; // nr of iterations for generation
  gRandom = new TRandom3(0);
  
  // define aux vals for plotting
  double m_min[] = {2.94, 3.0, 3.21};
  double m_max[] = {2.95, 3.2, 3.26};
  
  TF1 *fp3 = new TF1("fp3", "bkg_exp(x,[0],[1])", m_min[0], m_max[2]);
  fp3->SetParNames("m_bkg", "b_bkg");
  
  // get fit parameters
  TFile *inBG = new TFile("files/mfit.root");
  int n_par = 2;
  TGraphErrors** g_par = new TGraphErrors*[n_par];
  for(int i = 0; i < n_par; i++) {
    inBG->GetObject(Form("fit_%s", fp3->GetParName(i)), g_par[i]);
  }
  inBG->Close();

  // get pT binning from the fit pars
  int n_pt = g_par[0]->GetN();
  double ptBins[n_pt+1];
  for(int i = 0; i < n_pt; i++) {
    ptBins[i] = g_par[0]->GetX()[i]-g_par[0]->GetEX()[i];
  }
  ptBins[n_pt] = g_par[0]->GetX()[n_pt-1]+g_par[0]->GetEX()[n_pt-1];
  
  // prepare mass histograms
  TH1D **h_d1d = new TH1D*[n_pt];
  TFile *fin = new TFile("../../PR_fit/files/mStore.root");
  for(int ip = 0; ip < n_pt; ip++) {
    fin->GetObject(Form("mH%.0f", ptBins[ip]), h_d1d[ip]);
    h_d1d[ip]->SetDirectory(0);
  }
  fin->Close();

  // cycle over all pt bins
  TH1F **h_fbg = new TH1F*[n_pt];
  double a_par[n_par];
  TFile *fout = new TFile("files/fbgDists.root", "recreate");
  for(int i_pt = 0; i_pt < n_pt; i_pt++) {
    cout << "running pt bin " << ptBins[i_pt] << " - " << ptBins[i_pt+1] << endl;
    h_fbg[i_pt] = new TH1F(Form("h_fbg_%d", i_pt), "f_BG", 100, 0, 0.1);
    for(int i_gen = 0; i_gen < ngen; i_gen++) {
      for(int i_par = 0; i_par < n_par; i_par++) {
	// generate acc'd to gaussian
	double mean = g_par[i_par]->GetY()[i_pt];
	double sigma = g_par[i_par]->GetEY()[i_pt];
	a_par[i_par] = gRandom->Gaus(mean, sigma);
      }
      // get the bkg fraction in the signal region (3.0 - 3.2 GeV)
      fp3->SetParameters(a_par[0], a_par[1]);
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
