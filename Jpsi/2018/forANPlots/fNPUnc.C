double gPI = TMath::Pi();

// define negative exponential only for positive x
double pos_exp(double x, double ld)
{
  if(x > 0) return exp(-x/ld);
  else return 0;
}


void fNPUnc()
{
  int ngen = 1e5; // nr of iterations for generation
  gRandom = new TRandom3(0);
  
  // PART 1 : the f_bg unc

  // define aux vals for plotting
  double lowPlot = -0.1;
  double pr_lim = 0.05;
  double lowt = -0.05;
  double hit = 0.5;

  // define the resolution (=PR) function
  TF1 *fres = new TF1("fres", "[0]*([1]*TMath::Gaus(x, [2],[3]) + (1.-[1])*TMath::Gaus(x, [2], [4]))", 5*lowt, 5*hit);

  // define the NP function by convolution
  TF1 *fexp = new TF1("fexp", "pos_exp(x,[0])", 5*lowt, 5*hit);
  TF1Convolution *fcNP = new TF1Convolution(fres, fexp, 5*lowt, 5*hit);
  fcNP->SetRange(5*lowt, 5*hit);
  fcNP->SetNofPointsFFT(1000);
  TF1 *fNP = new TF1("fNP", *fcNP, lowPlot, hit, fcNP->GetNpar());
  fNP->SetParNames("N_NP", "f", "mu", "sig1", "sig2", "t");

  // get fit parameters
  TFile *inNP = new TFile("../PR_fit/files/ltfit.root");
  int n_par = fNP->GetNpar(); 
  TGraphErrors** g_par = new TGraphErrors*[n_par];
  for(int i = 0; i < n_par; i++) {
    inNP->GetObject(Form("fit_bf_%s", fNP->GetParName(i)), g_par[i]);
  }
  inNP->Close();
  double mults[] = {1, 100., 1e3, 1e3, 1e3, 1e3};

  int n_pt = g_par[0]->GetN();
  double ptBins[n_pt+1];
  for(int i = 0; i < n_pt; i++) {
    ptBins[i] = g_par[0]->GetX()[i]-g_par[0]->GetEX()[i];
  }
  ptBins[n_pt] = g_par[0]->GetX()[n_pt-1]+g_par[0]->GetEX()[n_pt-1];
  
  // prepare mass histograms
  TH1D **h_d1d = new TH1D*[n_pt];
  TFile *fin = new TFile("../PR_fit/files/ltStore.root");
  for(int ip = 0; ip < n_pt; ip++) {
    fin->GetObject(Form("ltH%.0f", ptBins[ip]), h_d1d[ip]);
    h_d1d[ip]->SetDirectory(0);
  }
  fin->Close();

  // cycle over all pt bins
  double epsrel = 1e-6, error = 0;
  TH1F **h_fnp = new TH1F*[n_pt];
  double a_par[n_par];
  TFile *fout = new TFile("fnp_gen.root", "recreate");
  for(int i_pt = 0; i_pt < n_pt; i_pt++) {
    cout << "running pt bin " << ptBins[i_pt] << " - " << ptBins[i_pt+1] << endl;
    h_fnp[i_pt] = new TH1F(Form("h_fnp_%d", i_pt), "f_NP", 100, 0.15, 0.3);

    // get total events in PR region for that pT bin
    double min_bin = h_d1d[i_pt]->GetXaxis()->FindBin(-(pr_lim-1e-6));
    double max_bin = h_d1d[i_pt]->GetXaxis()->FindBin(pr_lim-1e-6);
    double evt_all = h_d1d[i_pt]->Integral(min_bin, max_bin, "width");

    // generation cycle
    for(int i_gen = 0; i_gen < ngen; i_gen++) {
      for(int i_par = 0; i_par < n_par; i_par++) {
	// generate acc'd to gaussian
	double mean = g_par[i_par]->GetY()[i_pt];
	double sigma = g_par[i_par]->GetEY()[i_pt];
	a_par[i_par] = gRandom->Gaus(mean, sigma);
	if(i_par == 0) a_par[i_par] *= (ptBins[i_pt+1]-ptBins[i_pt]);
	a_par[i_par] /= mults[i_par];
	// fill tf1 value with generated par
	fNP->SetParameter(i_par, a_par[i_par]);      
      }
      
      // get the bkg fraction in the pr region
      double evt_bkg = fNP->IntegralOneDim(-pr_lim, pr_lim, epsrel, 1., error);
      h_fnp[i_pt]->Fill(evt_bkg/evt_all);
    }
    h_fnp[i_pt]->Write();
  }
  fout->Close();
  
}
