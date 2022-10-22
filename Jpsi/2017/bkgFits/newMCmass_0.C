double gPI = TMath::Pi();
//pt bins defined globally for access from functions
const int nPtBins = 19;
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

// fit function parser - called by TF2
// parameters: N, f, mu, sig1, sig2, n, alpha
double cb_func(double *x, double *par)
{
  // get m, pt and corresp pt bin
  double m = x[0], pt = x[1];
  int pt_bin;
  for(int i = 0; i < nPtBins; i++)
    if(ptBins[i] < pt && ptBins[i+1] > pt)
      pt_bin = i;

  double N = par[pt_bin];
  double f = par[nPtBins+pt_bin];

  double mu = par[2*nPtBins+pt_bin];
  double sig1 = par[3*nPtBins+pt_bin];
  double sig2 = par[4*nPtBins+pt_bin];
  
  double n = par[5*nPtBins+pt_bin];
  double alpha = par[6*nPtBins+pt_bin];
  
  double func = f * cb_exp(m, N, sig1, mu, n, alpha) + (1.-f) * cb_exp(m, N, sig2, mu, n, alpha);
  return func;
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


// MAIN
void newMCmass_0()
{
  // PART 1 : FILLING THE MASS HISTO
  // prepare binning and histograms for plots
  for(int i = 0; i < 10; i++) ptBins[i] = 25 + 2.5*i;
  for(int i = 0; i < 6; i++) ptBins[i+10] = 50 + 5.*i;
  for(int i = 0; i < 2; i++) ptBins[i+16] = 80 + 10.*i;
  for(int i = 0; i < 2; i++) ptBins[i+18] = 100 + 20.*i;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // prepare mass histograms
  int mbins = 80;
  double lowm = 2.9, him = 3.3;
  TH1D **h_m1d = new TH1D*[nPtBins];
  for(int ip = 0; ip < nPtBins; ip++)
    h_m1d[ip] = new TH1D(Form("mH%.0f", ptBins[ip]), Form("Full MC M(#mu#mu) (%.1f < p_{T} < %.1f)",  ptBins[ip], ptBins[ip+1]), mbins, lowm, him);
  
  TH2D *h_m2d = new TH2D("h_m2d", "Full MC M(#mu#mu)", mbins, lowm, him, nPtBins, ptBins);
 
  cout << "all MC mass histograms initialized" << endl;

  // filling all the histos at once    
  // open and read the data tree
  TFile *fin1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCO17_cos.root");
  TTree *tree1 = (TTree*)fin1->Get("MC_cos");
  TFile *fin2 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCh17_cos.root");
  TTree *tree2 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCvh17_cos.root");
  TTree *tree3 = (TTree*)fin3->Get("MC_cos");
  
  // MC 1
  Double_t mc_pt, mc_lt, mc_m, mc_y;  
  
  tree1->SetBranchAddress("dimPt", &mc_pt);
  tree1->SetBranchAddress("Rap", &mc_y);
  tree1->SetBranchAddress("Mass", &mc_m);
  tree1->SetBranchAddress("lt", &mc_lt);
    
  // cycle over data , fill the lifetime histogram
  int mEvt = tree1->GetEntries();
  for(int i = 0; i < mEvt; i++)
    {
      tree1->GetEntry(i);
      if(mc_pt > ptBins[0] && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005) {
	for(int i_p = 0; i_p < nPtBins; i_p++)
	  if(ptBins[i_p] < 46 && mc_pt > ptBins[i_p] && mc_pt < ptBins[i_p+1])
	    h_m1d[i_p]->Fill(mc_m);
      }
    }
  fin1->Close();

  // MC 2
  tree2->SetBranchAddress("dimPt", &mc_pt);
  tree2->SetBranchAddress("Rap", &mc_y);
  tree2->SetBranchAddress("Mass", &mc_m);
  tree2->SetBranchAddress("lt", &mc_lt);
    
  // cycle over data , fill the lifetime histogram
  mEvt = tree2->GetEntries();
  for(int i = 0; i < mEvt; i++)
    {
      tree2->GetEntry(i);
      if(mc_pt > 46 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005) {
	for(int i_p = 0; i_p < nPtBins; i_p++)
	  if(ptBins[i_p] > 46 && ptBins[i_p] < 66 && mc_pt > ptBins[i_p] && mc_pt < ptBins[i_p+1])
	    h_m1d[i_p]->Fill(mc_m);
      }
    }
  fin2->Close();

  // MC 3
  tree3->SetBranchAddress("dimPt", &mc_pt);
  tree3->SetBranchAddress("Rap", &mc_y);
  tree3->SetBranchAddress("Mass", &mc_m);
  tree3->SetBranchAddress("lt", &mc_lt);
    
  // cycle over data , fill the lifetime histogram
  mEvt = tree3->GetEntries();
  for(int i = 0; i < mEvt; i++)
    {
      tree3->GetEntry(i);
      if(mc_pt > 66 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005) {
	for(int i_p = 0; i_p < nPtBins; i_p++)
	  if(ptBins[i_p] > 66 && mc_pt > ptBins[i_p] && mc_pt < ptBins[i_p+1])
	    h_m1d[i_p]->Fill(mc_m);
      }
    }
  fin3->Close();

  TFile *fout_s = new TFile("files/mStore_MC.root", "recreate");
  for(int ip = 0; ip < nPtBins; ip++) {
    h_m1d[ip]->Write();	
  }
  fout_s->Close();
  
  cout << "all MC mass histograms filled" << endl << endl;

  // fill 2d histo
  for(int i = 0; i < nPtBins; i++) {
    for(int j = 0; j < mbins; j++) {
      h_m2d->SetBinContent(j+1, i+1, h_m1d[i]->GetBinContent(j+1));
      h_m2d->SetBinError(j+1, i+1, h_m1d[i]->GetBinError(j+1));
    }
  }

  double fit_i = 2.92;
  double fit_f = 3.28;

  // define 2d function for fitting
  TF2 *f_cb = new TF2("f_cb", cb_func, fit_i, fit_f, ptBins[0], ptBins[nPtBins], 7*nPtBins, 2);

  string par_n[] = {"N", "f", "mu", "sig1", "sig2", "n", "alpha"};
  double par_v[] = {1., 0.7, 3.097, 2e-2, 4e-2, 1., 2};
  // define parameters - all free
  for(int i = 0; i < nPtBins; i++) {
    f_cb->SetParName(i, Form("N_%d", i));
    f_cb->SetParameter(i, h_m1d[i]->GetMaximum()/20.);

    for(int j = 1; j < 7; j++) {
      f_cb->SetParName(j*nPtBins+i, Form("%s_%d", par_n[j].c_str(), i));
      f_cb->SetParameter(j*nPtBins+i, par_v[j]);
    }
  }


  // tf1 for plotting in the 1D bins
  TF1 *f_1d = new TF1("f_1d", "[1]*cb_exp(x,[0],[3],[2],[5],[6]) + (1.-[1]) * cb_exp(x,[0],[4],[2], [5], [6])", fit_i, fit_f);
  f_1d->SetParNames("N", "f", "mu", "sigma1", "sigma2", "n", "alpha");
  
  // separate parts of the fit function
  TF1 *fp1 = new TF1("fp1", "[1]*cb_exp(x,[0],[3],[2],[4],[5])", fit_i, fit_f);
  fp1->SetParNames("N", "f", "mu", "sigma1", "n", "alpha");
  TF1 *fp2 = new TF1("fp2", "(1.-[1]) * cb_exp(x,[0],[3],[2],[4],[5])", fit_i, fit_f);
  fp2->SetParNames("N", "f", "mu", "sigma2", "n", "alpha");

  // fit the 2d function to the mass:pT map
  TCanvas *c = new TCanvas("", "", 700, 700);

  TFile *fout = new TFile("files/MCfit_0.root", "recreate");
  h_m2d->Fit("f_cb", "R");
  
  double pt_val[nPtBins], pt_err[nPtBins];
  double pars[7][nPtBins], epars[7][nPtBins];
      
  // cycle over all pT bins
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    pt_val[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
    pt_err[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);
    
    // storing free parameters
    for(int j = 0; j < 7; j++) {
      pars[j][i_pt] = f_cb->GetParameter(j*nPtBins+i_pt);
      epars[j][i_pt] = f_cb->GetParError(j*nPtBins+i_pt);
    }
    
    c->SetLogy();
	  
    h_m1d[i_pt]->SetMaximum(h_m1d[i_pt]->GetMaximum()*1.2);
    //h_m1d[i_pt]->SetMinimum(h_m1d[i_pt]->GetMaximum()*1e-3);
    h_m1d[i_pt]->SetMinimum(h_m1d[i_pt]->GetMaximum()*1e-5);
    h_m1d[i_pt]->GetYaxis()->SetTitle(Form("Events per %.0f MeV", (him-lowm)/mbins*1000));
    h_m1d[i_pt]->GetYaxis()->SetTitleOffset(1.4);
    h_m1d[i_pt]->GetXaxis()->SetTitle(Form("M(#mu#mu) (GeV)"));
    h_m1d[i_pt]->SetMarkerStyle(20);
    h_m1d[i_pt]->SetMarkerSize(0.75);
    h_m1d[i_pt]->SetMarkerColor(kBlack);
    h_m1d[i_pt]->Draw("error");
    f_1d->SetLineColor(kBlue);
    f_1d->Draw("lsame");

    // initializing f_1d and plotting
    f_1d->SetParameters(pars[0][i_pt],
			pars[1][i_pt],
			pars[2][i_pt],
			pars[3][i_pt],
			pars[4][i_pt],
			pars[5][i_pt],
			pars[6][i_pt]);

    fp1->SetParameters(pars[0][i_pt], pars[1][i_pt], pars[2][i_pt], pars[3][i_pt], pars[5][i_pt], pars[6][i_pt]);
    fp1->SetLineColor(kRed);
    fp1->SetLineStyle(kDashed);
    fp1->Draw("lsame");
    fp2->SetParameters(pars[0][i_pt], pars[1][i_pt], pars[2][i_pt], pars[4][i_pt], pars[5][i_pt], pars[6][i_pt]);
    fp2->SetLineColor(kGreen);
    fp2->SetLineStyle(kDashed);
    fp2->Draw("lsame");
	  
    c->SaveAs(Form("plots/MCMass/fit_0/CB_pt%d.pdf", i_pt));
    c->Clear();
	  
    // calculating pulls
    double mv[mbins], pv[mbins];
    for(int i_m = 0 ; i_m < mbins; i_m++) {
      mv[i_m] = h_m1d[i_pt]->GetBinCenter(i_m+1);
      double fitv = f_1d->Eval(mv[i_m]);
      double datav = h_m1d[i_pt]->GetBinContent(i_m+1);
      double datau = h_m1d[i_pt]->GetBinError(i_m+1);
      if(datau > 0 && mv[i_m] > fit_i && mv[i_m] < fit_f) {
	pv[i_m] = (datav-fitv)/datau;
      }
      else
	pv[i_m] = 0;
    }
	  
    c->SetLogy(0);

    double lowmp = 2.94, himp = 3.2;
	  
    // plotting the puls
    TH1F *fl = c->DrawFrame(lowmp, -6, himp, 6);
    fl->SetXTitle("M(#mu#mu) (GeV)");
    fl->SetYTitle("pulls");
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("MC J/#psi Pulls (%.1f < p_{T} < %.1f)", ptBins[i_pt], ptBins[i_pt+1]));
	  
    TGraph *g_pull = new TGraph(mbins, mv, pv);
    g_pull->SetLineColor(kBlack);
    g_pull->SetMarkerColor(kBlack);
    g_pull->SetMarkerStyle(20);
    g_pull->Draw("p");
	  
    TLine *zero = new TLine(lowmp, 0, himp, 0);
    zero->SetLineStyle(kDashed);
    zero->Draw();
	  
    c->SaveAs(Form("plots/MCMass/fit_0/pulls_pt%d_dep.pdf", i_pt));
    c->Clear();

    // clean up parameters for plotting
    pars[0][i_pt] /= (ptBins[i_pt+1]-ptBins[i_pt]);
    epars[0][i_pt] /= (ptBins[i_pt+1]-ptBins[i_pt]);

    pars[1][i_pt] *= 100.;
    epars[1][i_pt] *= 100.;
    
    pars[2][i_pt] *= 1000.;
    epars[2][i_pt] *= 1000.;
    pars[3][i_pt] *= 1000.;
    epars[3][i_pt] *= 1000.;
    pars[4][i_pt] *= 1000.;
    epars[4][i_pt] *= 1000.;
  }

  // plotting the free parameters
  double parmin[] = {4.6e-3, 0, 3.09, 0.018, 0.03,  0.6, 1.9};
  double parmax[] = {5e-3,   1, 3.10, 0.03,  0.055, 2.0, 2.4};
  string partit[] = {"N", "f", "#mu", "#sigma_{1}", "#sigma_{2}", "n", "#alpha"};
  string parax[] = {"N", "f", "#mu (GeV)", "#sigma_{1} (GeV)", "#sigma_{2} (GeV)", "n", "#alpha"};
  string parlab[] = {"N", "f", "mu", "sig1", "sig2", "n", "alpha"};

  for(int i_p = 0; i_p < 7; i_p++) {
   
    TH1F *fl = c->DrawFrame(ptBins[0]-5, parmin[i_p], ptBins[nPtBins]+5, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(partit[i_p].c_str());

    TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[i_p], pt_err, epars[i_p]);
    g_par->SetLineColor(kBlue);
    g_par->SetMarkerColor(kBlue);
    g_par->Draw("p");
    g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
    
    c->Clear();
  }

  TLine *l_chiN = new TLine(ptBins[0], f_cb->GetChisquare()/f_cb->GetNDF(), ptBins[nPtBins], f_cb->GetChisquare()/f_cb->GetNDF());
  l_chiN->Write(Form("fit_chiN"));

  TLine *l_chi = new TLine(ptBins[0], f_cb->GetChisquare(), ptBins[nPtBins], f_cb->GetChisquare());
  l_chi->Write(Form("fit_chi"));
  TLine *l_ndf = new TLine(ptBins[0], f_cb->GetNDF(), ptBins[nPtBins], f_cb->GetNDF());
  l_ndf->Write(Form("fit_ndf"));

  // output fit parameters as a table
  ofstream ftex;
  ftex.open("text_output/mfit_MC_0.tex");
  ftex << "\\begin{tabular}{c||c|c|c|c|c|c|c}\n";
  ftex << "$\\pt$ (GeV) & $N$ & $f$ (\\%) & $\\mu$ (MeV) & $\\sigma_1$ (MeV) & $\\sigma_2$ (MeV) & $n$ & $\\alpha$ \\\\\n";
  ftex << "\\hline\n";

  for(int i = 0; i < nPtBins; i++) {
    // pT bin
    ftex << Form("$[%.0f, %.0f]$", ptBins[i], ptBins[i+1]);
    for(int i_p = 0; i_p < 7; i_p++) {
      double val = pars[i_p][i], unc = epars[i_p][i];
      int p_norm = 1.; 
      if(unc < 1 ) 
	p_norm = ceil(-log10(unc))+1;	
      ftex << " & " <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc;
    }
    ftex <<  "\\\\\n";
  }
  ftex << "\\end{tabular}\n";
  ftex.close();

  fout->Close();
      
  c->Destructor(); 
}
