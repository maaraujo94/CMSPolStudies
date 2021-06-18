double gPI = TMath::Pi();
//pt bins defined globally for access from functions
const int nPtBins = 18;
double ptBins[nPtBins+1];
int DO_FILL = 0;

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

// crystal ball function parser - called by TF2
// parameters: f, N(per bin), mu, sig1(per bin), sig2(per bin), n(per bin), alpha (per bin)
double cb_func(double *x, double *par)
{
  // get m, pt and corresp pt bin
  double m = x[0], pt = x[1];
  int pt_bin;
  for(int i = 0; i < nPtBins; i++)
    if(ptBins[i] < pt && ptBins[i+1] > pt)
      pt_bin = i;

  double f = par[0];
  double mu = par[nPtBins+1];
  
  double N = par[1+pt_bin];
  double sig1 = par[2+nPtBins+pt_bin];
  double sig2 = par[2+2*nPtBins+pt_bin];
  double n = par[2+3*nPtBins+pt_bin];
  double alpha = par[2+4*nPtBins+pt_bin];
  
  double func = f * cb_exp(m, N, sig1, mu, n, alpha) + (1.-f) * cb_exp(m, N, sig2, mu, n, alpha);
  return func;
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


// MAIN
void MCmass_mixed()
{
  // PART 1 : FILLING THE MASS HISTO
  // prepare binning and histograms for plots
  for(int i = 0; i < 1; i++) ptBins[i] = i+25.;
  for(int i = 0; i < 7; i++) ptBins[i+1] = 26.+i*2.;
  for(int i = 0; i < 2; i++) ptBins[i+8] = 40.+i*3;
  for(int i = 0; i < 1; i++) ptBins[i+10] = 46.+i*4;
  for(int i = 0; i < 4; i++) ptBins[i+11] = 50.+i*5.;
  ptBins[14] = 66;
  for(int i = 0; i < 4; i++) ptBins[i+15] = 70.+i*10.;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // prepare mass histograms
  int mbins = 80;
  double lowm = 2.9, him = 3.3;
  TH1D **h_m1d = new TH1D*[nPtBins];
  for(int ip = 0; ip < nPtBins; ip++)
    h_m1d[ip] = new TH1D(Form("mH%.0f", ptBins[ip]), Form("2018 MC M(#mu#mu) (%.0f < p_{T} < %.0f)",  ptBins[ip], ptBins[ip+1]), mbins, lowm, him);
  
  TH2D *h_m2d = new TH2D("h_m2d", "2018 MC M(#mu#mu)", mbins, lowm, him, nPtBins, ptBins);
 
  cout << "all MC mass histograms initialized" << endl;

  if(DO_FILL == 1) {
    // filling all the histos at once    
    // open and read the data tree
    TFile *fin1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MC18_cos.root");
    TTree *tree1 = (TTree*)fin1->Get("MC_cos");
    TFile *fin2 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCh18_cos.root");
    TTree *tree2 = (TTree*)fin2->Get("MC_cos");
    TFile *fin3 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCvh18_cos.root");
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
	if(mc_pt > 25 && mc_pt < 46 && mc_lt > -0.1 && mc_lt < 0.1) {
	  for(int i_p = 0; i_p < nPtBins; i_p++)
	    if(mc_pt > ptBins[i_p] && mc_pt < ptBins[i_p+1])
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
	if(mc_pt > 46 && mc_pt < 66 && mc_lt > -0.1 && mc_lt < 0.1) {
	  for(int i_p = 0; i_p < nPtBins; i_p++)
	    if(mc_pt > ptBins[i_p] && mc_pt < ptBins[i_p+1])
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
	if(mc_pt > 66 && mc_pt < 100 && mc_lt > -0.1 && mc_lt < 0.1) {
	  for(int i_p = 0; i_p < nPtBins; i_p++)
	    if(mc_pt > ptBins[i_p] && mc_pt < ptBins[i_p+1])
	      h_m1d[i_p]->Fill(mc_m);
	}
      }
    fin3->Close();

    TFile *fout = new TFile("files/mStore_MC.root", "recreate");
    for(int ip = 0; ip < nPtBins; ip++) {
      h_m1d[ip]->Write();	
    }
    fout->Close();
  }
  else if (DO_FILL == 0) {
    TFile *fin = new TFile("files/mStore_MC.root");
    for(int ip = 0; ip < nPtBins; ip++) {
      fin->GetObject(Form("mH%.0f", ptBins[ip]), h_m1d[ip]);
      h_m1d[ip]->SetDirectory(0);
    }
    fin->Close();
  }
  
  cout << "all psi(2S) MC mass histograms filled" << endl << endl;

  // scale 1d histos and fill 2d histo
  for(int i = 0; i < nPtBins; i++) {
    //h_m1d[i]->Scale(1./h_m1d[i]->Integral());
    for(int j = 0; j < mbins; j++) {
      h_m2d->SetBinContent(j+1, i+1, h_m1d[i]->GetBinContent(j+1));
      h_m2d->SetBinError(j+1, i+1, h_m1d[i]->GetBinError(j+1));
    }
  }

  double fit_i = 2.92;
  double fit_f = 3.28;

  // define 2d function for fitting
  TF2 *f_cb = new TF2("f_cb", cb_func, fit_i, fit_f, 25, 100, 5*nPtBins+2, 2);
  // define constant parameters - f, mu
  f_cb->SetParName(0, "f");
  f_cb->SetParameter(0, 0.7);
  f_cb->SetParName(nPtBins+1, "mu");
  f_cb->SetParameter(nPtBins+1, 3.097);
  // define free parameters - N, alpha, n
  for(int i = 0; i < nPtBins; i++) {
    f_cb->SetParName(i+1, Form("N_%d", i));
    //f_cb->SetParameter(i+1, 5e-3);  
    f_cb->SetParameter(i+1, h_m1d[i]->Integral()/100.);  

    f_cb->SetParName(i+2+nPtBins, Form("sig1_%d", i));
    f_cb->SetParameter(i+2+nPtBins, 2.e-2);

    f_cb->SetParName(i+2+2*nPtBins, Form("sig2_%d", i));
    f_cb->SetParameter(i+2+2*nPtBins, 4.e-2);
    
    f_cb->SetParName(i+2+3*nPtBins, Form("n_%d", i));
    f_cb->SetParameter(i+2+3*nPtBins, 2.0);

    f_cb->SetParName(i+2+4*nPtBins, Form("alpha_%d", i));
    f_cb->SetParameter(i+2+4*nPtBins, 2.);
  }
  // tf1 for plotting in the 1D bins
  TF1 *f_1d = new TF1("f_1d", "[0]*cb_exp(x,[1],[3],[2],[5],[6]) + (1.-[0]) * cb_exp(x,[1],[4],[2],[5],[6])", fit_i, fit_f);
  f_1d->SetParNames("f", "N", "mu", "sigma1", "sigma2", "n", "alpha");
  
  // separate parts of the fit function
  TF1 *fp1 = new TF1("fp1", "[0]*cb_exp(x,[1],[3],[2],[4],[5])", fit_i, fit_f);
  fp1->SetParNames("f", "N", "mu", "sigma1", "n", "alpha");
  TF1 *fp2 = new TF1("fp2", "(1.-[0]) * cb_exp(x,[1],[3],[2],[4],[5])", fit_i, fit_f);
  fp2->SetParNames("f", "N", "mu", "sigma2", "n", "alpha");

  // fit the 2d function to the mass:pT map
  TCanvas *c = new TCanvas("", "", 700, 700);

  TFile *fout = new TFile("files/mfit_MC_mix.root", "recreate");

  for(int ifit = 0; ifit < 1; ifit++)
    {
      if(ifit == 1) { // fix n to its central value 1.2
	// define constant parameters - f, mu
	f_cb->SetParameter(0, 0.7);
	f_cb->SetParameter(nPtBins+1, 3.097);
	// define linear parameters - sigma1, sigma2
	f_cb->SetParameter(nPtBins+2, 8.e-5);
	f_cb->SetParameter(nPtBins+3, 2.e-2);
	f_cb->SetParameter(nPtBins+4, 1.e-4);
	f_cb->SetParameter(nPtBins+5, 3.e-2);
	// define free parameters - N, alpha, n
	for(int i = 0; i < nPtBins; i++) {
	  f_cb->SetParameter(i+1, h_m1d[i]->GetMaximum()/20.);  
	  f_cb->FixParameter(nPtBins+6+i, 1.2);
	  f_cb->SetParameter(i+6+2*nPtBins, 2.);
	}
      }
      if(ifit == 2) { // fix alpha to its central value 2.15
	// define constant parameters - f, mu
	f_cb->SetParameter(0, 0.7);
	f_cb->SetParameter(nPtBins+1, 3.097);
	// define linear parameters - sigma1, sigma2
	f_cb->SetParameter(nPtBins+2, 8.e-5);
	f_cb->SetParameter(nPtBins+3, 2.e-2);
	f_cb->SetParameter(nPtBins+4, 1.e-4);
	f_cb->SetParameter(nPtBins+5, 3.e-2);
	// define free parameters - N, alpha, n
	for(int i = 0; i < nPtBins; i++) {
	  f_cb->SetParameter(i+1, h_m1d[i]->GetMaximum()/20.);  
	  f_cb->FixParameter(i+6+2*nPtBins, 2.15);
	}
      }

      h_m2d->Fit("f_cb", "R");
  
      double pt_val[nPtBins], pt_err[nPtBins];
      double pars[5][nPtBins], epars[5][nPtBins];
      
      // cycle over all pT bins
      for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
	pt_val[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
	pt_err[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);

	// storing free parameters
	pars[0][i_pt] = f_cb->GetParameter(i_pt+1);
	pars[1][i_pt] = f_cb->GetParameter(i_pt+2+nPtBins);
	pars[2][i_pt] = f_cb->GetParameter(i_pt+2+2*nPtBins);
	pars[3][i_pt] = f_cb->GetParameter(i_pt+2+3*nPtBins);
	pars[4][i_pt] = f_cb->GetParameter(i_pt+2+4*nPtBins);
	epars[0][i_pt] = f_cb->GetParError(i_pt+1);
	epars[1][i_pt] = f_cb->GetParError(i_pt+2+nPtBins);
	epars[2][i_pt] = f_cb->GetParError(i_pt+2+2*nPtBins);
	epars[3][i_pt] = f_cb->GetParError(i_pt+2+3*nPtBins);
	epars[4][i_pt] = f_cb->GetParError(i_pt+2+4*nPtBins);


	if(ifit == 2) { // only plot in the last fit
	  c->SetLogy();
	  
	  h_m1d[i_pt]->SetMaximum(h_m1d[i_pt]->GetMaximum()*1.2);
	  h_m1d[i_pt]->SetMinimum(h_m1d[i_pt]->GetMaximum()*8e-3);
	  h_m1d[i_pt]->SetMarkerStyle(20);
	  h_m1d[i_pt]->SetMarkerSize(0.75);
	  h_m1d[i_pt]->SetMarkerColor(kBlack);
	  h_m1d[i_pt]->Draw("error");
	  f_1d->SetLineColor(kBlue);
	  f_1d->Draw("lsame");

	  // initializing f_1d and plotting
	  f_1d->SetParameters(f_cb->GetParameter(0),
			      pars[0][i_pt],
			      f_cb->GetParameter(nPtBins+1),
			      pars[1][i_pt],
			      pars[2][i_pt],
			      pars[3][i_pt],
			      pars[4][i_pt]);

	  fp1->SetParameters(f_cb->GetParameter(0), pars[0][i_pt], f_cb->GetParameter(nPtBins+1), pars[1][i_pt], pars[3][i_pt], pars[4][i_pt]);
	  fp1->SetLineColor(kRed);
	  fp1->SetLineStyle(kDashed);
	  fp1->Draw("lsame");
	  fp2->SetParameters(f_cb->GetParameter(0), pars[0][i_pt], f_cb->GetParameter(nPtBins+1), pars[2][i_pt], pars[3][i_pt], pars[4][i_pt]);
	  fp2->SetLineColor(kGreen);
	  fp2->SetLineStyle(kDashed);
	  fp2->Draw("lsame");
	  
	  c->SaveAs(Form("plots/MCMass/CB_pt%d_mix.pdf", i_pt));
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
	  
	  // plotting the puls
	  TH1F *fl = c->DrawFrame(lowm, -6, him, 6);
	  fl->GetYaxis()->SetTitleOffset(1.3);
	  fl->GetYaxis()->SetLabelOffset(0.01);
	  fl->SetTitle(Form("MC J/#psi Pulls (%.0f < p_{T} < %.0f)", ptBins[i_pt], ptBins[i_pt+1]));
	  
	  TGraph *g_pull = new TGraph(mbins, mv, pv);
	  g_pull->SetLineColor(kBlack);
	  g_pull->SetMarkerColor(kBlack);
	  g_pull->SetMarkerStyle(20);
	  g_pull->Draw("p");
	  
	  TLine *zero = new TLine(lowm, 0, him, 0);
	  zero->SetLineStyle(kDashed);
	  zero->Draw();
	  
	  c->SaveAs(Form("plots/MCMass/pulls_pt%d_mix.pdf", i_pt));
	  c->Clear();
	}
	pars[0][i_pt] /= (ptBins[i_pt+1]-ptBins[i_pt]);
	epars[0][i_pt] /= (ptBins[i_pt+1]-ptBins[i_pt]);
	pars[1][i_pt] *= 1e3;
	epars[1][i_pt] *= 1e3;
	pars[2][i_pt] *= 1e3;
	epars[2][i_pt] *= 1e3;	  
      }

      // plotting the free parameters
      double parmin[] = {0, 4.6e-3, 3.09, 0.018, 0.03,  0.6, 1.9};
      double parmax[] = {1, 5e-3,   3.10, 0.03,  0.055, 2.0, 2.4};
      string partit[] = {"f", "N", "#mu", "#sigma_{1}", "#sigma_{2}", "n", "#alpha"};
      string parax[] = {"f", "N", "#mu (GeV)", "#sigma_{1} (GeV)", "#sigma_{2} (GeV)", "n", "#alpha"};
      string parlab[] = {"f", "N", "mu", "sig1", "sig2", "n", "alpha"};

      int ct_p = 0, ct_pos = 0;
      int p_pos[] = {0, nPtBins+1};
  
      for(int i_p = 0; i_p < 7; i_p++) {
   
	TH1F *fl = c->DrawFrame(ptBins[0]-5, parmin[i_p], ptBins[nPtBins]+5, parmax[i_p]);
	fl->SetXTitle("p_{T} (GeV)");
	fl->SetYTitle(parax[i_p].c_str());
	fl->GetYaxis()->SetTitleOffset(1.3);
	fl->GetYaxis()->SetLabelOffset(0.01);
	fl->SetTitle(partit[i_p].c_str());

	// free parameters: N, n, alpha
	if(i_p != 0 && i_p != 2) {
	  TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[ct_p], pt_err, epars[ct_p]);
	  g_par->SetLineColor(kBlue);
	  g_par->SetMarkerColor(kBlue);
	  g_par->Draw("p");
	  g_par->Write(Form("%s_mix", parlab[i_p].c_str()));
	  ct_p++;
	}
	// constant parameters - mu and f
	else {
	  double par_val[nPtBins], par_err[nPtBins], mult = 1.;
	  if(i_p == 0) mult = 100.;
	  else mult = 1e3;
	  for(int i = 0; i < nPtBins; i++) {
	    par_val[i] = f_cb->GetParameter(p_pos[ct_pos])*mult;
	    par_err[i] = f_cb->GetParError(p_pos[ct_pos])*mult;
	  }
	  TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, par_val, pt_err, par_err);
	  g_par->SetLineColor(kBlue);
	  g_par->SetMarkerColor(kBlue);
	  g_par->Draw("p");
	  g_par->Write(Form("%s_mix", parlab[i_p].c_str()));
	  ct_pos++;
	}
    
	c->Clear();
      }
      TLine *l_chi = new TLine(ptBins[0], f_cb->GetChisquare()/f_cb->GetNDF(), ptBins[nPtBins], f_cb->GetChisquare()/f_cb->GetNDF());
      l_chi->Write(Form("chiN_mix"));
      
      ofstream ftex;
      ftex.open(Form("text_output/mfit_MC_mix.tex"));
      ftex << "\\begin{tabular}{c||c|c|c|c|c}\n";
      ftex << "$\\pt$ (GeV) & $N$ & $\\sigma_1$ (MeV) & $\\sigma_2$ (MeV) & $n$ & $\\alpha$ \\\\\n";
      ftex << "\\hline\n";

      for(int i = 0; i < nPtBins; i++) {
	// pT bin
	ftex << Form("$[%.0f, %.0f]$", ptBins[i], ptBins[i+1]);
	for(int i_p = 0; i_p < 5; i_p++) {
	  double mult = 1.;
	  
	  double val = pars[i_p][i]*mult, unc = epars[i_p][i]*mult;
	  int p_norm = 1.; 
	  if(unc != 0) 
	    p_norm = ceil(-log10(unc))+1;	
	  ftex << " & " <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc;
	}
	ftex <<  "\\\\\n";
      }
      ftex << "\\end{tabular}\n";
      ftex.close();

      ofstream fout2;
      fout2.open(Form("text_output/mfit_MC_mixA.tex"));
      fout2 << "\\begin{tabular}{c|c||c}\n";
      fout2 << " $f$ $(\\%)$ & $\\mu$ $(MeV)$ & $\\chi^2/$ndf \\\\\n";
      fout2 << "\\hline\n";

      // f
      double val = f_cb->GetParameter(p_pos[0])*100.;
      double unc = f_cb->GetParError(p_pos[0])*100.;
      int p_norm = ceil(-log10(unc))+1;	
      fout2 <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
      // mu
      val = f_cb->GetParameter(p_pos[1])*1000.;
      unc = f_cb->GetParError(p_pos[1])*1000.;
      p_norm = ceil(-log10(unc))+1;	
      fout2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
      // chi^2
      fout2 << setprecision(0) << f_cb->GetChisquare() << "/" << f_cb->GetNDF() << "\\\\\n";
      fout2 << "\\end{tabular}\n";
      fout2.close();


      cout << f_cb->GetChisquare() << "/" << f_cb->GetNDF() << endl;
    }
  fout->Close();
  
  c->Destructor(); 
}
