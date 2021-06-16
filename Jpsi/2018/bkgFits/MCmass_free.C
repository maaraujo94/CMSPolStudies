int DO_FILL = 1;
double gPI = TMath::Pi();

// crystal ball function
double cb_func(double m, double N, double sig, double m0, double n, double alpha)
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

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


// MAIN
void MCmass_free()
{
  // PART 1 : FILLING THE MASS HISTOS
  // prepare binning and histograms for plots
  const int nPtBins = 18;
  double ptBins[nPtBins+1];
  for(int i = 0; i < 1; i++) ptBins[i] = i+25.;
  for(int i = 0; i < 7; i++) ptBins[i+1] = 26.+i*2.;
  for(int i = 0; i < 2; i++) ptBins[i+8] = 40.+i*3;
  for(int i = 0; i < 1; i++) ptBins[i+10] = 46.+i*4;
  for(int i = 0; i < 4; i++) ptBins[i+11] = 50.+i*5.;
  for(int i = 0; i < 4; i++) ptBins[i+15] = 70.+i*10.;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // prepare mass histograms
  int mbins = 80;
  double lowm = 2.9, him = 3.3;
  TH1D **mHist = new TH1D*[nPtBins];
  for(int ip = 0; ip < nPtBins; ip++) {
    mHist[ip] = new TH1D(Form("mH%.0f", ptBins[ip]), Form("2018 MC M(#mu#mu) (%.0f < p_{T} < %.0f GeV)", ptBins[ip], ptBins[ip+1]), mbins, lowm, him);
  }
 
  cout << "all psi(2S) MC mass histograms initialized" << endl;

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
	      mHist[i_p]->Fill(mc_m);
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
	      mHist[i_p]->Fill(mc_m);
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
	      mHist[i_p]->Fill(mc_m);
	}
      }
    fin3->Close();

    TFile *fout = new TFile("files/mStore_MC.root", "recreate");
    for(int ip = 0; ip < nPtBins; ip++) {
      mHist[ip]->Write();	
    }
    fout->Close();
  }
  else if (DO_FILL == 0) {
    TFile *fin = new TFile("files/mStore_MC.root");
    for(int ip = 0; ip < nPtBins; ip++) {
      fin->GetObject(Form("mH%.0f", ptBins[ip]), mHist[ip]);
      mHist[ip]->SetDirectory(0);
    }
    fin->Close();
  }

  cout << "all psi(2S) MC mass histograms filled" << endl << endl;

  // scale histos for easier normalization
  for(int i = 0; i < nPtBins; i++)
    mHist[i]->Scale(1./mHist[i]->Integral());

  double fit_i = 2.92;
  double fit_f = 3.28;

  // PART 2 : DOING THE FIT 
  // define the fit function
  TF1 *f_cb = new TF1("f_cb", "[1]*cb_func(x,[0],[3],[2],[5],[6]) + (1.-[1]) * cb_func(x,[0],[4],[2],[5],[6])", fit_i, fit_f);
  f_cb->SetParNames("N", "f", "mu", "sigma1", "sigma2", "n", "alpha");
  double alpha= 2.1, n = 1.2;
  double mu = 3.097, sig1 = 0.025, sig2 = 0.04;
  double f = 0.7;

  TH1D *h_pull = new TH1D("h_pull", "MC J/#psi Pulls", mbins, lowm, him);

  TCanvas *c = new TCanvas("", "", 700, 700);

  // output for latex
  ofstream fout;
  fout.open(Form("text_output/mfit_MC_free.tex"));
  fout << "\\begin{tabular}{c||c|c|c|c|c|c|c||c}\n";
  fout << "$\\pt$ (GeV) & $N$ $(\\times1e5)$ & $f$ (\\%) & $\\mu$ (MeV) & $\\sigma_1$ (MeV) & $\\sigma_2$ (MeV) & $n$ &  $\\alpha$ & $\\chi^2$/ndf \\\\\n";
  fout << "\\hline\n";

  double pars[8][nPtBins], epars[8][nPtBins];
  double pt_val[nPtBins], pt_err[nPtBins];

  // cycle fitting over all pT bins
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    pt_val[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
    pt_err[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);
    
    double N = mHist[i_pt]->GetMaximum()/10.;
    f_cb->SetParameters(N, f, mu, sig1, sig2, n, alpha);
    
    // fitting, plotting histo
    c->SetLogy();
    mHist[i_pt]->Fit("f_cb", "R");
    mHist[i_pt]->Draw("error");

    // store the fit parameters in output latex table
    fout << Form("$[%.0f, %.0f]$", ptBins[i_pt], ptBins[i_pt+1]) << " & ";
    for(int i_p = 0; i_p < 7; i_p++) {
      pars[i_p][i_pt] = f_cb->GetParameter(i_p);
      epars[i_p][i_pt] = f_cb->GetParError(i_p);
      if(f_cb->GetParError(i_p) != 0) {
	double mult = 1.;
	if(i_p == 0) mult = 1e5;
	else if(i_p == 1) mult = 100.;
	else if((i_p > 1 && i_p < 5)) mult = 1000.;
	int p_norm = ceil(-log10(f_cb->GetParError(i_p)*mult))+1;	
	fout <<  setprecision(p_norm) << fixed << f_cb->GetParameter(i_p)*mult << " $\\pm$ " << f_cb->GetParError(i_p)*mult << " & " << endl;
      }
      else {
	int p_norm = 3;//ceil(-log10(f_cb->GetParameter(i_p)))+1;	
	fout << setprecision(p_norm) << fixed << f_cb->GetParameter(i_p) << " & " << endl;
      }
    }
    pars[7][i_pt] = f_cb->GetChisquare()/f_cb->GetNDF();//TMath::Prob(f_cb->GetChisquare(), f_cb->GetNDF());
    epars[7][i_pt] = 0;
    fout << Form("%.0f", f_cb->GetChisquare()) << "/" << Form("%d", f_cb->GetNDF()) << "\\\\\n";
    
    c->SaveAs(Form("plots/MCMass/CB_pt%d_free.pdf", i_pt));
    c->Clear();

    // get the pulls distribution
    double mv[mbins], pv[mbins];
    for(int i_m = 0 ; i_m < mbins; i_m++) {
      mv[i_m] = mHist[i_pt]->GetBinCenter(i_m+1);
      double fitv = f_cb->Eval(mv[i_m]);
      double datav = mHist[i_pt]->GetBinContent(i_m+1);
      double datau = mHist[i_pt]->GetBinError(i_m+1);
      if(datau > 0 && mv[i_m] > fit_i && mv[i_m] < fit_f) {
	pv[i_m] = (fitv-datav)/datau;
      }
      else
	pv[i_m] = 0;
    }

    c->SetLogy(0);

    // plot the pulls
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
    
    c->SaveAs(Form("plots/MCMass/pulls_pt%d_free.pdf", i_pt));
    c->Clear();
    
  }
  fout << "\\end{tabular}\n";
  fout.close();

  // PART 3 : PLOT FIT PARAMS
  // define aux arrays
  double parmin[] = {4.7e-3, 0, 3.09, 0.018,  0.03,  0.6, 1.9};
  double parmax[] = {5e-3,   1, 3.1,  0.03,   0.055, 2.0, 2.4};
  string partit[] = {"N", "f", "#mu", "#sigma_{1}", "#sigma_{2}", "n", "#alpha"};
  string parax[] = {"N", "f", "#mu (GeV)", "#sigma_{1} (GeV)", "#sigma_{2} (GeV)", "n", "#alpha"};
  string parlab[] = {"N", "f", "mu", "sig1", "sig2", "n", "alpha"};
  
  TFile *fout2 = new TFile("files/mfit_MC_free.root", "recreate");
  
  // cycle for plotting over all fit parameters
  for(int i_p = 0; i_p < 7; i_p++) {
    TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[i_p], pt_err, epars[i_p]);

    TH1F *fl = c->DrawFrame(ptBins[0]-5, parmin[i_p], ptBins[nPtBins]+5, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(partit[i_p].c_str());
    
    g_par->SetLineColor(kBlack);
    g_par->SetMarkerColor(kBlack);
    g_par->Draw("p");

    g_par->Write(Form("%s_free", parlab[i_p].c_str()));

    c->Clear();
  }
  TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[7], pt_err, epars[7]);

  TH1F *fl = c->DrawFrame(ptBins[0]-5, 0, ptBins[nPtBins]+5, 10);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#chi^{2}/ndf");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("#chi^{2} / ndf");
    
  g_par->SetLineColor(kBlack);
  g_par->SetMarkerColor(kBlack);
  g_par->SetMarkerStyle(20);
  g_par->Draw("p");

  g_par->Write("chiN_free");
      
  c->Clear();
  
  c->Destructor();
  fout2->Close();
}
