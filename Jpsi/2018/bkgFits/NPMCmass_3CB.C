double gPI = TMath::Pi();
int DO_FILL = 1;

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
// parameters: N, f1, mu, sig1, sig2, n, alpha, fG, sigG, f3, sig3
double cb_func(double *x, double *par)
{
  // get m, pt and corresp pt bin
  double m = x[0];

  double N = par[0];
  double f1 = par[1]; // f is constant in pt, only take the first value and use in all cases

  double mu = par[2]; // mu is constant in pt
  double sig1 = par[3]; 
  double sig2 = par[4]; // sigmas linear in pt
  
  double n = par[5]; // n is constant in pt
  double alpha = par[6]; // alpha is constant in pt

  double fG = par[7];
  double sigG = par[8];

  double f3 = par[9];
  double sig3 = par[10];

  double func = f1 * cb_exp(m, N, sig1, mu, n, alpha) + (1.-f1-f3-fG) * cb_exp(m, N, sig2, mu, n, alpha) + f3 * cb_exp(m, N, sig3, mu, n, alpha) + fG * g_exp(m, N, sigG, mu);
  return func;
}

// MAIN
void NPMCmass_3CB()
{
  // prepare mass histograms
  // Full lt range and [100, 500] micron interval only
  int mbins = 40;
  double lowm = 2.9, him = 3.3;
  TH1D *h_mF = new TH1D("mH_Full","2018 NP MC M(#mu#mu)", mbins, lowm, him);
  TH1D *h_mNP = new TH1D("mH_NP","2018 NP MC M(#mu#mu)", mbins, lowm, him);
  // same at high pT
  TH1D *h_mF_hp = new TH1D("mH_Full_hp","2018 NP MC M(#mu#mu), p_{T} > 50 GeV", mbins, lowm, him);
  TH1D *h_mNP_hp = new TH1D("mH_NP_hp","2018 NP MC M(#mu#mu), p_{T} > 50 GeV", mbins, lowm, him);
  
  if(DO_FILL == 1) {
    // filling the histo
    // open and read the data tree
    TFile *fin1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCNP18_cos.root");
    TTree *tree1 = (TTree*)fin1->Get("MC_cos");

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
	if(mc_pt > 25 && mc_pt < 120 && abs(mc_y) < 1.2) {
	  h_mF->Fill(mc_m);
	  if(mc_pt > 50)
	    h_mF_hp->Fill(mc_m);
	  if(mc_lt > 0.01 && mc_lt < 0.05) {
	    h_mNP->Fill(mc_m);
	    if(mc_pt > 50)
	      h_mNP_hp->Fill(mc_m);
	  }
	}
      }
    fin1->Close();

    TFile *fout = new TFile("files/mStore_MCNP.root", "recreate");
    h_mF->Write();	
    h_mNP->Write();	
    h_mF_hp->Write();	
    h_mNP_hp->Write();	
    
    fout->Close();
  }
  else if (DO_FILL == 0) {
    TFile *fin = new TFile("files/mStore_MCNP.root");
    fin->GetObject("mH_Full", h_mF);
    h_mF->SetDirectory(0);
    fin->GetObject("mH_NP", h_mNP);
    h_mNP->SetDirectory(0);
    fin->GetObject("mH_Full_hp", h_mF_hp);
    h_mF_hp->SetDirectory(0);
    fin->GetObject("mH_NP_hp", h_mNP_hp);
    h_mNP_hp->SetDirectory(0);
    
    fin->Close();
  }
  
  cout << "all MC mass histograms filled" << endl << endl;

  double fit_i = 2.92;
  double fit_f = 3.28;

  // define function for fitting
  TF1 *f_cb = new TF1("f_cb", cb_func, fit_i, fit_f, 11, 1);

  string par_n[] = {"N", "f1", "mu", "sig1", "sig2", "n", "alpha", "fG", "sigG", "f3", "sig3"};
  double par_v[] = {1., 0.4, 3.097, 2e-2, 3e-2, 2., 2., 1e-2, 8e-2, 0.1, 4e-2};
  // define parameters - all free
  for(int i = 0; i < 11; i++) {
    f_cb->SetParName(i, par_n[i].c_str());
  }

  // separate parts of the fit function
  TF1 *fcb1 = new TF1("fcb1", "[1]*cb_exp(x,[0],[3],[2],[4],[5])", fit_i, fit_f);
  fcb1->SetParNames("N", "f1", "mu", "sigma1", "n", "alpha");
  TF1 *fcb2 = new TF1("fcb2", "(1.-[1]-[6]-[7]) * cb_exp(x,[0],[3],[2],[4],[5])", fit_i, fit_f);
  fcb2->SetParNames("N", "f1", "mu", "sigma2", "n", "alpha", "fG", "f3");
  TF1 *fcb3 = new TF1("fcb3", "[1]*cb_exp(x,[0],[3],[2],[4],[5])", fit_i, fit_f);
  fcb3->SetParNames("N", "f3", "mu", "sigma3", "n", "alpha");
  TF1 *fg = new TF1("fg", "[1]*g_exp(x,[0],[3],[2])", fit_i, fit_f);
  fg->SetParNames("N", "fG", "mu", "sigmaG");

  // run and plot the fit
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLeftMargin(0.12);
  
  TFile *fout = new TFile("files/MCNPfit.root", "recreate");

  // first running the full fit
  f_cb->SetParameters(par_v);
  f_cb->SetParameter(0, h_mF->GetMaximum()/20.);
  f_cb->SetLineColor(kBlue);
  f_cb->SetNpx(1000);
  h_mF->Fit("f_cb", "R");

  // storing parameters
  double pars[11], epars[11];
  for(int j = 0; j < 11; j++) {
    pars[j] = f_cb->GetParameter(j);
    epars[j] = f_cb->GetParError(j);
  }

  h_mF->SetMaximum(h_mF->GetMaximum()*1.1);
  h_mF->SetMinimum(0);
  h_mF->SetStats(0);
  h_mF->GetYaxis()->SetTitle(Form("Events per %.0f MeV", (him-lowm)/mbins*1000));
  h_mF->GetYaxis()->SetTitleOffset(1.8);
  h_mF->GetXaxis()->SetTitle(Form("M(#mu#mu) (GeV)"));
  h_mF->SetMarkerStyle(20);
  h_mF->SetMarkerSize(0.75);
  h_mF->SetMarkerColor(kBlack);
  h_mF->Draw("error");

  // separate parts of the fit function
  fcb1->SetParameters(pars[0], pars[1], pars[2], pars[3], pars[5], pars[6]);
  fcb1->SetLineColor(kRed);
  fcb1->SetLineStyle(kDashed);
  fcb1->Draw("lsame");
  fcb2->SetParameters(pars[0], pars[1], pars[2], pars[4], pars[5], pars[6], pars[7], pars[9]);
  fcb2->SetLineColor(kGreen);
  fcb2->SetLineStyle(kDashed);
  fcb2->Draw("lsame");
  fcb3->SetParameters(pars[0], pars[9], pars[2], pars[10], pars[5], pars[6]);
  fcb3->SetLineColor(kBlack);
  fcb3->SetLineStyle(kDashed);
  fcb3->Draw("lsame");
  fg->SetParameters(pars[0], pars[7], pars[2], pars[8]);
  fg->SetLineColor(kViolet);
  fg->SetLineStyle(kDashed);
  fg->Draw("lsame");

  TLatex lc;
  lc.SetTextSize(0.04);
  lc.DrawLatex(3.15, 700, Form("Free n=%.2f", pars[5]));

  c->SaveAs("plots/NPMCmass/NP3_fit.pdf");
  c->Clear();
  
  // calculating pulls
  double mv[mbins], pv[mbins];
  for(int i_m = 0 ; i_m < mbins; i_m++) {
    mv[i_m] = h_mF->GetBinCenter(i_m+1);
    double fitv = f_cb->Eval(mv[i_m]);
    double datav = h_mF->GetBinContent(i_m+1);
    double datau = h_mF->GetBinError(i_m+1);
    if(datau > 0 && mv[i_m] > fit_i && mv[i_m] < fit_f) {
      pv[i_m] = (datav-fitv)/datau;
    }
    else
      pv[i_m] = 0;
  }
  
  c->SetLogy(0);

  double lowmp = fit_i, himp = fit_f;
  
  // plotting the pulls
  TH1F *fl = c->DrawFrame(lowmp, -7, himp, 7);
  fl->SetXTitle("M(#mu#mu) (GeV)");
  fl->SetYTitle("pulls");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("MC mass fit pulls");
  
  TGraph *g_pull = new TGraph(mbins, mv, pv);
  g_pull->SetLineColor(kBlack);
  g_pull->SetMarkerColor(kBlack);
  g_pull->SetMarkerStyle(20);
  g_pull->Draw("p");
  
  TLine *zero = new TLine(lowmp, 0, himp, 0);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLine *plim1 = new TLine(lowmp, -5, himp, -5);
  plim1->SetLineStyle(kDotted);
  plim1->Draw("lsame");
  TLine *plim2 = new TLine(lowmp, -3, himp, -3);
  plim2->SetLineStyle(kDotted);
  plim2->Draw("lsame");
  TLine *plim3 = new TLine(lowmp, 3, himp, 3);
  plim3->SetLineStyle(kDotted);
  plim3->Draw("lsame");
  TLine *plim4 = new TLine(lowmp, 5, himp, 5);
  plim4->SetLineStyle(kDotted);
  plim4->Draw("lsame");
	  
  c->SaveAs("plots/NPMCmass/NP3_pulls.pdf");
  c->Clear();

  ofstream ftable;
  ftable.open("text_output/NP_fit.txt" );
  for(int i = 0; i < 11; i++) {
    if(i == 0)
      ftable << pars[i]/(120.-25.) << "\t " << epars[i]/(120.-25.) << "\t ";
    else
      ftable << pars[i] << "\t " << epars[i] << "\t ";
  }
  ftable << f_cb->GetChisquare() << "\t " << f_cb->GetNDF() << "\n";
  ftable.close();

  
  // then running just the NP range fit
  f_cb->SetParameters(par_v);
  f_cb->SetParameter(0, h_mNP->GetMaximum()/20.);
  h_mNP->Fit("f_cb", "R");

  // storing parameters
  for(int j = 0; j < 11; j++) {
    pars[j] = f_cb->GetParameter(j);
    epars[j] = f_cb->GetParError(j);
  }
  
  h_mNP->SetMaximum(h_mNP->GetMaximum()*1.1);
  h_mNP->SetMinimum(0);
  h_mNP->SetStats(0);
  h_mNP->GetYaxis()->SetTitle(Form("Events per %.0f MeV", (him-lowm)/mbins*1000));
  h_mNP->GetYaxis()->SetTitleOffset(1.8);
  h_mNP->GetXaxis()->SetTitle(Form("M(#mu#mu) (GeV)"));
  h_mNP->SetMarkerStyle(20);
  h_mNP->SetMarkerSize(0.75);
  h_mNP->SetMarkerColor(kBlack);
  h_mNP->Draw("error");

  // separate parts of the fit function
  fcb1->SetParameters(pars[0], pars[1], pars[2], pars[3], pars[5], pars[6]);
  fcb1->SetLineColor(kRed);
  fcb1->SetLineStyle(kDashed);
  fcb1->Draw("lsame");
  fcb2->SetParameters(pars[0], pars[1], pars[2], pars[4], pars[5], pars[6], pars[7], pars[9]);
  fcb2->SetLineColor(kGreen);
  fcb2->SetLineStyle(kDashed);
  fcb2->Draw("lsame");
  fcb3->SetParameters(pars[0], pars[9], pars[2], pars[10], pars[5], pars[6]);
  fcb3->SetLineColor(kBlack);
  fcb3->SetLineStyle(kDashed);
  fcb3->Draw("lsame");
  fg->SetParameters(pars[0], pars[7], pars[2], pars[8]);
  fg->SetLineColor(kViolet);
  fg->SetLineStyle(kDashed);
  fg->Draw("lsame");

  TLatex lcR;
  lcR.SetTextSize(0.04);
  lcR.DrawLatex(3.15, 6000, Form("Free n=%.2f", pars[5]));
  lcR.DrawLatex(3.15, 5500, Form("#chi^{2}/ndf = %.0f/%d", f_cb->GetChisquare(), f_cb->GetNDF()));

  c->SaveAs("plots/NPMCmass/NPR3_fit.pdf");
  c->Clear();
  
  // calculating pulls
  for(int i_m = 0 ; i_m < mbins; i_m++) {
    mv[i_m] = h_mNP->GetBinCenter(i_m+1);
    double fitv = f_cb->Eval(mv[i_m]);
    double datav = h_mNP->GetBinContent(i_m+1);
    double datau = h_mNP->GetBinError(i_m+1);
    if(datau > 0 && mv[i_m] > fit_i && mv[i_m] < fit_f) {
      pv[i_m] = (datav-fitv)/datau;
    }
    else
      pv[i_m] = 0;
  }
  
  c->SetLogy(0);

  // plotting the pulls
  TH1F *flNP = c->DrawFrame(lowmp, -7, himp, 7);
  flNP->SetXTitle("M(#mu#mu) (GeV)");
  flNP->SetYTitle("pulls");
  flNP->GetYaxis()->SetTitleOffset(1.3);
  flNP->GetYaxis()->SetLabelOffset(0.01);
  flNP->SetTitle("MC mass fit pulls");
  
  TGraph *g_pullNP = new TGraph(mbins, mv, pv);
  g_pullNP->SetLineColor(kBlack);
  g_pullNP->SetMarkerColor(kBlack);
  g_pullNP->SetMarkerStyle(20);
  g_pullNP->Draw("p");
  
  zero->Draw();

  plim1->Draw("lsame");
  plim2->Draw("lsame");
  plim3->Draw("lsame");
  plim4->Draw("lsame");

  c->SaveAs("plots/NPMCmass/NPR3_pulls.pdf");
  c->Clear();
      
  c->Destructor();

  ofstream ftableNP;
  ftableNP.open("text_output/NPR_fit.txt" );
  for(int i = 0; i < 11; i++) {
    if(i == 0)
      ftableNP << pars[i]/(120.-25.) << "\t " << epars[i]/(120.-25.) << "\t ";
    else
      ftableNP << pars[i] << "\t " << epars[i] << "\t ";
  }
  ftableNP << f_cb->GetChisquare() << "\t " << f_cb->GetNDF() << "\n";
  ftableNP.close();

}
