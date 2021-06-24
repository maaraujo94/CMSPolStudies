//pt bins defined globally for access from functions
const int nPtBins = 7;
double ptBins[nPtBins+1];

// functions to access within other functions
TF1 *fres;
TF1 *fNP;
TF1 *fbkL, *fbkR, *fbkD;

// define negative exponential only for positive x
double pos_exp(double x, double ld)
{
  if(x > 0) return exp(-x/ld);
  else return 0;
}

// define final fit function summing the PR and NP contributions
// parameters: N_PR (per bin), N_NP (per bin), f (constant), mu (constant), sigma1 (per bin), sigma2 (per bin), lambda (per bin)
double func_sum(double *xx, double *pp)
{
  // get pt bin
  int pt_bin;
  for(int i = 0; i < nPtBins; i++)
    if(ptBins[i] < xx[1] && ptBins[i+1] > xx[1])
      pt_bin = i;

  double f = pp[2*nPtBins];
  double mu = pp[2*nPtBins+1];

  double N_PR = pp[pt_bin];
  double N_NP = pp[nPtBins+pt_bin];
  double sig1 = pp[2+2*nPtBins+pt_bin];
  double sig2 = pp[2+3*nPtBins+pt_bin];
  double ld = pp[2+4*nPtBins+pt_bin];

  double inp[] = {xx[0]};
  double par_1[] = {N_PR, f, mu, sig1, sig2};
  double par_2[] = {N_NP, f, mu, sig1, sig2, ld};

  return fres->EvalPar(inp, par_1) + fNP->EvalPar(inp, par_2);
}

// define final fit function summing the PR and NP contributions
double sum_1d(double *xx, double *pp)
{
  double lt = xx[0];
  double N_PR = pp[0], N_NP = pp[1], f2 = pp[2], mu = pp[3], sig1 = pp[4], sig2 = pp[5], ld = pp[6];

  double inp[] = {lt};
  double par_1[] = {N_PR, f2, mu, sig1, sig2};
  double par_2[] = {N_NP, f2, mu, sig1, sig2, ld};
  
  return fres->EvalPar(inp, par_1) + fNP->EvalPar(inp, par_2);
}


void ltBkg()
{
  // prepare binning and histograms for plots
  for(int i=0; i<3; i++) ptBins[i] = 7.*i+25.;
  for(int i=0; i<4; i++) ptBins[i+3] = 46.+10.*i;
  ptBins[7] = 100;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  TH1D **h_d1d = new TH1D*[nPtBins];  
  TFile *fin = new TFile("files/ltStore.root");
  for(int ip = 0; ip < nPtBins; ip++) {
    fin->GetObject(Form("ltH%.0f", ptBins[ip]), h_d1d[ip]);
    h_d1d[ip]->SetDirectory(0);
  }
  fin->Close();

  int tbins = h_d1d[0]->GetNbinsX();
  double lowt = h_d1d[0]->GetXaxis()->GetBinLowEdge(1);
  double hit = h_d1d[0]->GetXaxis()->GetBinUpEdge(tbins);
  double wbin = (hit-lowt)/(double)tbins;

  // define aux vals for plotting
  double pr_lim = 0.1;
  double np_lim = 0.14;
  double lowPlot = -0.15;

  // Fill 2d histo
  TH2D *h_d2d = new TH2D("h_d2d", "2018 data c#tau", tbins, lowt, hit, nPtBins, ptBins);

  for(int i = 0; i < nPtBins; i++) {
    for(int j = 0; j < tbins; j++) {
      h_d2d->SetBinContent(j+1, i+1, h_d1d[i]->GetBinContent(j+1));
      h_d2d->SetBinError(j+1, i+1, h_d1d[i]->GetBinError(j+1));
    }
  }

  // define the resolution (=PR) function
  fres = new TF1("fres", "[0]*((1-[1])*TMath::Gaus(x, [2],[3]) + [1]*TMath::Gaus(x, [2], [4]))", 5*lowt, 5*hit);
  
  // define the NP function by convolution
  TF1 *fexp = new TF1("fexp", "pos_exp(x,[0])", 5*lowt, 5*hit);
  TF1Convolution *fcNP = new TF1Convolution(fres, fexp, 5*lowt, 5*hit);
  fcNP->SetRange(5*lowt, 5*hit);
  fcNP->SetNofPointsFFT(1000);
  fNP = new TF1("fNP", *fcNP, lowt, hit, fcNP->GetNpar());

  TF2 *fitS = new TF2("fitS", func_sum, lowt, hit, ptBins[0], ptBins[nPtBins], 2+5*nPtBins, 2);

  // define constant parameters - f, mu
  fitS->SetParName(2*nPtBins, "f");
  fitS->SetParameter(2*nPtBins, 0.2);
  fitS->FixParameter(2*nPtBins, 0.14);

  fitS->SetParName(2*nPtBins+1, "mu");
  fitS->SetParameter(2*nPtBins+1, 0);
  fitS->FixParameter(2*nPtBins+1, 0);
  // define free parameters - N, sigma, lambda
  for(int i = 0; i < nPtBins; i++) {
    fitS->SetParName(i, Form("N_PR_%d", i));
    fitS->SetParameter(i, h_d1d[i]->GetMaximum()/2.);
    
    fitS->SetParName(i+nPtBins, Form("N_NP_%d", i));
    fitS->SetParameter(i+nPtBins, h_d1d[i]->GetMaximum()*2.);
    
    fitS->SetParName(i+2*nPtBins+2, Form("sig1_%d", i));
    fitS->SetParameter(i+2*nPtBins+2, 1e-2);

    fitS->SetParName(i+3*nPtBins+2, Form("sig2_%d", i));
    fitS->SetParameter(i+3*nPtBins+2, 2e-2);

    fitS->SetParName(i+4*nPtBins+2, Form("lambda_%d", i));
    fitS->SetParameter(i+4*nPtBins+2, 0.4);
    //fitS->FixParameter(i+4*nPtBins+2, 0.355);
  }
  // fit the 2d function to the lifetime:pT map
  TCanvas *c = new TCanvas("", "", 700, 700);
  h_d2d->Fit("fitS", "V");

  // tf1 for plotting in the 1D bins
  // separate parts of the fit function - given by fres and fNP
  TF1 *f_1d = new TF1("f_1d", sum_1d, lowt, hit, 7);
  f_1d->SetParNames("N_PR", "N_NP", "f", "mu", "sigma1", "sigma2", "lambda");

  double pt_val[nPtBins], pt_err[nPtBins];
  double pars[5][nPtBins], epars[5][nPtBins];
  double fracNP[nPtBins], efz[nPtBins];

  // cycle over all pT bins
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    pt_val[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
    pt_err[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);

    // storing free parameters
    pars[0][i_pt] = fitS->GetParameter(i_pt);
    pars[1][i_pt] = fitS->GetParameter(i_pt+nPtBins);
    pars[2][i_pt] = fitS->GetParameter(i_pt+2+2*nPtBins);
    pars[3][i_pt] = fitS->GetParameter(i_pt+2+3*nPtBins);
    pars[4][i_pt] = fitS->GetParameter(i_pt+2+4*nPtBins);
    epars[0][i_pt] = fitS->GetParError(i_pt);
    epars[1][i_pt] = fitS->GetParError(i_pt+nPtBins);
    epars[2][i_pt] = fitS->GetParError(i_pt+2+2*nPtBins);
    epars[3][i_pt] = fitS->GetParError(i_pt+2+3*nPtBins);
    epars[4][i_pt] = fitS->GetParError(i_pt+2+4*nPtBins);

    // initializing f_1d and plotting
    f_1d->SetParameters(pars[0][i_pt],
			pars[1][i_pt],
			fitS->GetParameter(2*nPtBins),
			fitS->GetParameter(1+2*nPtBins),
			pars[2][i_pt],
			pars[3][i_pt],
			pars[4][i_pt]);


    c->SetLogy();
     
    h_d1d[i_pt]->SetMaximum(h_d1d[i_pt]->GetMaximum()*1.2);
    h_d1d[i_pt]->SetMinimum(h_d1d[i_pt]->GetMaximum()*5e-4);

    TH1F *fh = c->DrawFrame(lowPlot, h_d1d[i_pt]->GetMinimum(), hit, h_d1d[i_pt]->GetMaximum());
    fh->SetXTitle("c#tau (mm)");
    fh->SetYTitle(Form("Events / %.0f #mum", wbin*1000.));
    fh->GetYaxis()->SetTitleOffset(1.3);
    fh->GetYaxis()->SetLabelOffset(0.01);
    fh->SetTitle(h_d1d[i_pt]->GetTitle());

    h_d1d[i_pt]->SetMarkerStyle(20);
    h_d1d[i_pt]->SetMarkerSize(0.75);
    h_d1d[i_pt]->SetMarkerColor(kBlack);
    h_d1d[i_pt]->SetLineColor(kBlack);
    h_d1d[i_pt]->Draw("error same");
    
    f_1d->SetLineColor(kBlue);
    f_1d->Draw("lsame");

    // draw fit contributions
    fres->SetParameters(f_1d->GetParameter(0), f_1d->GetParameter(2), f_1d->GetParameter(3), f_1d->GetParameter(4), f_1d->GetParameter(5));
    fres->SetLineStyle(kDashDotted);
    fres->SetLineColor(kGreen);
    fres->Draw("lsame");
    fNP->SetParameters(f_1d->GetParameter(1), f_1d->GetParameter(2), f_1d->GetParameter(3), f_1d->GetParameter(4), f_1d->GetParameter(5), f_1d->GetParameter(6));
    fNP->SetLineStyle(kDashDotted);
    fNP->SetLineColor(kViolet);
    fNP->Draw("lsame");

    // aux lines for the 2.5 sigma and 4 sigma limits
    TLine *lsig1 = new TLine(-pr_lim, h_d1d[i_pt]->GetMinimum(), -pr_lim, h_d1d[i_pt]->GetMaximum());
    lsig1->SetLineStyle(kDashed);
    lsig1->Draw("lsame");
    TLine *lsig2 = new TLine(pr_lim, h_d1d[i_pt]->GetMinimum(), pr_lim, h_d1d[i_pt]->GetMaximum());
    lsig2->SetLineStyle(kDashed);
    lsig2->Draw("lsame");
    TLine *lsig3 = new TLine(np_lim, h_d1d[i_pt]->GetMinimum(), np_lim, h_d1d[i_pt]->GetMaximum());
    lsig3->SetLineStyle(kDashed);
    lsig3->Draw("lsame");

    c->SaveAs(Form("plots/lifetime/fit_pt%d.pdf", i_pt));
    c->Clear();

    // get the NP fraction in the signal region (+- 100 mum)
    double evt_NP = fNP->Integral(-pr_lim, pr_lim);
    
    double min_bin = h_d1d[i_pt]->GetXaxis()->FindBin(-pr_lim);
    double max_bin = h_d1d[i_pt]->GetXaxis()->FindBin(pr_lim);
    double evt_all = h_d1d[i_pt]->Integral(min_bin, max_bin, "width");
    fracNP[i_pt] = evt_NP / evt_all;
    efz[i_pt] = 0;

    // calculating pulls
    double tv[tbins], pv[tbins], dv[tbins];
    for(int i_t = 0 ; i_t < tbins; i_t++) {
      tv[i_t] = h_d1d[i_pt]->GetBinCenter(i_t+1);
      double fitv = f_1d->Eval(tv[i_t]);
      double datav = h_d1d[i_pt]->GetBinContent(i_t+1);
      double datau = h_d1d[i_pt]->GetBinError(i_t+1);
      if(datau > 0 ) pv[i_t] = (datav-fitv)/datau;
      else pv[i_t] = 0;
      if(fitv > 0 ) dv[i_t] = (datav-fitv)/fitv *100.;
      else dv[i_t] = 0;
    }
    
    c->SetLogy(0);
    
    // plotting the pulls

    TH1F *fp = c->DrawFrame(lowPlot, -15, hit, 15);
    fp->SetXTitle("c#tau (mm)");
    fp->SetYTitle("pulls");
    fp->GetYaxis()->SetTitleOffset(1.3);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("Lifetime fit pulls (%.0f < p_{T} < %.0f GeV)", ptBins[i_pt], ptBins[i_pt+1]));
  
    TGraph *g_pull = new TGraph(tbins, tv, pv);
    g_pull->SetLineColor(kBlack);
    g_pull->SetMarkerColor(kBlack);
    g_pull->SetMarkerStyle(20);
    g_pull->SetMarkerSize(.75);
    g_pull->Draw("psame");

    // aux lines - pull=0 and sigma limits
    TF1 *fcons = new TF1("fcons", "[0]", lowPlot, hit);
    fcons->SetParameter(0, 0);
    fcons->SetLineStyle(kDashed);
    fcons->SetLineColor(kBlack);
    fcons->Draw("lsame");

    TLine *psig1 = new TLine(-pr_lim, -15, -pr_lim, 15);
    psig1->SetLineStyle(kDashed);
    psig1->Draw("lsame");
    TLine *psig2 = new TLine(pr_lim, -15, pr_lim, 15);
    psig2->SetLineStyle(kDashed);
    psig2->Draw("lsame");
    TLine *psig3 = new TLine(np_lim, -15, np_lim, 15);
    psig3->SetLineStyle(kDashed);
    psig3->Draw("lsame");
  
    c->SaveAs(Form("plots/lifetime/pulls_pt%d.pdf", i_pt));
    c->Clear();

    // plotting the devs
    TH1F *fd = c->DrawFrame(lowPlot, -100, hit, 100);
    fd->SetXTitle("c#tau (mm)");
    fd->SetYTitle("relative difference (%)");
    fd->GetYaxis()->SetTitleOffset(1.3);
    fd->GetYaxis()->SetLabelOffset(0.01);
    fd->SetTitle(Form("Lifetime fit rel. difference (%.0f < p_{T} < %.0f GeV)", ptBins[i_pt], ptBins[i_pt+1]));

    TGraph *g_dev = new TGraph(tbins, tv, dv);
    g_dev->SetLineColor(kBlack);
    g_dev->SetMarkerColor(kBlack);
    g_dev->SetMarkerStyle(20);
    g_dev->SetMarkerSize(.75);
    g_dev->Draw("psame");
  
    // aux lines - pull = 0 and sigma limits
    fcons->Draw("lsame");

    TLine *dsig1 = new TLine(-pr_lim, -100, -pr_lim, 100);
    dsig1->SetLineStyle(kDashed);
    dsig1->Draw("lsame");
    TLine *dsig2 = new TLine(pr_lim, -100, pr_lim, 100);
    dsig2->SetLineStyle(kDashed);
    dsig2->Draw("lsame");
    TLine *dsig3 = new TLine(np_lim, -100, np_lim, 100);
    dsig3->SetLineStyle(kDashed);
    dsig3->Draw("lsame");

    c->SaveAs(Form("plots/lifetime/devs_pt%d.pdf", i_pt));
    c->Clear();
  }
  
  // storing the free parameters
  TFile *fout = new TFile("files/tfit.root", "recreate");
  string parlab[] = {"N_PR", "N_NP", "f", "mu", "sig1", "sig2", "lambda"};
  int ct_p = 0, ct_pos = 0;
  int p_pos[] = {2*nPtBins, 2*nPtBins+1};  

  for(int i_p = 0; i_p < 7; i_p++) {
    // free parameters: N, sig, lambda
    if(i_p < 2 || i_p > 3) {
      TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[ct_p], pt_err, epars[ct_p]);
      g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
      ct_p++;
    }
    // constant parameters: f, mu
    else {
      double par_val[nPtBins], par_err[nPtBins];
      for(int i = 0; i < nPtBins; i++) {
	par_val[i] = fitS->GetParameter(p_pos[ct_pos]);
	par_err[i] = fitS->GetParError(p_pos[ct_pos]);
      }
      TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, par_val, pt_err, par_err);
      g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
      ct_pos++;
    }
  }
  TGraphErrors *g_fNP = new TGraphErrors(nPtBins, pt_val, fracNP, pt_err, efz);
  g_fNP->Write("fit_fNP");
  
  TLine *l_chi = new TLine(ptBins[0], fitS->GetChisquare()/fitS->GetNDF(), ptBins[nPtBins], fitS->GetChisquare()/fitS->GetNDF());
  l_chi->Write("fit_chiN");
  fout->Close();

  ofstream ftex;
  ftex.open(Form("text_output/tfit_res.tex"));
  ftex << "\\begin{tabular}{c||c|c|c|c|c||c}\n";
  ftex << "$\\pt$ (GeV) & $N_{PR}$ & $N_{NP}$ & $\\sigma_1$ ($\\mu$m) & $\\sigma_2$ ($\\mu$m)  & $\\lambda$ ($\\mu$m) & $f_{NP}$ (\\%) \\\\\n";
  ftex << "\\hline\n";

  for(int i = 0; i < nPtBins; i++) {
    // pT bin
    ftex << Form("$[%.0f, %.0f]$", ptBins[i], ptBins[i+1]);
    for(int i_p = 0; i_p < 5; i_p++) {
      double mult = 1.;
      if(i_p) < 2) mult = 1./(ptBins[i+1]-ptBins[i]);
      else if(i_p >= 2) mult = 1e3;
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
    ftex << " & " << setprecision(2) << fixed << fracNP[i]*100.;
    ftex <<  "\\\\\n";
  }
  ftex << "\\end{tabular}\n";
  ftex.close();

  ofstream fout2;
  fout2.open(Form("text_output/tfit_resA.tex"));
  fout2 << "\\begin{tabular}{c|c||c}\n";
  fout2 << " $f$ $(\\%)$ & $\\mu$ $($\\mu$m) & $\\chi^2/$ndf \\\\\n";
  fout2 << "\\hline\n";

  // f
  double val = fitS->GetParameter(p_pos[0])*100.;
  double unc = fitS->GetParError(p_pos[0])*100.;
  int p_norm = 1.;
  if(unc < 1) p_norm = ceil(-log10(unc))+1;	
  fout2 <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
  // mu
  val = fitS->GetParameter(p_pos[1])*1000.;
  unc = fitS->GetParError(p_pos[1])*1000.;
  p_norm = 1.;
  if(unc < 1) p_norm = ceil(-log10(unc))+1;	
  fout2 << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
  // chi^2
  fout2 << setprecision(0) << fitS->GetChisquare() << "/" << fitS->GetNDF() << "\\\\\n";
  fout2 << "\\end{tabular}\n";
  fout2.close();


  cout << fitS->GetChisquare() << "/" << fitS->GetNDF() << endl;
  
  c->Destructor();

}
