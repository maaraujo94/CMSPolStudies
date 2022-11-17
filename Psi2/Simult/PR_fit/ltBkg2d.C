//pt bins defined globally for access from functions
#import "../ptbins.C"

// functions to access within other functions
TF1 *fres;
TF1 *fNP;

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
  fres->SetParameters(N_PR, f, mu, sig1, sig2);
  fNP->SetParameters(N_NP, f, mu, sig1, sig2, ld);

  return fres->Eval(xx[0]) + fNP->Eval(xx[0]);
}

// define final fit function summing the PR and NP contributions
double sum_1d(double *xx, double *pp)
{
  double lt = xx[0];
  double N_PR = pp[0], N_NP = pp[1], f2 = pp[2], mu = pp[3], sig1 = pp[4], sig2 = pp[5], ld = pp[6];

  double inp[] = {lt};
  fres->SetParameters(N_PR, f2, mu, sig1, sig2);
  fNP->SetParameters(N_NP, f2, mu, sig1, sig2, ld);

  return fres->Eval(lt) + fNP->Eval(lt);
}

void ltBkg2d()
{
  // prepare binning and histograms for plots
  TH2D *h_d2d = new TH2D();  
  TFile *fin = new TFile("files/ltStore.root");
  fin->GetObject("ltH", h_d2d);
  h_d2d->SetDirectory(0);
  fin->Close();

  int tbins = h_d2d->GetNbinsX();
  double lowt = h_d2d->GetXaxis()->GetBinLowEdge(1);
  double hit = h_d2d->GetXaxis()->GetBinUpEdge(tbins);
  double wbin = (hit-lowt)/(double)tbins;

  // Make 1d histos
  TH1D **h_d1d = new TH1D*[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    h_d1d[i] = h_d2d->ProjectionX(Form("ltH%.0f", ptBins[i]), i+1, i+1);
    h_d1d[i]->SetTitle(Form("Run 2 data c#tau (%.1f < p_{T} < %.1f GeV)", ptBins[i], ptBins[i+1]));
  }

  // define aux vals for plotting
  double pr_lim = 0.05;
  double np_lim = 0.1;
  double lowPlot = -0.1;

  // define the resolution (=PR) function
  fres = new TF1("fres", "[0]*([1]*TMath::Gaus(x, [2],[3]) + (1.-[1])*TMath::Gaus(x, [2], [4]))", 5*lowt, 5*hit);
  
  // define the NP function by convolution
  TF1 *fexp = new TF1("fexp", "pos_exp(x,[0])", 5*lowt, 5*hit);
  TF1Convolution *fcNP = new TF1Convolution(fres, fexp, 5*lowt, 5*hit);
  fcNP->SetRange(5*lowt, 5*hit);
  fcNP->SetNofPointsFFT(1000);
  fNP = new TF1("fNP", *fcNP, lowt, hit, fcNP->GetNpar());

  TF2 *fitS = new TF2("fitS", func_sum, lowt, hit, ptBins[0], ptBins[nPtBins], 2+5*nPtBins, 2);

  // define constant parameters - f, mu
  fitS->SetParName(2*nPtBins, "f");
  fitS->SetParameter(2*nPtBins, 0.75);
  fitS->SetParName(2*nPtBins+1, "mu");
  fitS->SetParameter(2*nPtBins+1, 0);
  
  // define free parameters - N, sigma, lambda
  for(int i = 0; i < nPtBins; i++) {
    fitS->SetParName(i, Form("N_PR_%d", i));
    fitS->SetParameter(i, h_d1d[i]->GetMaximum());
    
    fitS->SetParName(i+nPtBins, Form("N_NP_%d", i));
    fitS->SetParameter(i+nPtBins, h_d1d[i]->GetMaximum()*5.);
    
    fitS->SetParName(i+2*nPtBins+2, Form("sig1_%d", i));
    fitS->SetParameter(i+2*nPtBins+2, 1e-2);

    fitS->SetParName(i+3*nPtBins+2, Form("sig2_%d", i));
    fitS->SetParameter(i+3*nPtBins+2, 2.5e-2);

    fitS->SetParName(i+4*nPtBins+2, Form("lambda_%d", i));
    fitS->SetParameter(i+4*nPtBins+2, 0.35);
  }
  // fit the 2d function to the lifetime:pT map
  TCanvas *c = new TCanvas("", "", 700, 700);
  cout << "starting fit" << endl;
  TFitResultPtr fitres = h_d2d->Fit("fitS", "SVR");

  cout << "output right after fit" << endl;
  
  // tf1 for plotting in the 1D bins
  // separate parts of the fit function - given by fres and fNP
  TF1 *f_1d = new TF1("f_1d", sum_1d, lowt, hit, 7);
  f_1d->SetParNames("N_PR", "N_NP", "f", "mu", "sigma1", "sigma2", "lambda");

  double pt_val[nPtBins], pt_err[nPtBins];
  double pars[7][nPtBins], epars[7][nPtBins];
  double fracNP[nPtBins], efz[nPtBins];

  // cycle over all pT bins
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    cout << i_pt << endl;
    
    pt_val[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
    pt_err[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);

    // storing all parameters - N_PR, N_NP, f, mu, sig_1, sig_2, ld
    pars[0][i_pt] = fitS->GetParameter(i_pt);
    pars[1][i_pt] = fitS->GetParameter(i_pt+nPtBins);
    pars[2][i_pt] = fitS->GetParameter(2*nPtBins);
    pars[3][i_pt] = fitS->GetParameter(2*nPtBins+1);
    pars[4][i_pt] = fitS->GetParameter(i_pt+2+2*nPtBins);
    pars[5][i_pt] = fitS->GetParameter(i_pt+2+3*nPtBins);
    pars[6][i_pt] = fitS->GetParameter(i_pt+2+4*nPtBins);
    epars[0][i_pt] = fitS->GetParError(i_pt);
    epars[1][i_pt] = fitS->GetParError(i_pt+nPtBins);
    epars[2][i_pt] = fitS->GetParError(2*nPtBins);
    epars[3][i_pt] = fitS->GetParError(2*nPtBins+1);
    epars[4][i_pt] = fitS->GetParError(i_pt+2+2*nPtBins);
    epars[5][i_pt] = fitS->GetParError(i_pt+2+3*nPtBins);
    epars[6][i_pt] = fitS->GetParError(i_pt+2+4*nPtBins);

    // initializing f_1d and plotting
    f_1d->SetParameters(pars[0][i_pt],
			pars[1][i_pt],
			pars[2][i_pt],
			pars[3][i_pt],
			pars[4][i_pt],
			pars[5][i_pt],
			pars[6][i_pt]);


    c->SetLogy();
     
    h_d1d[i_pt]->SetMaximum(h_d1d[i_pt]->GetMaximum()*1.2);
    h_d1d[i_pt]->SetMinimum(h_d1d[i_pt]->GetMaximum()*5e-4);

    TH1F *fh = c->DrawFrame(lowPlot, h_d1d[i_pt]->GetMinimum(), hit, h_d1d[i_pt]->GetMaximum());
    fh->SetXTitle("c#tau (mm)");
    fh->SetYTitle(Form("Events per %.0f #mum", wbin*1000.));
    fh->GetYaxis()->SetTitleOffset(1.3);
    fh->GetYaxis()->SetLabelOffset(0.01);
    fh->SetTitle(h_d1d[i_pt]->GetTitle());

    h_d1d[i_pt]->SetMarkerStyle(20);
    h_d1d[i_pt]->SetMarkerColor(kBlack);
    h_d1d[i_pt]->SetLineColor(kBlack);
    h_d1d[i_pt]->SetMarkerSize(0.75);
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

    c->SaveAs(Form("plots/lifetime2d/fit_%d.pdf", i_pt));
    c->Clear();

    // get the NP fraction in the signal region (+- 100 mum)
    double evt_NP = fNP->Integral(-pr_lim, pr_lim);
    
    double min_bin = h_d1d[i_pt]->GetXaxis()->FindBin(-pr_lim+1e-6);
    double max_bin = h_d1d[i_pt]->GetXaxis()->FindBin(pr_lim-1e-6);
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
      if(fitv > 0 ) dv[i_t] = (datav-fitv)/fitv ;
      else dv[i_t] = 0;
    }
    
    // plotting the pulls
    c->SetLogy(0);
    
    TH1F *fp = c->DrawFrame(lowPlot, -7, hit, 7);
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

    TLine *psig1 = new TLine(-pr_lim, -7, -pr_lim, 7);
    psig1->SetLineStyle(kDashed);
    psig1->Draw("lsame");
    TLine *psig2 = new TLine(pr_lim, -7, pr_lim, 7);
    psig2->SetLineStyle(kDashed);
    psig2->Draw("lsame");
    TLine *psig3 = new TLine(np_lim, -7, np_lim, 7);
    psig3->SetLineStyle(kDashed);
    psig3->Draw("lsame");

    TLine *plim1 = new TLine(lowPlot, -5, hit, -5);
    plim1->SetLineStyle(kDotted);
    plim1->Draw("lsame");
    TLine *plim2 = new TLine(lowPlot, -3, hit, -3);
    plim2->SetLineStyle(kDotted);
    plim2->Draw("lsame");
    TLine *plim3 = new TLine(lowPlot, 3, hit, 3);
    plim3->SetLineStyle(kDotted);
    plim3->Draw("lsame");
    TLine *plim4 = new TLine(lowPlot, 5, hit, 5);
    plim4->SetLineStyle(kDotted);
    plim4->Draw("lsame");

    c->SaveAs(Form("plots/lifetime2d/pulls_%d.pdf", i_pt));
    c->Clear();

    // plotting the devs
    TH1F *fd = c->DrawFrame(lowPlot, -1, hit, 1);
    fd->SetXTitle("c#tau (mm)");
    fd->SetYTitle("deviation");
    fd->GetYaxis()->SetTitleOffset(1.3);
    fd->GetYaxis()->SetLabelOffset(0.01);
    fd->SetTitle(Form("Lifetime fit deviations (%.0f < p_{T} < %.0f GeV)", ptBins[i_pt], ptBins[i_pt+1]));

    TGraph *g_dev = new TGraph(tbins, tv, dv);
    g_dev->SetLineColor(kBlack);
    g_dev->SetMarkerColor(kBlack);
    g_dev->SetMarkerStyle(20);
    g_dev->SetMarkerSize(.75);
    g_dev->Draw("psame");
  
    // aux lines - pull = 0 and sigma limits
    fcons->Draw("lsame");

    TLine *dsig1 = new TLine(-pr_lim, -1, -pr_lim, 1);
    dsig1->SetLineStyle(kDashed);
    dsig1->Draw("lsame");
    TLine *dsig2 = new TLine(pr_lim, -1, pr_lim, 1);
    dsig2->SetLineStyle(kDashed);
    dsig2->Draw("lsame");
    TLine *dsig3 = new TLine(np_lim, -1, np_lim, 1);
    dsig3->SetLineStyle(kDashed);
    dsig3->Draw("lsame");

    c->SaveAs(Form("plots/lifetime2d/devs_pt%d.pdf", i_pt));
    c->Clear();
  }
  
  // storing the free parameters
  TFile *fout = new TFile("files/ltfitres2d.root", "recreate");
  string parlab[] = {"N_PR", "N_NP", "f", "mu", "sig1", "sig2", "lambda"};

  for(int i_p = 0; i_p < 7; i_p++) {
    TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[i_p], pt_err, epars[i_p]);
    g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
  }
  TGraphErrors *g_fNP = new TGraphErrors(nPtBins, pt_val, fracNP, pt_err, efz);
  g_fNP->Write("fit_fNP");
  
  fitres->SetName("fitres");
  fitres->Write();
  fout->Close();

  double mult[] = {1., 1., 100., 1e3, 1e3, 1e3, 1e3};

  ofstream ftex;
  ftex.open(Form("text_output/ltfit_res.tex"));
  ftex << "\\begin{tabular}{c||c|c|c|c|c|c|c||c}\n";
  ftex << "$\\pt$ (GeV) & $N_{PR}$ & $N_{NP}$ & f (\\%) & $\\mu$ ($\\mu$m) & $\\sigma_1$ ($\\mu$m) & $\\sigma_2$ ($\\mu$m)  & t$_{NP}$ ($\\mu$m) & $f_{NP}$ (\\%) \\\\\n";
  ftex << "\\hline\n";

  for(int i = 0; i < nPtBins; i++) {
    // pT bin
    ftex << Form("$[%.0f, %.0f]$", ptBins[i], ptBins[i+1]);
    for(int i_p = 0; i_p < 7; i_p++) {
      // plot fixed values - f (2), mu (3)
      if((i_p == 2 || i_p == 3) && i == 0) {
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
      else if((i_p == 2 || i_p == 3) && i != 0)
	ftex << " & ";
      else {
	double val = pars[i_p][i]*mult[i_p], unc = epars[i_p][i]*mult[i_p];
	if(i_p < 2) {
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
    }
    ftex << " & " << setprecision(2) << fixed << fracNP[i]*100.;
    ftex <<  "\\\\\n";
  }
  ftex << "\\end{tabular}\n";
  ftex.close();

  cout << fitS->GetChisquare() << "/" << fitS->GetNDF() << endl;
  
  c->Destructor();

}
