//pt bins defined globally for access from functions
#import "../../ptbins.C"

// functions to access within other functions
TF1 *fres;
TF1 *fNP_SR, *fNP_bkg;
TF1 *fexp_bkg;

// define negative exponential only for positive x
double pos_exp_p(double x, double ld)
{
  if(x > 0) return exp(-x/ld);
  else return 0;
}

// define final fit function summing the PR and NP contributions
// parameters: N_PR, N_NP, f1, f2, mu, sigma1, sig2/sig1, sig3/sig1, tNP
double func_sum(double *xx, double *pp)
{
  // get pt bin
  int pt_bin;
  for(int i = 0; i < nPtBins; i++)
    if(ptBins[i] < xx[1] && ptBins[i+1] > xx[1])
      pt_bin = i;

  // direct pars
  double N_PR = pp[pt_bin];
  double N_NP = pp[nPtBins+pt_bin];
  double f1 = pp[2*nPtBins]; // constant in pT
  double f2 = pp[3*nPtBins]; // constant in pT
  double mu = pp[4*nPtBins]; // constant in pT
  double sig1 = pp[5*nPtBins+pt_bin];
  double sigR21 = pp[6*nPtBins]; // constant in pT
  double sigR31 = pp[7*nPtBins]; // constant in pT
  double tnp_SR = pp[8*nPtBins+pt_bin];
  double N_bkg = pp[9*nPtBins+pt_bin]; // fixed in each bin
  double tnp = pp[10*nPtBins+pt_bin]; // fixed in each bin

  // indirect pars
  double sig2 = sigR21*sig1;
  double sig3 = sigR31*sig1;

  fres->SetParameters(N_PR, f1, f2, mu, sig1, sig2, sig3);
  fNP_SR->SetParameters(N_NP, f1, f2, mu, sig1, sig2, sig3, tnp_SR);
  fNP_bkg->SetParameters(1, f1, f2, mu, sig1, sig2, sig3, tnp);

  // scale N_bkg par
  fexp_bkg->SetParameters(N_bkg, tnp);
  double N_bkg_conv = fexp_bkg->Eval(0.4)/fNP_bkg->Eval(0.4);
  fNP_bkg->SetParameter(0, N_bkg_conv);
    
  return fres->Eval(xx[0]) + fNP_SR->Eval(xx[0]) + fNP_bkg->Eval(xx[0]);
}

double nbkg_eval(double N_PR, double N_NP, double f1, double f2, double mu, double sig1, double sig2, double sig3, double tnp_SR, double N_bkg, double tnp)
{
  fNP_bkg->SetParameters(1, f1, f2, mu, sig1, sig2, sig3, tnp);

  // scale N_bkg par
  fexp_bkg->SetParameters(N_bkg, tnp);
  double N_bkg_conv = fexp_bkg->Eval(0.4)/fNP_bkg->Eval(0.4);

  return N_bkg_conv;
}

// define final fit function summing the PR and NP contributions
double sum_1d(double *xx, double *pp)
{
  double lt = xx[0];
  double N_PR = pp[0], N_NP = pp[1], f1 = pp[2], f2 = pp[3], mu = pp[4], sig1 = pp[5], sig2 = pp[6], sig3 = pp[7], tnp_SR = pp[8], N_bkg = pp[9], tnp = pp[10];

  fres->SetParameters(N_PR, f1, f2, mu, sig1, sig2, sig3);
  fNP_SR->SetParameters(N_NP, f1, f2, mu, sig1, sig2, sig3, tnp_SR);
  fNP_bkg->SetParameters(N_bkg, f1, f2, mu, sig1, sig2, sig3, tnp);

  return fres->Eval(lt) + fNP_SR->Eval(lt) + fNP_bkg->Eval(lt);
}

// define 1d function summing only the NP contributions
double sum_NP(double *xx, double *pp)
{
  double lt = xx[0];
  double N_PR = pp[0], N_NP = pp[1], f1 = pp[2], f2 = pp[3], mu = pp[4], sig1 = pp[5], sig2 = pp[6], sig3 = pp[7], tnp_SR = pp[8], N_bkg = pp[9], tnp = pp[10];

  fNP_SR->SetParameters(N_NP, f1, f2, mu, sig1, sig2, sig3, tnp_SR);
  fNP_bkg->SetParameters(N_bkg, f1, f2, mu, sig1, sig2, sig3, tnp);

  return fNP_SR->Eval(lt) + fNP_bkg->Eval(lt);
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


void ltBkg2d()
{
  // prepare binning and histograms for plots
  TH2D *h_d2d = new TH2D();  
  TFile *fin = new TFile("../../PR_fit/files/ltStore.root");
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
  fres = new TF1("fres", "[0]*([1]*TMath::Gaus(x, [3],[4]) + [2]*TMath::Gaus(x, [3], [5])+ (1.-[1]-[2])*TMath::Gaus(x, [3], [6]))", 5*lowt, 5*hit);
  
  // define the NP function by convolution
  TF1 *fexp = new TF1("fexp", "pos_exp_p(x,[0])", 5*lowt, 5*hit);
  TF1Convolution *fcNP = new TF1Convolution(fres, fexp, 5*lowt, 5*hit);
  fcNP->SetRange(5*lowt, 5*hit);
  fcNP->SetNofPointsFFT(1000);
 fNP_SR = new TF1("fNP_SR", *fcNP, lowt, hit, fcNP->GetNpar());

  // define the NP bkg function by convolution
  TF1 *fexp2 = new TF1("fexp2", "pos_exp_p(x,[0])", 5*lowt, 5*hit);
  TF1Convolution *fcNP_bkg = new TF1Convolution(fres, fexp2, 5*lowt, 5*hit);
  fcNP_bkg->SetRange(5*lowt, 5*hit);
  fcNP_bkg->SetNofPointsFFT(1000);
  fNP_bkg = new TF1("fNP_bkg", *fcNP_bkg, lowt, hit, fcNP_bkg->GetNpar());

  // define the bkg exp for N_bkg adjustment
  fexp_bkg = new TF1("fexp_bkg", "[0]*pos_exp_p(x,[1])", 5*lowt, 5*hit);
  
  // get the fixed values of the SB part
  TFile *fin_SB = new TFile("../../SBLtFits/files/store_SB.root");
  string parsav[] = {"N_NP", "t_NP"};
  double par_bkg[2][nPtBins], epar_bkg[2][nPtBins];
  for(int i = 0; i < 2; i++){
    TGraphErrors *g_par = (TGraphErrors*)fin_SB->Get(Form("g_%s", parsav[i].c_str()));
    for(int j = 0; j < nPtBins; j++) {
      par_bkg[i][j] = g_par->GetY()[j];
      epar_bkg[i][j] = g_par->GetEY()[j];
    }
  }
  fin_SB->Close();
  
  // SYST CHECK - run N-unc(N)
  for(int j = 0; j < nPtBins; j++) {
    par_bkg[0][j] -= epar_bkg[0][j];
  }
  
  TF2 *fitS = new TF2("fitS", func_sum, lowt, hit, ptBins[0], ptBins[nPtBins], 11*nPtBins, 2);
  string par_n[] = {"N_PR", "N_NP", "f1", "f2", "mu", "sigma1", "sig2/sig1", "sig3/sig1", "t_NPSR", "N_bkg", "t_NP"};
  double par_v[] =  {1, 1, 0.55, 0.4, 0, 1e-2, 1.5, 3.5, 0.3, 1, 1};

  // define the parameters
  for(int i = 0; i < nPtBins; i++) {
    // normalizations
    fitS->SetParName(i, Form("%s_%d", par_n[0].c_str(), i));
    fitS->SetParameter(i, h_d1d[i]->GetMaximum());
    fitS->SetParName(i+nPtBins, Form("%s_%d", par_n[1].c_str(), i));
    fitS->SetParameter(i+nPtBins, h_d1d[i]->GetMaximum());

    for(int i_p = 2; i_p < 11; i_p++) {
      fitS->SetParName(i+i_p*nPtBins, Form("%s_%d", par_n[i_p].c_str(), i));
      fitS->SetParameter(i+i_p*nPtBins, par_v[i_p]);
      // setting the constant parameters f1, f2, mu, sig2/sig1, sig3/sig1
      if((i_p < 5 || i_p == 6 || i_p == 7) && i > 0) fitS->FixParameter(i+i_p*nPtBins, par_v[i_p]);
      // fixing the parameters from the bkg part
      if(i_p > 8) fitS->FixParameter(i+i_p*nPtBins, par_bkg[i_p-9][i]);
    }
  }

  // fit the 2d function to the lifetime:pT map
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.13);

  TFitResultPtr fitres = h_d2d->Fit("fitS", "SVR");

  // tf1 for plotting in the 1D bins
  // separate parts of the fit function - given by fres and fNP
  TF1 *f_1d = new TF1("f_1d", sum_1d, lowt, hit, 11);

  double pt_val[nPtBins], pt_err[nPtBins];
  double pars[11][nPtBins], epars[11][nPtBins];
  double rpars[11][nPtBins];

  // cycle over all pT bins
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    pt_val[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
    pt_err[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);

    // storing parameters
    for(int j = 0; j < 11; j++) {
      if(j == 2 || j == 3 || j == 4 || j == 6 || j == 7) { // constant pars
	pars[j][i_pt] = fitS->GetParameter(j*nPtBins);
	epars[j][i_pt] = fitS->GetParError(j*nPtBins);
      }
      else {
	pars[j][i_pt] = fitS->GetParameter(j*nPtBins + i_pt);
	epars[j][i_pt] = fitS->GetParError(j*nPtBins + i_pt);
      }
      // sigma ratios turned into sigmas
      if(j == 6 || j == 7) rpars[j][i_pt] = pars[5][i_pt]*pars[j][i_pt];
      else rpars[j][i_pt] = pars[j][i_pt];
    }

    // getting corrected N_bkg
    double N_bkg_conv = nbkg_eval(rpars[0][i_pt], rpars[1][i_pt], rpars[2][i_pt], rpars[3][i_pt], rpars[4][i_pt], rpars[5][i_pt], rpars[6][i_pt], rpars[7][i_pt], rpars[8][i_pt], rpars[9][i_pt], rpars[10][i_pt]);

    // initializing f_1d and plotting
    for(int i = 0; i < 11; i++)
      f_1d->SetParameter(i, rpars[i][i_pt]);
    f_1d->SetParameter(9, N_bkg_conv);

    c->SetTopMargin(0.015);
    c->SetLogy();
     
    h_d1d[i_pt]->SetMaximum(h_d1d[i_pt]->GetMaximum()*1.2);
    h_d1d[i_pt]->SetMinimum(h_d1d[i_pt]->GetMaximum()*5e-3);

    TH1F *fh = c->DrawFrame(lowPlot, h_d1d[i_pt]->GetMinimum(), hit, h_d1d[i_pt]->GetMaximum());
    fh->SetXTitle("c#tau (mm)");
    fh->SetYTitle(Form("Events per %.0f #mum", wbin*1000.));
    fh->GetYaxis()->SetTitleOffset(1.8);
    fh->GetYaxis()->SetLabelOffset(0.01);
    fh->GetXaxis()->SetLabelOffset(0.015);
    fh->GetXaxis()->SetTitleOffset(1.3);

    h_d1d[i_pt]->SetMarkerStyle(20);
    h_d1d[i_pt]->SetMarkerColor(kBlack);
    h_d1d[i_pt]->SetLineColor(kBlack);
    h_d1d[i_pt]->SetMarkerSize(0.75);
    h_d1d[i_pt]->Draw("error same");
    
    f_1d->SetLineColor(kBlue);
    f_1d->Draw("lsame");

    // draw fit contributions
    fres->SetParameters(rpars[0][i_pt], rpars[2][i_pt], rpars[3][i_pt], rpars[4][i_pt], rpars[5][i_pt], rpars[6][i_pt], rpars[7][i_pt]);
    fres->SetLineStyle(kDashed);
    fres->SetLineColor(kGreen+3);
    fres->Draw("lsame");
    // NP = sum of NP psi + NP bkg terms
    TF1 *fNP = new TF1("fnp", sum_NP, lowt, hit, 13);
    for(int i = 0; i < 11; i++)
      fNP->SetParameter(i, rpars[i][i_pt]);
    fNP->SetParameter(9, N_bkg_conv);
    fNP->SetLineStyle(kDashed);
    fNP->SetLineColor(kViolet);
    fNP->Draw("lsame");

   // plot each NP term
    // NP psi
    fNP_SR->SetParameters(rpars[1][i_pt], rpars[2][i_pt], rpars[3][i_pt], rpars[4][i_pt], rpars[5][i_pt], rpars[6][i_pt], rpars[7][i_pt], rpars[8][i_pt]);
    fNP_SR->SetLineStyle(kDashed);
    fNP_SR->SetLineColor(kRed+1);
    fNP_SR->Draw("lsame");
    // NP bkg
    fNP_bkg->SetParameters(N_bkg_conv, rpars[2][i_pt], rpars[3][i_pt], rpars[4][i_pt], rpars[5][i_pt], rpars[6][i_pt], rpars[7][i_pt], rpars[10][i_pt]);
    fNP_bkg->SetLineStyle(kDashed);
    fNP_bkg->SetLineColor(kYellow+3);
    fNP_bkg->Draw("lsame");

    // individual Gaussians
    TF1 **fGs = new TF1*[3];
    fGs[0] = new TF1("fGs_0", "[0]*[1]*TMath::Gaus(x, [2],[3])", 5*lowt, 5*hit);
    fGs[0]->SetParameters(f_1d->GetParameter(0), f_1d->GetParameter(2), f_1d->GetParameter(4), f_1d->GetParameter(5));
    fGs[0]->SetLineStyle(kDashDotted);
    fGs[0]->SetLineColor(kRed+2);
    fGs[0]->Draw("lsame");

    fGs[1] = new TF1("fGs_1", "[0]*[1]*TMath::Gaus(x, [2],[3])", 5*lowt, 5*hit);
    fGs[1]->SetParameters(f_1d->GetParameter(0), f_1d->GetParameter(3), f_1d->GetParameter(4), f_1d->GetParameter(6));
    fGs[1]->SetLineStyle(kDashDotted);
    fGs[1]->SetLineColor(kOrange+2);
    fGs[1]->Draw("lsame");

    fGs[2] = new TF1("fGs_2", "[0]*[1]*TMath::Gaus(x, [2],[3])", 5*lowt, 5*hit);
    fGs[2]->SetParameters(f_1d->GetParameter(0), 1.-f_1d->GetParameter(2)-f_1d->GetParameter(3), f_1d->GetParameter(4), f_1d->GetParameter(7));
    fGs[2]->SetLineStyle(kDashDotted);
    fGs[2]->SetLineColor(kYellow+3);
    fGs[2]->Draw("lsame");

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

    TLatex lc;
    lc.SetTextSize(0.03);

    // draw CMS text
    double xp = getPos(lowPlot, hit, 0.35, 0);
    double yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.95, 1);
    lc.DrawLatex(xp, yp, "CMS");
    // draw L
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.9, 1);
    lc.DrawLatex(xp, yp, "#bf{L = 103.3 fb^{-1}}");
    // draw sqrt(s)
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.85, 1);
    lc.DrawLatex(xp, yp, "#bf{#sqrt{s} = 13 TeV}");
    // draw pT
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.75, 1);
    lc.DrawLatex(xp, yp, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV}", ptBins[i_pt], ptBins[i_pt+1]));
    // draw y
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.7, 1);
    lc.DrawLatex(xp, yp, "#bf{|#it{y}| < 1.2}");
    // draw the chi^2/ndf
    xp = getPos(lowPlot, hit, 0.7, 0);
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.9, 1);
    lc.DrawLatex(xp, yp, Form("#bf{#chi^{2}/ndf = %.0f / %d}", fitS->GetChisquare(), fitS->GetNDF()));
    // draw the state
    lc.SetTextSize(0.04);
    xp = getPos(lowPlot, hit, 0.6, 0);
    yp = getPos(h_d1d[i_pt]->GetMinimum(), h_d1d[i_pt]->GetMaximum(), 0.5, 1);
    lc.DrawLatex(xp, yp, "#bf{J/#psi}");

    TLegend *leg = new TLegend(0.424, 0.1, 0.724, 0.55);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(kWhite,0);
    leg->AddEntry(h_d1d[i_pt], "Data", "pl");
    leg->AddEntry(f_1d, "Total fit", "l");
    leg->AddEntry(fres, "Prompt cont.", "l");
    for(int i = 0; i < 3; i++){
      leg->AddEntry(fGs[i], Form("G_{%d}", i+1), "l");
    }
    leg->AddEntry(fNP, "Non-prompt cont.", "l");
    leg->AddEntry(fNP_SR, "Non-prompt #psi(2S) cont.", "l");
    leg->AddEntry(fNP_bkg, "Non-prompt bkg cont.", "l");
    leg->Draw();

    c->SaveAs(Form("plots/lifetime2d/fit/fit_%d.pdf", i_pt));
    c->Clear();

    // calculating pulls
    double tv[tbins], pv[tbins], dv[tbins];
    for(int i_t = 0 ; i_t < tbins; i_t++) {
      tv[i_t] = h_d1d[i_pt]->GetBinCenter(i_t+1);
      double fitv = f_1d->Eval(tv[i_t]);
      double datav = h_d1d[i_pt]->GetBinContent(i_t+1);
      double datau = h_d1d[i_pt]->GetBinError(i_t+1);
      if(datau > 0 ) pv[i_t] = (datav-fitv)/datau;
      else pv[i_t] = 0;
      if(fitv > 0 ) dv[i_t] = (datav-fitv)/fitv * 100;
      else dv[i_t] = 0;
    }
    
    // plotting the pulls
    c->SetTopMargin(0.1);
    c->SetLogy(0);
    
    TH1F *fp = c->DrawFrame(lowPlot, -9, hit, 9);
    fp->SetXTitle("c#tau (mm)");
    fp->SetYTitle("pulls");
    fp->GetYaxis()->SetTitleOffset(1.3);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("Lifetime fit pulls (%.1f < p_{T} < %.1f GeV)", ptBins[i_pt], ptBins[i_pt+1]));
  
    TGraph *g_pull = new TGraph(tbins, tv, pv);
    g_pull->SetLineColor(kBlack);
    g_pull->SetMarkerColor(kBlack);
    g_pull->SetMarkerStyle(20);
    g_pull->Draw("psame");

    // aux lines - pull=0 and sigma limits
    TF1 *fcons = new TF1("fcons", "[0]", lowPlot, hit);
    fcons->SetParameter(0, 0);
    fcons->SetLineStyle(kDashed);
    fcons->SetLineColor(kBlack);
    fcons->Draw("lsame");

    TLine *psig1 = new TLine(-pr_lim, -9, -pr_lim, 9);
    psig1->SetLineStyle(kDashed);
    psig1->Draw("lsame");
    TLine *psig2 = new TLine(pr_lim, -9, pr_lim, 9);
    psig2->SetLineStyle(kDashed);
    psig2->Draw("lsame");
    TLine *psig3 = new TLine(np_lim, -9, np_lim, 9);
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

    c->SaveAs(Form("plots/lifetime2d/fit/pulls_%d.pdf", i_pt));
    c->Clear();

    // plotting the devs
    TH1F *fd = c->DrawFrame(lowPlot, -15, hit, 15);
    fd->SetXTitle("c#tau (mm)");
    fd->SetYTitle("deviation");
    fd->GetYaxis()->SetTitleOffset(1.3);
    fd->GetYaxis()->SetLabelOffset(0.01);
    fd->SetTitle(Form("Lifetime fit deviations (%.1f < p_{T} < %.1f GeV)", ptBins[i_pt], ptBins[i_pt+1]));

    TGraph *g_dev = new TGraph(tbins, tv, dv);
    g_dev->SetLineColor(kBlack);
    g_dev->SetMarkerColor(kBlack);
    g_dev->SetMarkerStyle(20);
    g_dev->Draw("psame");
  
    // aux lines - pull = 0 and sigma limits
    fcons->Draw("lsame");

    TLine *dsig1 = new TLine(-pr_lim, -15, -pr_lim, 15);
    dsig1->SetLineStyle(kDashed);
    dsig1->Draw("lsame");
    TLine *dsig2 = new TLine(pr_lim, -15, pr_lim, 15);
    dsig2->SetLineStyle(kDashed);
    dsig2->Draw("lsame");
    TLine *dsig3 = new TLine(np_lim, -15, np_lim, 15);
    dsig3->SetLineStyle(kDashed);
    dsig3->Draw("lsame");

    c->SaveAs(Form("plots/lifetime2d/fit/devs_pt%d.pdf", i_pt));
    c->Clear();
  }
  
  // storing the free parameters
  TFile *fout = new TFile("files/ltfitres2d.root", "recreate");
  string parlab[] = {"N_PR", "N_NP", "f1", "f2", "mu", "sig1", "sigR21", "sigR31", "t_NPSR", "N_bkg", "t_NP"};

  for(int i_p = 0; i_p < 11; i_p++) {
    TGraphErrors *g_par = new TGraphErrors(nPtBins, pt_val, pars[i_p], pt_err, epars[i_p]);
    g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
  }
  
  fitres->SetName("fitres");
  fitres->Write();
  fout->Close();

  cout << fitS->GetChisquare() << "/" << fitS->GetNDF() << endl;
  
  c->Destructor();

}
