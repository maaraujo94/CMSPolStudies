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
double sum_1d(double *xx, double *pp)
{
  double lt = xx[0];
  double N_PR = pp[0], N_NP = pp[1], f2 = pp[2], mu = pp[3], sig1 = pp[4], sig2 = pp[5], ld = pp[6];

  double inp[] = {lt};
  fres->SetParameters(N_PR, f2, mu, sig1, sig2);
  fNP->SetParameters(N_NP, f2, mu, sig1, sig2, ld);

  return fres->Eval(lt) + fNP->Eval(lt);
}

// macro to plot the lifetime dist + fit 
void plotLtClean()
{
  // get the lifetime distributions
  // prepare binning and histograms for plots
  TH2D *h_d2d = new TH2D();  
  TFile *fin = new TFile("../PR_fit/files/ltStore.root");
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
    h_d1d[i]->SetTitle(Form("Run 2 data c#tau (%.0f < p_{T} < %.0f GeV)", ptBins[i], ptBins[i+1]));
  }

  // define aux vals for plotting
  double pr_lim = 0.05;
  double np_lim = 0.1;
  double lowPlot = -0.1;

  // get the fit parameters
  TFile *fin2 = new TFile("../PR_fit/files/ltfitres2d.root");
  string parlab[] = {"N_PR", "N_NP", "f", "mu", "sig1", "sig2", "lambda"};
  double pars[7][nPtBins];
  for(int i_p = 0; i_p < 7; i_p++) {
    TGraphErrors *g_par = (TGraphErrors*)fin2->Get(Form("fit_%s", parlab[i_p].c_str()));
    for(int i = 0; i < nPtBins; i++) 
      pars[i_p][i] = g_par->GetY()[i];
  }
  fin2->Close();
  
  // define the resolution (=PR) function
  fres = new TF1("fres", "[0]*([1]*TMath::Gaus(x, [2],[3]) + (1.-[1])*TMath::Gaus(x, [2], [4]))", 5*lowt, 5*hit);
  
  // define the NP function by convolution
  TF1 *fexp = new TF1("fexp", "pos_exp(x,[0])", 5*lowt, 5*hit);
  TF1Convolution *fcNP = new TF1Convolution(fres, fexp, 5*lowt, 5*hit);
  fcNP->SetRange(5*lowt, 5*hit);
  fcNP->SetNofPointsFFT(1000);
  fNP = new TF1("fNP", *fcNP, lowt, hit, fcNP->GetNpar());
  
  // tf1 and th1 for plotting in the 1D bins
  // separate parts of the fit function - given by fres and fNP
  TF1 *f_1d = new TF1("f_1d", sum_1d, lowt, hit, 7);
  f_1d->SetParNames("N_PR", "N_NP", "f", "mu", "sigma1", "sigma2", "lambda");

  double pt_val[nPtBins], pt_err[nPtBins];

  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.13);

  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    pt_val[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
    pt_err[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);

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
    h_d1d[i_pt]->SetMinimum(f_1d->Eval(-pr_lim));
    //    h_d1d[i_pt]->SetMinimum(h_d1d[i_pt]->GetMaximum()*5e-3);

    TH1F *fh = c->DrawFrame(lowPlot, h_d1d[i_pt]->GetMinimum(), hit, h_d1d[i_pt]->GetMaximum());
    fh->SetXTitle("c#tau (mm)");
    fh->SetYTitle(Form("Events per %.0f #mum", wbin*1000.));
    fh->GetYaxis()->SetTitleOffset(1.8);
    fh->GetYaxis()->SetLabelOffset(0.01);
    fh->SetTitle(h_d1d[i_pt]->GetTitle());

    h_d1d[i_pt]->SetMarkerStyle(20);
    h_d1d[i_pt]->SetMarkerColor(kBlack);
    h_d1d[i_pt]->SetLineColor(kBlack);
    h_d1d[i_pt]->SetMarkerSize(0.75);
    h_d1d[i_pt]->Draw("error same");

    f_1d->SetNpx(300);
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

    c->SaveAs(Form("plots/lifetime/fit_%d.pdf", i_pt));
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
    fp->SetTitle(Form("Lifetime fit pulls (%.1f < p_{T} < %.1f GeV)", ptBins[i_pt], ptBins[i_pt+1]));
  
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

    c->SaveAs(Form("plots/lifetime/pulls_%d.pdf", i_pt));
    c->Clear();

  }
  c->Destructor();
}
