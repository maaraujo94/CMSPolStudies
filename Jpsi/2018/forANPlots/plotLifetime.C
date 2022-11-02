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
void plotLifetime()
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
    h_d1d[i]->SetTitle(Form("2018 data c#tau (%.1f < p_{T} < %.1f GeV)", ptBins[i], ptBins[i+1]));
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

  int i_pt = 14;

  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.13);
  
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

  
  h_d1d[i_pt]->SetMaximum(6000);
  h_d1d[i_pt]->SetMinimum(0);

  TH1F *fh = c->DrawFrame(lowPlot, h_d1d[i_pt]->GetMinimum(), hit, h_d1d[i_pt]->GetMaximum());
  fh->SetXTitle("c#tau (mm)");
  fh->SetYTitle(Form("Events per %.0f #mum", wbin*1000.));
  fh->GetYaxis()->SetTitleOffset(1.8);
  fh->GetYaxis()->SetLabelOffset(0.01);
  fh->SetTitle("");

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
  //fres->Draw("lsame");
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

  TLatex lc;
  lc.SetTextSize(0.04);
  lc.DrawLatex(0.2, 4500, "2018 data");
  lc.DrawLatex(0.2, 4000, Form("%.0f < p_{T} < %.0f GeV", ptBins[i_pt], ptBins[i_pt+1]));
  lc.Draw();
  
  c->SaveAs(Form("plots/ltfit.pdf"));
  c->Clear();
  c->Destructor();
}
