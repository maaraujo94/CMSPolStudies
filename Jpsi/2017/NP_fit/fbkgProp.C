// macro to propagate the uncertainties to f_bkg analytically
// then fit the function and get the unc band of that

// background function
double bkg_exp(double m, double p1, double p2)
{
  return p1 * exp( - m / p2 );
}

void fbkgProp()
{
  // PART 1: get the uncertainty of the f_bkg
  
  // get all needed input
  // mass region ranges
  double m_min[] = {2.94, 3.0, 3.21};
  double m_max[] = {2.95, 3.2, 3.26};

  // bkg function
  TF1 *f_exp = new TF1("f_exp", "bkg_exp(x,[0],[1])", m_min[0], m_max[2]);
  f_exp->SetParNames("NB", "lambda");
  
  // get fit parameters - need to know which params are being used
  TFile *inBG = new TFile("files/mfit.root");
  int n_par = 2;
  TGraphErrors** g_par = new TGraphErrors*[n_par];
  for(int i = 0; i < n_par; i++) {
    inBG->GetObject(Form("fit_%s", f_exp->GetParName(i)), g_par[i]);
  }
  TFitResult* fitres = (TFitResult*)inBG->Get("fitres");
  inBG->Close();

  // get pT binning from the fit pars
  int n_pt = g_par[0]->GetN();
  double ptBins[n_pt+1];
  for(int i = 0; i < n_pt; i++) {
    ptBins[i] = g_par[0]->GetX()[i]-g_par[0]->GetEX()[i];
  }
  ptBins[n_pt] = g_par[0]->GetX()[n_pt-1]+g_par[0]->GetEX()[n_pt-1];
  
  // prepare mass histograms
  TH1D **h_d1d = new TH1D*[n_pt];
  TFile *fin = new TFile("files/mStore.root");
  for(int ip = 0; ip < n_pt; ip++) {
    fin->GetObject(Form("mH%.0f", ptBins[ip]), h_d1d[ip]);
    h_d1d[ip]->SetDirectory(0);
  }
  fin->Close();

  // fbkg = integral / evt_all (in signal region)
  TH1D *h_fbkg = new TH1D("h_fbkg", "2017 f_{bkg}", n_pt, ptBins);
  double ln = 10000;
  const int n_p = 2;
  double fit_v[n_p], dpar[n_p];
  double cov[n_p][n_p];

  for(int i_pt = 0; i_pt < n_pt; i_pt++) {
    // N_histo has no uncertainty: just histogram integral
    double min_bin = h_d1d[i_pt]->GetXaxis()->FindBin(m_min[1]+1e-6);
    double max_bin = h_d1d[i_pt]->GetXaxis()->FindBin(m_max[1]-1e-6);
    double evt_all = h_d1d[i_pt]->Integral(min_bin, max_bin, "width");

    // define input parameters for integral - varies by pT bin
    for(int i = 0; i < n_p; i++) {
      fit_v[i] = g_par[i]->GetY()[i_pt];
      for(int j = 0; j < n_p; j++) {
	cov[i][j] = fitres->GetCovarianceMatrix()[(7+i)*n_pt+i_pt][(7+j)*n_pt+i_pt];
      }
    }
    // integral is a function of NB and lambda
    f_exp->SetParameters(fit_v);
    double fv = f_exp->Integral(m_min[1], m_max[1]), fe = 0;

    // get the function deviation for each parameter
    for(int i = 0; i < n_p; i++) {
      f_exp->SetParameter(i, fit_v[i] + sqrt(cov[i][i])/ln);
      dpar[i] = (f_exp->Integral(m_min[1], m_max[1])-fv)/(sqrt(cov[i][i])/ln);
      f_exp->SetParameter(i, fit_v[i]);
    }
    
    for(int i = 0; i < n_p; i++) 
      for(int j = 0; j < n_p; j++) {
	fe += dpar[i]*dpar[j]*cov[i][j];
      }
    fe = sqrt(fe);

    // fill pT bin
    h_fbkg->SetBinContent(i_pt+1, fv/evt_all*100.);
    h_fbkg->SetBinError(i_pt+1, fe/evt_all*100.);
  }

  // plotting in pT
  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1F *fr1 = c->DrawFrame(ptBins[0]-5, 0, ptBins[n_pt]+5, 15);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f_{bkg} (%)");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("2017 f_{bkg} vs p_{T}");

  h_fbkg->SetStats(0);
  h_fbkg->SetMarkerStyle(20);
  h_fbkg->SetMarkerSize(.5);
  h_fbkg->SetMarkerColor(kBlack);
  h_fbkg->SetLineColor(kBlack);
  h_fbkg->Draw("e1 same");
  
  c->SaveAs("plots/fBG_unc.pdf");
 
  // PART 2: fitting f_bkg

  // the fit function
  TF1 *f_fit1 = new TF1("fit_SB", "[0]", 0, 125);
  f_fit1->SetParName(0, "f_bkg");
  f_fit1->SetParameter(0, 7);
  f_fit1->SetLineColor(kBlue);
  TFitResultPtr fitbgres = h_fbkg->Fit("fit_SB", "S0");
  f_fit1->SetRange(0, 125);
  f_fit1->Draw("same");

  c->SaveAs("plots/fBG_fit.pdf");
  c->Clear();

  // PART 3: generate f_bkg histo

  // now to generate an uncertainty band over the 17 pT bins
  // constant so it's just the uncertainty on f_bkg
  // get binning from the stored data histos
  TFile *infile = new TFile("../PR_fit/files/histoStore.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("dataH_ab"), hist);
  hist->SetDirectory(0);
  infile->Close();

  // get the binning
  int nBinsX = hist->GetNbinsX(), nBinsY = hist->GetNbinsY();
  const double *yBins = hist->GetYaxis()->GetXbins()->GetArray();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/nBinsX;
  
  // f_bkg(pT) but generating 2d map so it's easier to apply uncertainties
  TH2D *h_fbkg2d = new TH2D("h_fbkg2d", "2017 f_{bkg}", nBinsX, minX, maxX, nBinsY, yBins);
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
      h_fbkg2d->SetBinContent(i_cos+1, i_pt+1, f_fit1->GetParameter(0));
      h_fbkg2d->SetBinError(i_cos+1, i_pt+1, f_fit1->GetParError(0));
    }
  }

  // plotting the 1d projection into pT
  TH1D* h_fbkgpt = h_fbkg2d->ProjectionY("h_fbkgpd", 1, 1);

  h_fbkgpt->SetMinimum(0);
  h_fbkgpt->SetMaximum(15);
  h_fbkgpt->GetXaxis()->SetTitle("p_{T} (GeV)");
  h_fbkgpt->GetYaxis()->SetTitle("f_{bkg} (%)");
  h_fbkgpt->GetYaxis()->SetTitleOffset(1.3);
  h_fbkgpt->GetYaxis()->SetLabelOffset(0.01);
  h_fbkgpt->SetTitle("2017 f_{bkg} vs p_{T}");
  h_fbkgpt->SetStats(0);
  h_fbkgpt->SetFillColorAlpha(kBlue, 0.5);
  h_fbkgpt->Draw("e3");
  f_fit1->Draw("same");
  h_fbkg->Draw("e0 same");
  
  c->SaveAs("plots/fBG_band.pdf");
  c->Clear();
  c->Destructor();

  // scale fractions down from percentage
  h_fbkg->Scale(1./100.);
  f_fit1->SetParameter(0, f_fit1->GetParameter(0)/100.);
  h_fbkg2d->Scale(1./100.);

  TFile *fout = new TFile("files/bkgFrac.root", "recreate");
  h_fbkg->SetName("fbkg_unc");
  h_fbkg->Write();
  f_fit1->Write();
  h_fbkg2d->SetName("h_fbkg");
  h_fbkg2d->Write();
  fout->Close();

}
