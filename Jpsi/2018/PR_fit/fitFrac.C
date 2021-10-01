// macro to get the f_bkg uncertainty and compare to direct fit results
// then fit the background fraction f_bkg
void fitFrac()
{
  // PART 1: get f_bkg w unc and compare to direct results
  
  // get the fractions obtained directly from the fits
  TFile *inSBo = new TFile("files/mfit.root");
  TGraphErrors *fit_fBG = (TGraphErrors*)inSBo->Get("fit_fBG");
  inSBo->Close();

  // fix f_SB to span 0 to 100
  int nbg = fit_fBG->GetN();
  double *fy = fit_fBG->GetY();
  double *fey = fit_fBG->GetEY();
  for(int i = 0; i < nbg; i++) {
    fy[i] *= 100.;
    fey[i] *= 100.;
  }

  // get the histos from the generation and obtain mean, std dev
  TH1F **h_fbg = new TH1F*[nbg];
  double bgy[nbg], bgey[nbg];
  TFile *inSB = new TFile("files/fbgDists.root");
  for(int i = 0; i < nbg; i++) {
    inSB->GetObject(Form("h_fbg_%d", i), h_fbg[i]);
    bgy[i] = h_fbg[i]->GetMean() * 100.;
    bgey[i] = h_fbg[i]->GetStdDev() * 100.;
  }
  inSB->Close();
  TGraphErrors *fSB = new TGraphErrors(nbg, fit_fBG->GetX(), bgy, fit_fBG->GetEX(), bgey);

  // plotting the comparisons
  TCanvas *c = new TCanvas("", "", 900, 900);

  // f_bg first
  double x_min = fit_fBG->GetX()[0]-fit_fBG->GetEX()[0]-5;
  double x_max = fit_fBG->GetX()[nbg-1]+fit_fBG->GetEX()[nbg-1]+5;
  TH1F *ffbg = c->DrawFrame(x_min, 0, x_max, 15);
  ffbg->SetXTitle("p_{T} (GeV)");
  ffbg->SetYTitle("f_{BG} (%)");
  ffbg->GetYaxis()->SetTitleOffset(1.3);
  ffbg->GetYaxis()->SetLabelOffset(0.01);
  ffbg->SetTitle("2018 f_{BG} vs p_{T}");
  
  fit_fBG->SetMarkerStyle(20);
  fit_fBG->SetMarkerSize(.75);
  fit_fBG->SetMarkerColor(kBlack);
  fit_fBG->SetLineColor(kBlack);
  fit_fBG->Draw("psame");

  fSB->SetMarkerStyle(24);
  fSB->SetMarkerSize(.75);
  fSB->SetMarkerColor(kBlue);
  fSB->SetLineColor(kBlue);
  fSB->Draw("psame");

  c->SaveAs(Form("plots/fBG_comp.pdf"));
  c->Clear();

  // PART 2: fitting f_bkg

  // fit and plot f_SB
  TH1F *fr1 = c->DrawFrame(20., 0, 125, 15);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f_{bkg} (%)");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("2018 f_{bkg}");

  fSB->SetLineColor(kBlack);
  fSB->SetMarkerColor(kBlack);
  fSB->SetMarkerStyle(20);
  fSB->Draw("p");

  // the fit function
  TF1 *f_fit1 = new TF1("fit_SB", "[0]*(1-exp(-[1]*(x-[2])))", 0, 125);
  f_fit1->SetParNames("M", "a", "mu");
  f_fit1->SetParameters(10, 0.01, 1.);
  f_fit1->FixParameter(2, 0);
  f_fit1->SetLineColor(kBlue);
  TFitResultPtr fitres = fSB->Fit("fit_SB", "S");
  f_fit1->SetRange(0, 125);
  f_fit1->Draw("same");

  c->SaveAs("plots/fSB_fit.pdf");
  c->Clear();

  // now to generate an uncertainty band over the 17 pT bins
  double M_v = f_fit1->GetParameter(0), a_v = f_fit1->GetParameter(1);
  double M_e = f_fit1->GetParError(0), a_e = f_fit1->GetParError(1);
  double cov = fitres->GetCovarianceMatrix()[0][1];

  // get binning from the stored data histos
  TFile *infile = new TFile("files/histoStore.root");
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
  TH2D *h_fbkg = new TH2D("h_fbkg", "2018 f_{bkg}", nBinsX, minX, maxX, nBinsY, yBins);
  double ln = 10000;
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    double pt = 0.5*(yBins[i_pt+1]+yBins[i_pt]);
    double xpar[] = {pt};
    
    double par_M[] = {M_v+M_e/ln, a_v, 0};
    double par_a[] = {M_v, a_v+a_e/ln, 0};
    double dM = (f_fit1->EvalPar(xpar, par_M)-f_fit1->Eval(pt))/(M_e/ln);
    double da = (f_fit1->EvalPar(xpar, par_a)-f_fit1->Eval(pt))/(a_e/ln);

    for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
      h_fbkg->SetBinContent(i_cos+1, i_pt+1, f_fit1->Eval(pt));
      h_fbkg->SetBinError(i_cos+1, i_pt+1, sqrt( pow(dM * M_e, 2) + pow(da * a_e, 2) + 2*dM*da*cov ));
    }
  }

  // plotting the 1d projection into pT
  TH1D* h_fbkgpt = h_fbkg->ProjectionY("h_fbkgpd", 1, 1);
  
  h_fbkgpt->SetStats(0);
  h_fbkgpt->SetFillColorAlpha(kBlue, 0.5);
  h_fbkgpt->Draw("e3");
  f_fit1->Draw("same");
  fSB->Draw("p");

  c->SaveAs("plots/fBG_band.pdf");
  c->Clear();
  
  c->Destructor();

  // scale f_bkg down from percentage
  f_fit1->SetParameter(0, f_fit1->GetParameter(0)/100.);
  h_fbkg->Scale(1./100.);
  
  TFile *fout = new TFile("files/bkgFrac.root", "recreate");
  fSB->Write();
  f_fit1->Write();
  h_fbkg->Write();
  fout->Close();
}
