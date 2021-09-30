// macro to fit the background fractions fNP and fSB
void fitFrac()
{ 
  // get background fractions
  TFile *fin1 = new TFile("files/fbkg_unc.root");
  TGraphErrors *fSB = (TGraphErrors*)fin1->Get("f_bkg");
  fin1->Close();

  // fit and plot f_SB
  TCanvas *c = new TCanvas("", "", 900, 900);

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

  TF1 *f_fit1 = new TF1("fit_SB", "[0]*(1-exp(-[1]*(x-[2])))", 0, 125);
  f_fit1->SetParNames("M", "a", "mu");
  f_fit1->SetParameters(10, 0.01, 1.);
  f_fit1->FixParameter(2, 0);
  f_fit1->SetLineColor(kBlue);
  TFitResultPtr fitres = fSB->Fit("fit_SB", "VS");
  f_fit1->SetRange(0, 125);
  f_fit1->Draw("same");

  c->SaveAs("fSB_fit.pdf");
  c->Clear();


  //f_fit1->SetParameter(0, f_fit1->GetParameter(0)/100.);

  // now to generate an uncertainty band over the 17 pT bins
  double M_v = f_fit1->GetParameter(0), a_v = f_fit1->GetParameter(1);
  double M_e = f_fit1->GetParError(0), a_e = f_fit1->GetParError(1);
  double cov = fitres->GetCovarianceMatrix()[0][1];

  // get binning from the stored data histos
  TFile *infile = new TFile("../../PR_fit/files/histoStore.root");
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

  TH1D* h_fbkgpt = h_fbkg->ProjectionY("h_fbkgpd", 1, 1);
  
  h_fbkgpt->SetStats(0);
  h_fbkgpt->SetFillColorAlpha(kBlue, 0.5);
  h_fbkgpt->Draw("e3");
  f_fit1->Draw("same");
  fSB->Draw("p");

  c->SaveAs("band_comp.pdf");
  c->Clear();
  
  c->Destructor();

  f_fit1->SetParameter(0, f_fit1->GetParameter(0)/100.);
  h_fbkg->Scale(1./100.);
  
  TFile *fout = new TFile("files/bkgFrac.root", "recreate");
  fSB->Write();
  f_fit1->Write();
  h_fbkg->Write();
  fout->Close();
}
