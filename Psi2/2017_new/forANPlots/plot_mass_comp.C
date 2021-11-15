void plot_mass_comp()
{
  // prepare binning and histograms for plots
  const int nPtBins = 7;
  double ptBins[nPtBins+1];
  for(int i=0; i<3; i++) ptBins[i] = 7.*i+25.;
  for(int i=0; i<4; i++) ptBins[i+3] = 46.+10.*i;
  ptBins[7] = 120;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // get the PR window mass dists
  TH1D **h_mPR = new TH1D*[nPtBins];
  TFile *fPR = new TFile("../PR_fit/files/mStore.root");
  for(int i = 0; i< nPtBins; i++) {
    h_mPR[i] = (TH1D*)fPR->Get(Form("mH%.0f", ptBins[i]));
    h_mPR[i]->SetDirectory(0);
  }
  fPR->Close();

  // get the NP window mass dists
  TH1D **h_mNP = new TH1D*[nPtBins];
  TFile *fNP = new TFile("../NP_fit/files/mStore.root");
  for(int i = 0; i< nPtBins; i++) {
    h_mNP[i] = (TH1D*)fNP->Get(Form("mH%.0f", ptBins[i]));
    h_mNP[i]->SetDirectory(0);
  }
  fNP->Close();

  int bin_min = h_mPR[0]->FindBin(3.21+1e-5);
  int bin_max = h_mPR[0]->FindBin(3.28-1e-5);
  TCanvas *c = new TCanvas("", "", 900, 900);
  for(int i = 0; i < nPtBins; i++) {
    h_mPR[i]->Scale(1./h_mPR[i]->Integral());
    h_mPR[i]->SetMinimum(0);
    h_mPR[i]->SetMaximum(h_mPR[i]->GetMaximum()*1.7);
    h_mPR[i]->SetStats(0);
    h_mPR[i]->SetLineColor(kBlack);
    h_mPR[i]->SetMarkerColor(kBlack);
    h_mPR[i]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
    h_mPR[i]->GetYaxis()->SetTitle("Events (a.u.)");
    h_mPR[i]->GetYaxis()->SetTitleOffset(1.5);
    h_mPR[i]->SetMarkerStyle(20);
    h_mPR[i]->SetMarkerSize(0.75);
    h_mPR[i]->Draw("histo");

    h_mNP[i]->Scale(1./h_mNP[i]->Integral());
    h_mNP[i]->SetLineColor(kBlue);
    h_mNP[i]->SetMarkerColor(kBlue);
    h_mNP[i]->SetMarkerStyle(20);
    h_mNP[i]->SetMarkerSize(0.75);
    h_mNP[i]->Draw("histo same");

    TLegend *leg = new TLegend(0.65, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(h_mPR[i], "PR window", "pl");
    leg->AddEntry(h_mNP[i], "NP window", "pl");
    leg->Draw();

    c->SaveAs(Form("plots/mass/comp_pt%d.pdf", i));
    c->Clear();
  }
}
