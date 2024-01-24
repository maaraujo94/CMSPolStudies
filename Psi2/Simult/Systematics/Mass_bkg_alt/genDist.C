// macro to generate the sideband costh dists in the final binning, with unc

void genDist()
{
  // get binning from the stored data histos
  TFile *infile = new TFile("../../PR_fit/files/histoStore.root");
  TH2D *h_LSB = (TH2D*)infile->Get("PRLH");
  TH2D *h_RSB = (TH2D*)infile->Get("PRRH");
  h_LSB->SetDirectory(0);
  h_RSB->SetDirectory(0);
  infile->Close();

  int nBinsX = h_LSB->GetNbinsX(), nBinsY = h_LSB->GetNbinsY();
  const double *yBins = h_LSB->GetYaxis()->GetXbins()->GetArray();
  double minX = h_LSB->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_LSB->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/nBinsX;

    // get fL from storage
  TFile *infL = new TFile("files/store_fL.root");
  double *fL = ((TGraphErrors*)infL->Get("g_fL"))->GetY();
  infL->Close();

  // get the 1d sb histos and scale to nr events
  double n_s[nBinsY];
  TH1D **h_LSB1d = new TH1D*[nBinsY];
  TH1D **h_RSB1d = new TH1D*[nBinsY];
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    h_LSB1d[i_pt] = h_LSB->ProjectionX(Form("h_LSB_%d", i_pt), i_pt+1, i_pt+1);
    h_RSB1d[i_pt] = h_RSB->ProjectionX(Form("h_RSB_%d", i_pt), i_pt+1, i_pt+1);
    n_s[i_pt] = fL[i_pt] * h_LSB1d[i_pt]->GetEntries() + (1.-fL[i_pt])*h_RSB1d[i_pt]->GetEntries();
    h_LSB1d[i_pt]->Scale(1./h_LSB1d[i_pt]->GetEntries());
    h_RSB1d[i_pt]->Scale(1./h_RSB1d[i_pt]->GetEntries());
  }

  // get the sideband histos by summing with proportion fL
  TH1D **h_SB = new TH1D*[nBinsY];
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    h_SB[i_pt] = new TH1D(Form("h_SB_%d", i_pt), Form("bkg |cos#theta| (%.0f < p_{T} < %.0f GeV)", yBins[i_pt], yBins[i_pt+1]), nBinsX, minX, maxX);
   
    h_SB[i_pt]->Sumw2();
    h_SB[i_pt]->Add(h_LSB1d[i_pt], h_RSB1d[i_pt], fL[i_pt], 1.-fL[i_pt]);

    h_SB[i_pt]->Scale(n_s[i_pt]);
    h_LSB1d[i_pt]->Scale(n_s[i_pt]);
    h_RSB1d[i_pt]->Scale(n_s[i_pt]);
  }

  cout << "all SB histos filled" << endl;
  
  TFile *fout = new TFile("files/bkgCosModel.root", "recreate");
  for(int i = 0; i < nBinsY; i++) {
    h_SB[i]->Write();
  }
  fout->Close();

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);
  
  h_SB[2]->SetStats(0);
  h_SB[2]->SetLineColor(kGreen+1);
  h_SB[2]->SetMinimum(0);
  h_SB[2]->SetMaximum(h_SB[2]->GetMaximum()*1.7);
  h_SB[2]->GetXaxis()->SetTitle("|cos#theta|");
  h_SB[2]->Draw();
  
  h_LSB1d[2]->SetLineColor(kBlack);
  h_LSB1d[2]->Draw("same");
  h_RSB1d[2]->SetLineColor(kBlue);
  h_RSB1d[2]->Draw("same");

  TLegend *leg = new TLegend(0.77, 0.65, 0.97, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_LSB1d[2], "LSB", "l");
  leg->AddEntry(h_RSB1d[2], "RSB", "l");
  leg->AddEntry(h_SB[2], "bkg", "l");
  leg->Draw();

  c->SaveAs("plots/SB_base_full_2.pdf");
  c->Clear();

  h_SB[5]->SetStats(0);
  h_SB[5]->SetLineColor(kGreen+1);
  h_SB[5]->SetMinimum(0);
  h_SB[5]->SetMaximum(h_SB[5]->GetMaximum()*1.7);
  h_SB[5]->GetXaxis()->SetTitle("|cos#theta|");
  h_SB[5]->Draw();
  
  h_LSB1d[5]->SetLineColor(kBlack);
  h_LSB1d[5]->Draw("same");
  h_RSB1d[5]->SetLineColor(kBlue);
  h_RSB1d[5]->Draw("same");

  leg->Draw();

  c->SaveAs("plots/SB_base_full_5.pdf");
  c->Clear();
  c->Destructor();
}
