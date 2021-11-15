// macro to generate the sideband costh dists in the final binning, with unc

void genDist()
{
  // get binning from the stored data histos
  TFile *infile = new TFile("../PR_fit/files/histoStore.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("dataH_ab"), hist);
  hist->SetDirectory(0);
  infile->Close();

  int nBinsX = hist->GetNbinsX(), nBinsY = hist->GetNbinsY();
  const double *yBins = hist->GetYaxis()->GetXbins()->GetArray();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/nBinsX;

  // get fit parameters from storage
  TFile *infL = new TFile("files/store_fL.root");
  double fL = ((TGraphErrors*)infL->Get("g_fL"))->GetY()[0];
  infL->Close();
  // get LSB and RSB histograms
  TFile *inSB = new TFile("files/bkgHist.root");
  TH2D *h_LSB = (TH2D*)inSB->Get("dataH0_ab");
  TH2D *h_RSB = (TH2D*)inSB->Get("dataH1_ab");
  h_LSB->SetDirectory(0);
  h_RSB->SetDirectory(0);
  inSB->Close();

  // get the 1d sb histos and scale to nr events
  TH1D **h_LSB1d = new TH1D*[nBinsY];
  TH1D **h_RSB1d = new TH1D*[nBinsY];
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    h_LSB1d[i_pt] = h_LSB->ProjectionX(Form("h_LSB_%d", i_pt), i_pt+1, i_pt+1);
    h_LSB1d[i_pt]->Scale(1./h_LSB1d[i_pt]->GetEntries());
    h_RSB1d[i_pt] = h_RSB->ProjectionX(Form("h_RSB_%d", i_pt), i_pt+1, i_pt+1);
    h_RSB1d[i_pt]->Scale(1./h_RSB1d[i_pt]->GetEntries());
  }

  // get the sideband histos by summing with proportion fL
  TH1D **h_SB = new TH1D*[nBinsY];
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    h_SB[i_pt] = new TH1D(Form("h_SB_%d", i_pt), Form("SB cos#theta (%.0f < p_{T} < %.0f GeV)", yBins[i_pt], yBins[i_pt+1]), nBinsX, minX, maxX);
    
    h_SB[i_pt]->Sumw2();
    h_SB[i_pt]->Add(h_LSB1d[i_pt], h_RSB1d[i_pt], fL, 1.-fL);
  }

  cout << "all SB histos filled" << endl;
  
  TFile *fout = new TFile("files/bkgCosModel.root", "recreate");
  for(int i = 0; i < nBinsY; i++) {
    h_SB[i]->Write();
  }
  fout->Close();

  double norm = 1.;
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  h_SB[0]->SetStats(0);
  h_SB[0]->SetLineColor(kGreen+1);
  h_SB[0]->Scale(norm);
  h_SB[0]->SetMinimum(0);
  //h_SB[0]->SetMaximum(cthRSB->Eval(0.95)*1.2);
  h_SB[0]->GetXaxis()->SetTitle("|cos#theta|");
  h_SB[0]->Draw();
  
  h_LSB1d[0]->SetLineColor(kBlack);
  h_LSB1d[0]->Draw("same");
  h_RSB1d[0]->SetLineColor(kBlue);
  h_RSB1d[0]->Draw("same");

  TLegend *leg = new TLegend(0.7, 0.2, 0.9, 0.45);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_SB[0], "SR", "l");
  leg->AddEntry(h_LSB1d[0], "LSB", "l");
  leg->AddEntry(h_RSB1d[0], "RSB", "l");
  leg->Draw();

  c->SaveAs("plots/SB_base.pdf");
  c->Clear();
  c->Destructor();


}
