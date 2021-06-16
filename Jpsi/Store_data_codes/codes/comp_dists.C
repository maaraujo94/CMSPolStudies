// plot the non-normalized pT dist of each sample for direct comparison
void comp_dists()
{
  TFile *fData = new TFile("../data18_cos.root");
  TFile *fMC1 = new TFile("../MC18_cos.root");
  TFile *fMC2 = new TFile("../MC18_hpt_cos.root");
  TFile *fMC3 = new TFile("../MC18_vhpt_cos.root");

  TTree *tData = (TTree*)fData->Get("data_cos");
  TTree *tMC1 = (TTree*)fMC1->Get("MC_cos");
  TTree *tMC2 = (TTree*)fMC2->Get("MC_cos");
  TTree *tMC3 = (TTree*)fMC3->Get("MC_cos");

  TH1D* hData1 = new TH1D("hData1", "2018 J/#psi p_{T}", 175, 25, 200);
  TH1D* hData2 = new TH1D("hData2", "2018 J/#psi p_{T}", 175, 25, 200);
  TH1D* hMC1 = new TH1D("hMC1", "2018 J/#psi p_{T}", 175, 25, 200);
  TH1D* hMC2 = new TH1D("hMC2", "2018 J/#psi p_{T}", 175, 25, 200);
  TH1D* hMC3 = new TH1D("hMC3", "2018 J/#psi p_{T}", 175, 25, 200);

  int nEvt;
  double pt, lts;

  nEvt = tData->GetEntries();
  tData->SetBranchAddress("JpsiPt", &pt);
  tData->SetBranchAddress("lts", &lts);
  for(int i = 0; i < nEvt; i++) {
    tData->GetEntry(i);
    if(lts<2.5)
      hData1->Fill(pt);
    else
      hData2->Fill(pt);
  }
  nEvt = tMC1->GetEntries();
  tMC1->SetBranchAddress("JpsiPt", &pt);
  for(int i = 0; i < nEvt; i++) {
    tMC1->GetEntry(i);
    hMC1->Fill(pt);
  }
  nEvt = tMC2->GetEntries();
  tMC2->SetBranchAddress("JpsiPt", &pt);
  for(int i = 0; i < nEvt; i++) {
    tMC2->GetEntry(i);
    hMC2->Fill(pt);
  }
  nEvt = tMC3->GetEntries();
  tMC3->SetBranchAddress("JpsiPt", &pt);
  for(int i = 0; i < nEvt; i++) {
    tMC3->GetEntry(i);
    hMC3->Fill(pt);
  }

  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLogy();
  
  hData1->SetStats(0);
  hData1->SetLineColor(1);
  hData1->GetXaxis()->SetTitle("p_{T} (GeV)");
  hData1->Draw();
  hData2->SetLineStyle(kDashed);
  hData2->SetLineColor(1);
  hData2->Draw("same");
  hMC1->SetLineColor(2);
  hMC1->Draw("same");
  hMC2->SetLineColor(3);
  hMC2->Draw("same");
  hMC3->SetLineColor(4);
  hMC3->Draw("same");

  TLegend *leg = new TLegend(0.6, 0.65, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(hData1, "Data (lts<2.5)", "l");
  leg->AddEntry(hData2, "Data (lts>2.5)", "l");
  leg->AddEntry(hMC1, "MC", "l");
  leg->AddEntry(hMC2, "MC (p_{T}>46GeV)", "l");
  leg->AddEntry(hMC3, "MC (p_{T}>66GeV)", "l");
  leg->Draw();
  
  c->SaveAs("pt_comp.pdf");
}
