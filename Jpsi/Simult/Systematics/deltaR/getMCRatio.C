void getMCRatio()
{
  // get the baseline MC
  TFile *fIn1 = new TFile("../../PR_fit/files/histoStore.root");
  TH2D *mc_1 = (TH2D*)fIn1->Get("mcH_ab");
  mc_1->SetDirectory(0);
  fIn1->Close();

  // get the loose cut MC
  TFile *fIn2 = new TFile("../../../Simult_dR2/PR_fit/files/histoStore.root");
  TH2D *mc_2 = (TH2D*)fIn2->Get("mcH_ab");
  mc_2->SetDirectory(0);
  fIn2->Close();

  // get the tight cut MC
  TFile *fIn3 = new TFile("../../../Simult_dR1/PR_fit/files/histoStore.root");
  TH2D *mc_3 = (TH2D*)fIn3->Get("mcH_ab");
  mc_3->SetDirectory(0);
  fIn3->Close();

  TH2D *ratio_L = (TH2D*)mc_2->Clone("ratio_L");
  ratio_L->Sumw2();
  ratio_L->SetTitle("Cut / Full MC (#DeltaR>0.15)");
  ratio_L->Divide(mc_1);
  TH2D *ratio_T = (TH2D*)mc_3->Clone("ratio_T");
  ratio_T->Sumw2();
  ratio_T->SetTitle("Cut / Full MC (#DeltaR>0.17)");
  ratio_T->Divide(mc_1);

  TFile *fout = new TFile("files/mc_ratio.root", "recreate");
  ratio_L->Write();
  ratio_T->Write();
  fout->Close();
}
