// code to plot the cut variables for the 2017/18 Psi(2S) MC
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y
*/

void getCuts_MC()
{
  // preparing histograms to be filled
  //2017
  TH1D *h7_muPpT  = new TH1D("h7_muPpT",  "2017 MC muon pT",            100, 0, 150);
  TH1D *h7_muNpT  = new TH1D("h7_muNpT",  "muonN pT",                     100, 0, 150);
  TH1D *h7_muPEta = new TH1D("h7_muPEta", "2017 MC muon eta",           100, -2, 2);
  TH1D *h7_muNEta = new TH1D("h7_muNEta", "muonN eta",                    100, -2, 2);
  TH1D *h7_JMass  = new TH1D("h7_JMass",  "2017 MC Psi2S mass",         100, 3.35, 4);
  TH1D *h7_JPt    = new TH1D("h7_JPt",    "2017 MC Psi2S pT",           100, 0, 200);
  TH1D *h7_Jy     = new TH1D("h7_Jy",     "2017 MC Psi2S y",            100, -2, 2);
  TH1D *h7_Jlts   = new TH1D("h7_Jlts",   "2017 MC Psi2S lifetime sig", 100, -5, 25);

  // 2018
  TH1D *h8_muPpT  = new TH1D("h8_muPpT",  "2018 MC muon pT",            100, 0, 150);
  TH1D *h8_muNpT  = new TH1D("h8_muNpT",  "muonN pT",                     100, 0, 150);
  TH1D *h8_muPEta = new TH1D("h8_muPEta", "2018 MC muon eta",           100, -2, 2);
  TH1D *h8_muNEta = new TH1D("h8_muNEta", "muonN eta",                    100, -2, 2);
  TH1D *h8_JMass  = new TH1D("h8_JMass",  "2018 MC Psi2S mass",         100, 3.35, 4);
  TH1D *h8_JPt    = new TH1D("h8_JPt",    "2018 MC Psi2S pT",           100, 0, 200);
  TH1D *h8_Jy     = new TH1D("h8_Jy",     "2018 MC Psi2S y",            100, -2, 2);
  TH1D *h8_Jlts   = new TH1D("h8_Jlts",   "2018 MC Psi2S lifetime sig", 100, -5, 25);

  Float_t  ct, ctErr;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;

  // 2017 tree
  TFile *fin7 = new TFile("/eos/user/m/maaraujo/Psi2SRun2/filtered-all-psi2s-mc-LOCAL17.root");
  TTree *tree7 = (TTree*)fin7->Get("psi2stree");

  // setting branch address  
  tree7->SetBranchAddress("muP_p4", &muP_p4);
  tree7->SetBranchAddress("muM_p4", &muM_p4);
  tree7->SetBranchAddress("mumu_p4", &mumu_p4);
  tree7->SetBranchAddress("ctpv", &ct);
  tree7->SetBranchAddress("ctpv_error", &ctErr);

  // aux vars for reading tree
  int nEvt = tree7->GetEntries();
  int perc = nEvt / 100;

  // reading tree, filling histograms
  for(int i = 0; i < nEvt; i++) {
    tree7->GetEntry(i);

    h7_muPpT->Fill(muP_p4->Pt());
    h7_muNpT->Fill(muM_p4->Pt());
    h7_muPEta->Fill(muP_p4->Eta());
    h7_muNEta->Fill(muM_p4->Eta());
    h7_JMass->Fill(mumu_p4->M());
    h7_JPt->Fill(mumu_p4->Pt());
    h7_Jy->Fill(mumu_p4->Rapidity());
    h7_Jlts->Fill(ct/ctErr);
 
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with MC" << endl; 
  }

  fin7->Close();
  
  TFile *fout7 = new TFile("Psi2_17_MC_cuts.root", "recreate");
  h7_muPpT->Write();
  h7_muNpT->Write();
  h7_muPEta->Write();
  h7_muNEta->Write();
  h7_JMass->Write();
  h7_JPt->Write();
  h7_Jy->Write();
  h7_Jlts->Write();
  fout7->Close();

    // 2018 tree
  TFile *fin8 = new TFile("/eos/user/m/maaraujo/Psi2SRun2/filtered-all-psi2s-mc-LOCAL18.root");
  TTree *tree8 = (TTree*)fin8->Get("psi2stree");

  // setting branch address  
  tree8->SetBranchAddress("muP_p4", &muP_p4);
  tree8->SetBranchAddress("muM_p4", &muM_p4);
  tree8->SetBranchAddress("mumu_p4", &mumu_p4);
  tree8->SetBranchAddress("ctpv", &ct);
  tree8->SetBranchAddress("ctpv_error", &ctErr);

  // aux vars for reading tree
  nEvt = tree8->GetEntries();
  perc = nEvt / 100;

  // reading tree, filling histograms
  for(int i = 0; i < nEvt; i++) {
    tree8->GetEntry(i);

    h8_muPpT->Fill(muP_p4->Pt());
    h8_muNpT->Fill(muM_p4->Pt());
    h8_muPEta->Fill(muP_p4->Eta());
    h8_muNEta->Fill(muM_p4->Eta());
    h8_JMass->Fill(mumu_p4->M());
    h8_JPt->Fill(mumu_p4->Pt());
    h8_Jy->Fill(mumu_p4->Rapidity());
    h8_Jlts->Fill(ct/ctErr);
 
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with MC" << endl; 
  }

  fin8->Close();

  TFile *fout8 = new TFile("Psi2_18_MC_cuts.root", "recreate");
  h8_muPpT->Write();
  h8_muNpT->Write();
  h8_muPEta->Write();
  h8_muNEta->Write();
  h8_JMass->Write();
  h8_JPt->Write();
  h8_Jy->Write();
  h8_Jlts->Write();
  fout8->Close();

}
