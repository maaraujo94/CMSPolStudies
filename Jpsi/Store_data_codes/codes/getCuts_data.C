// code to plot the cut variables for the 2018 Jpsi data
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y, ct/cterr
*/

void getCuts_data()
{
  // preparing histograms to be filled
  //2017
  TH1D *h7_muPpT  = new TH1D("h7_muPpT",  "2017 data muon pT",         100, 0, 150);
  TH1D *h7_muNpT  = new TH1D("h7_muNpT",  "muonN pT",                  100, 0, 150);
  TH1D *h7_muPEta = new TH1D("h7_muPEta", "2017 data muon eta",        100, -3, 3);
  TH1D *h7_muNEta = new TH1D("h7_muNEta", "muonN eta",                 100, -3, 3);
  TH1D *h7_JMass  = new TH1D("h7_JMass",  "2017 data dimuon mass",     100, 2.9, 3.3);
  TH1D *h7_JPt    = new TH1D("h7_JPt",    "2017 data dimuon pT",       100, 0, 200);
  TH1D *h7_Jy     = new TH1D("h7_Jy",     "2017 data dimuon y",        100, -2.5, 2.5);
  TH1D *h7_Jlt    = new TH1D("h7_Jlt",    "2017 data dimuon lifetime", 100, -0.01, 0.05);
  TH1D *h7_vP     = new TH1D("h7_vP",     "2017 data dimuon vProb",    100, 0, 1);

  // 2018
  TH1D *h8_muPpT  = new TH1D("h8_muPpT",  "2018 data muon pT",         100, 0, 150);
  TH1D *h8_muNpT  = new TH1D("h8_muNpT",  "muonN pT",                  100, 0, 150);
  TH1D *h8_muPEta = new TH1D("h8_muPEta", "2018 data muon eta",        100, -3, 3);
  TH1D *h8_muNEta = new TH1D("h8_muNEta", "muonN eta",                 100, -3, 3);
  TH1D *h8_JMass  = new TH1D("h8_JMass",  "2018 data dimuon mass",     100, 2.9, 3.3);
  TH1D *h8_JPt    = new TH1D("h8_JPt",    "2018 data dimuon pT",       100, 0, 200);
  TH1D *h8_Jy     = new TH1D("h8_Jy",     "2018 data dimuon y",        100, -2.5, 2.5);
  TH1D *h8_Jlt    = new TH1D("h8_Jlt",    "2018 data dimuon lifetime", 100, -0.01, 0.05);
  TH1D *h8_vP     = new TH1D("h8_vP",     "2018 data dimuon vProb",    100, 0, 1);

  Double_t  ct, vProb;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
  UInt_t trigger;

  // 2017 tree
  TFile *fin7 = new TFile("/eos/user/m/maaraujo/JpsiRun2/Data/filtered-all-psi-UL17.root");
  TTree *tree7 = (TTree*)fin7->Get("jpsitree");

  // setting branch address  
  tree7->SetBranchAddress("muP_p4", &muP_p4);
  tree7->SetBranchAddress("muM_p4", &muM_p4);
  tree7->SetBranchAddress("mumu_p4", &mumu_p4);
  tree7->SetBranchAddress("vProb", &vProb);
  tree7->SetBranchAddress("trigger", &trigger);
  tree7->SetBranchAddress("ctpv", &ct);

  // aux vars for reading tree
  int nEvt = tree7->GetEntries();
  int perc = nEvt / 100;
 
  // reading tree, filling histograms
  for(int i = 0; i < nEvt; i++) {
    tree7->GetEntry(i);

    if((trigger&16) == 16 && vProb > 0.01)
      {
	h7_muPpT->Fill(muP_p4->Pt());
	h7_muNpT->Fill(muM_p4->Pt());
	h7_muPEta->Fill(muP_p4->Eta());
	h7_muNEta->Fill(muM_p4->Eta());
	h7_JMass->Fill(mumu_p4->M());
	h7_JPt->Fill(mumu_p4->Pt());
	h7_Jy->Fill(mumu_p4->Rapidity());
	h7_Jlt->Fill(ct);
	h7_vP->Fill(vProb);
      }

    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2017 data" << endl; 
  }
  cout << endl;
  
  fin7->Close();
  
  TFile *fout7 = new TFile("data17_cuts.root", "recreate");
  h7_muPpT->Write();
  h7_muNpT->Write();
  h7_muPEta->Write();
  h7_muNEta->Write();
  h7_JMass->Write();
  h7_JPt->Write();
  h7_Jy->Write();
  h7_Jlt->Write();
  h7_vP->Write();
  fout7->Close();

  // 2018 tree
  TFile *fin8 = new TFile("/eos/user/m/maaraujo/JpsiRun2/Data/filtered-all-psi-UL18.root");
  TTree *tree8 = (TTree*)fin8->Get("jpsitree");

  // setting branch address  
  tree8->SetBranchAddress("muP_p4", &muP_p4);
  tree8->SetBranchAddress("muM_p4", &muM_p4);
  tree8->SetBranchAddress("mumu_p4", &mumu_p4);
  tree8->SetBranchAddress("vProb", &vProb);
  tree8->SetBranchAddress("trigger", &trigger);
  tree8->SetBranchAddress("ctpv", &ct);
  
  // aux vars for reading tree
  nEvt = tree8->GetEntries();
  perc = nEvt / 100;
 
  // reading tree, filling histograms
  for(int i = 0; i < nEvt; i++) {
    tree8->GetEntry(i);

    if((trigger&16) == 16 && vProb > 0.01)
      {
	h8_muPpT->Fill(muP_p4->Pt());
	h8_muNpT->Fill(muM_p4->Pt());
	h8_muPEta->Fill(muP_p4->Eta());
	h8_muNEta->Fill(muM_p4->Eta());
	h8_JMass->Fill(mumu_p4->M());
	h8_JPt->Fill(mumu_p4->Pt());
	h8_Jy->Fill(mumu_p4->Rapidity());
	h8_Jlt->Fill(ct);
	h8_vP->Fill(vProb);
      }

    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2018 data" << endl; 
  }
  cout << endl;
  
  fin8->Close();
  
  TFile *fout8 = new TFile("data18_cuts.root", "recreate");
  h8_muPpT->Write();
  h8_muNpT->Write();
  h8_muPEta->Write();
  h8_muNEta->Write();
  h8_JMass->Write();
  h8_JPt->Write();
  h8_Jy->Write();
  h8_Jlt->Write();
  h8_vP->Write();
  fout8->Close();

}
