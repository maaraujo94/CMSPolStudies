// code to plot the cut variables for the 2018 Psi(2S) data
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y, ct/cterr
*/

void getCuts_data()
{
  // tree for data
  TChain *tree = new TChain("psi2stree");

  tree->Add("/eos/user/m/maaraujo/Psi2SRun2/filtered-all-psi2s-UL17_18.root");

  // creating desired vars and setting branch address
  Double_t  ct, ctErr, vProb;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
  UInt_t trigger, run;
  
  tree->SetBranchAddress("muP_p4", &muP_p4);
  tree->SetBranchAddress("muM_p4", &muM_p4);
  tree->SetBranchAddress("mumu_p4", &mumu_p4);
  tree->SetBranchAddress("vProb", &vProb);
  tree->SetBranchAddress("trigger", &trigger);
  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("ctpv", &ct);
  tree->SetBranchAddress("ctpv_error", &ctErr);

  // aux vars for reading tree
  int nEvt = tree->GetEntries();
  int perc = nEvt / 100;

  // preparing histograms to be filled
  TH1D *h_muPpT  = new TH1D("h_muPpT",  "data muon pT",            100, 0, 150);
  TH1D *h_muNpT  = new TH1D("h_muNpT",  "muonN pT",                100, 0, 150);
  TH1D *h_muPEta = new TH1D("h_muPEta", "data muon eta",           100, -2, 2);
  TH1D *h_muNEta = new TH1D("h_muNEta", "muonN eta",               100, -2, 2);
  TH1D *h_JMass  = new TH1D("h_JMass",  "data Psi2S mass",         100, 3.35, 4);
  TH1D *h_JPt    = new TH1D("h_JPt",    "data Psi2S pT",           100, 0, 200);
  TH1D *h_Jy     = new TH1D("h_Jy",     "data Psi2S y",            100, -2, 2);
  TH1D *h_Jlts   = new TH1D("h_Jlts",   "data Psi2S lifetime sig", 100, -5, 25);
  TH1D *h_vP     = new TH1D("h_vP",     "data Psi2S vProb",        100, 0, 1);
  
  // reading tree, filling histograms
  for(int i = 0; i < nEvt; i++) {
    tree->GetEntry(i);

    if((trigger&8) == 8 && vProb > 0.01 && run > 313000)
      {
	h_muPpT->Fill(muP_p4->Pt());
	h_muNpT->Fill(muM_p4->Pt());
	h_muPEta->Fill(muP_p4->Eta());
	h_muNEta->Fill(muM_p4->Eta());
	h_JMass->Fill(mumu_p4->M());
	h_JPt->Fill(mumu_p4->Pt());
	h_Jy->Fill(mumu_p4->Rapidity());
	h_Jlts->Fill(ct/ctErr);
	h_vP->Fill(vProb);
      }

    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with data" << endl; 
  }

  TFile *outfile = new TFile("Psi2_data_cuts.root", "recreate");
  h_muPpT->Write();
  h_muNpT->Write();
  h_muPEta->Write();
  h_muNEta->Write();
  h_JMass->Write();
  h_JPt->Write();
  h_Jy->Write();
  h_Jlts->Write();
  h_vP->Write();
  outfile->Close();

}
