// code to plot the cut variables for the 2018 Jpsi MC (high pT)
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y
*/

void getCuts_MC_hpt()
{
  // tree for MC
  TChain *tree = new TChain("psi2stree");

  tree->Add("/eos/user/m/maaraujo/Psi2SRun2/filtered-all-psi2s-mc-LOCAL18-highpt.root");

  // creating desired vars and setting branch address
  Float_t  ct, ctErr;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
  
  tree->SetBranchAddress("muP_p4", &muP_p4);
  tree->SetBranchAddress("muM_p4", &muM_p4);
  tree->SetBranchAddress("mumu_p4", &mumu_p4);
  tree->SetBranchAddress("ctpv", &ct);
  tree->SetBranchAddress("ctpv_error", &ctErr);

  // aux vars for reading tree
  int nEvt = tree->GetEntries();
  int perc = nEvt / 100;

  // preparing histograms to be filled
  TH1D *h_muPpT  = new TH1D("h_muPpT",  "MC (high pT) muon pT",    100, 0, 150);
  TH1D *h_muNpT  = new TH1D("h_muNpT",  "muonN pT",      100, 0, 150);
  TH1D *h_muPEta = new TH1D("h_muPEta", "MC (high pT) muon eta",   100, -2, 2);
  TH1D *h_muNEta = new TH1D("h_muNEta", "muonN eta",     100, -2, 2);
  TH1D *h_JMass  = new TH1D("h_JMass",  "MC (high pT) Psi2S mass", 100, 3.35, 4);
  TH1D *h_JPt    = new TH1D("h_JPt",    "MC (high pT) Psi2S pT",   100, 0, 200);
  TH1D *h_Jy     = new TH1D("h_Jy",     "MC (high pT) Psi2S y",    100, -2, 2);
  TH1D *h_Jlts   = new TH1D("h_Jlts",   "MC (high pT) Psi2S lts",  100, -5, 25);

  // reading tree, filling histograms
  for(int i = 0; i < nEvt; i++) {
    tree->GetEntry(i);

    h_muPpT->Fill(muP_p4->Pt());
    h_muNpT->Fill(muM_p4->Pt());
    h_muPEta->Fill(muP_p4->Eta());
    h_muNEta->Fill(muM_p4->Eta());
    h_JMass->Fill(mumu_p4->M());
    h_JPt->Fill(mumu_p4->Pt());
    h_Jy->Fill(mumu_p4->Rapidity());
    h_Jlts->Fill(ct/ctErr);
 
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with MC" << endl; 
  }

  TFile *outfile = new TFile("Psi2_MC_hpt_cuts.root", "recreate");
  h_muPpT->Write();
  h_muNpT->Write();
  h_muPEta->Write();
  h_muNEta->Write();
  h_JMass->Write();
  h_JPt->Write();
  h_Jy->Write();
  h_Jlts->Write();
  outfile->Close();

}
