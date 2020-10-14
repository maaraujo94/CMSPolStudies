// code to plot the cut variables for the B->J/psi K sample
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y, ct/cterr
*/

void getCuts_BJK()
{
  // tree for data
  TChain *tree = new TChain("treeS");

  tree->Add("/eos/user/m/maaraujo/BtoPsiK/filtered-all-bupsik-psi25-DT18.root");

  // creating desired vars and setting branch address
  Double_t muPPt, muNPt, muPEta, muNEta, JpsiMass, JpsiPt, JpsiRap, Jpsict, JpsictErr, massB;
  TLorentzVector *mumu_p4 = 0, *muP_p4 = 0, *muM_p4 = 0;
  
  tree->SetBranchAddress("muP_p4", &muP_p4);
  tree->SetBranchAddress("muM_p4", &muM_p4);
  tree->SetBranchAddress("mumu_p4", &mumu_p4);
  tree->SetBranchAddress("ctpv", &Jpsict);
  tree->SetBranchAddress("ctpv_error", &JpsictErr);
  tree->SetBranchAddress("massB", &massB);
  
  // aux vars for reading tree
  int nEvt = tree->GetEntries();
  int perc = nEvt / 100;

  // preparing histograms to be filled
  TH1D *h_muPpT  = new TH1D("h_muPpT",  "data muon pT",            100, 0, 150);
  TH1D *h_muNpT  = new TH1D("h_muNpT",  "muonN pT",                100, 0, 150);
  TH1D *h_muPEta = new TH1D("h_muPEta", "data muon eta",           100, -3, 3);
  TH1D *h_muNEta = new TH1D("h_muNEta", "muonN eta",               100, -3, 3);
  TH1D *h_JMass  = new TH1D("h_JMass",  "data Jpsi mass",          100, 2.9, 3.3);
  TH1D *h_JPt    = new TH1D("h_JPt",    "data Jpsi pT",            100, 0, 200);
  TH1D *h_Jy     = new TH1D("h_Jy",     "data Jpsi y",             100, -2.5, 2.5);
  TH1D *h_Jlts   = new TH1D("h_Jlts",   "data Jpsi lifetime sig",  100, 0, 25);
  TH1D *h_BMass  = new TH1D("h_BMass",  "data B mass",             100, 5, 5.6);

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
    h_Jlts->Fill(abs(Jpsict/JpsictErr));
    h_BMass->Fill(massB);
    
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with data" << endl; 
  }

  TFile *outfile = new TFile("BtoJK_data_cuts.root", "recreate");
  h_muPpT->Write();
  h_muNpT->Write();
  h_muPEta->Write();
  h_muNEta->Write();
  h_JMass->Write();
  h_JPt->Write();
  h_Jy->Write();
  h_Jlts->Write();
  h_BMass->Write();
  outfile->Close();
  
}
