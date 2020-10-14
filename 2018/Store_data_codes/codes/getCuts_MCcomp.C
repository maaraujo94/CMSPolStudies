// code to compare the cut variables for the 2018 Jpsi MC in both samples (low and high pT - very high pT not included)
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y
*/

void getCuts_MCcomp()
{
  string name[2] = {"filtered-all-psi-mc-LOCAL18.root", "filtered-all-psi-mc-LOCAL18-highpt.root"};

  for(int i = 0; i < 2; i++)
    {
      // tree for MC
      TChain *tree = new TChain("jpsitree");

      tree->Add(Form("/eos/user/m/maaraujo/JpsiRun2/MC/%s", name[i].c_str()));

      // creating desired vars and setting branch address
      TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
  
      tree->SetBranchAddress("muP_p4", &muP_p4);
      tree->SetBranchAddress("muM_p4", &muM_p4);
      tree->SetBranchAddress("mumu_p4", &mumu_p4);
      
      // aux vars for reading tree
      int nEvt = tree->GetEntries();
      int perc = nEvt / 100;
      
      // preparing histograms to be filled
      TH1D *h_muPpT  = new TH1D(Form("h_muPpT_%d", i),  "data muon pT",            100, 0, 60);
      TH1D *h_muNpT  = new TH1D(Form("h_muNpT_%d", i),  "muonN pT",                100, 0, 60);
      TH1D *h_muPEta = new TH1D(Form("h_muPEta_%d", i), "data muon eta",           100, -1.5, 1.5);
      TH1D *h_muNEta = new TH1D(Form("h_muNEta_%d", i), "muonN eta",               100, -1.5, 1.5);
      TH1D *h_JMass  = new TH1D(Form("h_JMass_%d", i),  "data Jpsi mass",          100, 2.9, 3.3);
      TH1D *h_JPt    = new TH1D(Form("h_JPt_%d", i),    "data Jpsi pT",            100, 50, 55);
      TH1D *h_Jy     = new TH1D(Form("h_Jy_%d", i),     "data Jpsi y",             100, -1.5, 1.5);
      
      // reading tree, filling histograms
      for(int i = 0; i < nEvt; i++) {
	tree->GetEntry(i);
	if(mumu_p4->Pt() > 50 && mumu_p4->Pt() < 55)
	  {
	    
	    h_muPpT->Fill(muP_p4->Pt());
	    h_muNpT->Fill(muM_p4->Pt());
	    h_muPEta->Fill(muP_p4->Eta());
	    h_muNEta->Fill(muM_p4->Eta());
	    h_JMass->Fill(mumu_p4->M());
	    h_JPt->Fill(mumu_p4->Pt());
	    h_Jy->Fill(mumu_p4->Rapidity());
	    
	  }
	if((i+1)%perc == 0) cout << (i+1)/perc << "% done with MC" << endl; 
      }
      
      TFile *outfile = new TFile("Jpsi_MC_comp.root", "update");
      h_muPpT->Write();
      h_muNpT->Write();
      h_muPEta->Write();
      h_muNEta->Write();
      h_JMass->Write();
      h_JPt->Write();
      h_Jy->Write();
      outfile->Close();
    } 
}
