// code to compare the cut variables for the 2018 Jpsi MC in both samples (low and high pT - very high pT not included)
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y
*/

void getCuts_MCOldvNew()
{
  string name[2][2] = {{"filtered-all-psi-mc-LOCAL17.root", "filtered-all-psi-mc-UL17-25_46v1.root"},
		       {"filtered-all-psi-mc-LOCAL18.root", "filtered-all-psi-mc-UL18-25_46v2.root"}};
  string lbl[] = {"old", "new"};
  TFile *outE = new TFile("Jpsi_MC_comp.root", "recreate");
  outE->Close();
      
  // run for 2017 and 2018
  for(int i_y = 0; i_y < 2; i_y++) { 
    
    for(int i = 0; i < 2; i++)
      {
	// tree for MC
	TFile *fin = new TFile(Form("/eos/user/m/maaraujo/JpsiRun2/MC/%s", name[i_y][i].c_str()));
	TTree *tree = (TTree*)fin->Get("jpsitree");
    
	// creating desired vars and setting branch address
	TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
	float lt;
	
	tree->SetBranchAddress("muP_p4", &muP_p4);
	tree->SetBranchAddress("muM_p4", &muM_p4);
	tree->SetBranchAddress("mumu_p4", &mumu_p4);
	tree->SetBranchAddress("ctpv", &lt);
	
	// aux vars for reading tree
	int nEvt = tree->GetEntries();
	int perc = nEvt / 100;
      
	// preparing histograms to be filled
	TH1D *h_muPpT  = new TH1D(Form("h%d_muPpT_%s", i_y+7, lbl[i].c_str()),  "MC muon pT",      100, 0, 50);
	TH1D *h_muNpT  = new TH1D(Form("h%d_muNpT_%s", i_y+7, lbl[i].c_str()),  "muonN pT",        100, 0, 50);
	TH1D *h_muPEta = new TH1D(Form("h%d_muPEta_%s", i_y+7, lbl[i].c_str()), "MC muon eta",     120, -1.5, 1.5);
	TH1D *h_muNEta = new TH1D(Form("h%d_muNEta_%s", i_y+7, lbl[i].c_str()), "muonN eta",       120, -1.5, 1.5);
	TH1D *h_JMass  = new TH1D(Form("h%d_JMass_%s", i_y+7, lbl[i].c_str()),  "MC dimuon mass",  160, 2.9, 3.3);
	TH1D *h_JPt    = new TH1D(Form("h%d_JPt_%s", i_y+7, lbl[i].c_str()),    "MC dimuon pT",    120, 20, 50);
	TH1D *h_Jy     = new TH1D(Form("h%d_Jy_%s", i_y+7, lbl[i].c_str()),     "MC dimuon y",     120, -1.5, 1.5);
	TH1D *h_Jlt    = new TH1D(Form("h%d_Jlt_%s", i_y+7, lbl[i].c_str()),    "MC lifetime",     160, -0.02, 0.02);
	
	// reading tree, filling histograms
	for(int i = 0; i < nEvt; i++) {
	  tree->GetEntry(i);
	  if(mumu_p4->Pt() > 25 && mumu_p4->Pt() < 46)
	    {
	      h_muPpT->Fill(muP_p4->Pt());
	      h_muNpT->Fill(muM_p4->Pt());
	      h_muPEta->Fill(muP_p4->Eta());
	      h_muNEta->Fill(muM_p4->Eta());
	      h_JMass->Fill(mumu_p4->M());
	      h_JPt->Fill(mumu_p4->Pt());
	      h_Jy->Fill(mumu_p4->Rapidity());
	      h_Jlt->Fill(lt);
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
	h_Jlt->Write();
	outfile->Close();

	fin->Close();
      } 
  }
}
