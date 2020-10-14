// code that runs over MC (high pT), applies all cuts and saves a root file with Jpsi pT, y and mass + costheta_HX
// current cuts: based on trigger and distributions; no mass/ cut

// macro to calculate cos(theta)
double costh(TLorentzVector *p4_parent_lab, TLorentzVector *p4_daughter_lab)
{
  // boost that must be applied to get to parent frame
  TVector3 boostToParent = -p4_parent_lab->BoostVector();

  // apply boost to daughter particle
  TLorentzVector *p4_daughter_in_parent_restframe = p4_daughter_lab;
  p4_daughter_in_parent_restframe->Boost( boostToParent );

  // calculate theta as angle between daughter and parent
  double theta_HX = p4_daughter_in_parent_restframe->Angle(p4_parent_lab->Vect());
  double costheta_HX = cos(theta_HX);

  return costheta_HX;
}

void ang_MChpt()
{
  TChain *mcJ = new TChain("jpsitree");

  mcJ->Add("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-mc-LOCAL18-highpt.root");

  Double_t cosa, JpsiPt, mass, rap;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
  
  mcJ->SetBranchAddress("muP_p4", &muP_p4);
  mcJ->SetBranchAddress("muM_p4", &muM_p4);
  mcJ->SetBranchAddress("mumu_p4", &mumu_p4);
  
  int mEvt = mcJ->GetEntries();
  int perc = mEvt / 100;

  TFile *outfile = new TFile("MC18_hpt_cos.root", "recreate");
  TTree *newtree = new TTree("MC_cos", "");

  TBranch *cos_tree = newtree->Branch("costh", &cosa);
  TBranch *pT_tree = newtree->Branch("JpsiPt", &JpsiPt);
  newtree->Branch("JpsiMass", &mass);
  newtree->Branch("JpsiRap", &rap);
 
  for(int i = 0; i < mEvt; i++) {
    mcJ->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	abs(mumu_p4->Rapidity()) < 1.2  )
      {
	cosa = costh(mumu_p4, muP_p4);
	JpsiPt = mumu_p4->Pt();
	rap = abs(mumu_p4->Rapidity());
	mass = mumu_p4->M();

	newtree->Fill();
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with MC" << endl; 
  }

  outfile->Write();
  outfile->Close();
  
}
