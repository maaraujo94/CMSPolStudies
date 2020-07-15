// code that runs over MC, applies all cuts but mass and saves a root file with pT, y, mass and costheta_HX - for mass window studies
// current cuts: same as BPH-13-001 AN, no mass cut

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

void jpsi2012_massMC()
{
  TChain *mcJ = new TChain("jpsi_tuple");

  mcJ->Add("/eos/user/m/maaraujo/Jpsi2012/MC/onia2MuMu_tree_validation_flat_skim.root");

  Double_t muPPt, muNPt, muPEta, muNEta, muPPhi, muPMass, JpsiMass, JpsiPt, JpsiRap, JpsiEta, JpsiPhi, Jpsict, JpsictErr, vProb;
  Long64_t trigger;

  mcJ->SetBranchAddress("muPPt", &muPPt);
  mcJ->SetBranchAddress("muNPt", &muNPt);
  mcJ->SetBranchAddress("muPEta", &muPEta);
  mcJ->SetBranchAddress("muNEta", &muNEta);
  mcJ->SetBranchAddress("muPPhi", &muPPhi);
  mcJ->SetBranchAddress("muPMass", &muPMass);
  mcJ->SetBranchAddress("JpsiMass", &JpsiMass);
  mcJ->SetBranchAddress("JpsiPt", &JpsiPt);
  mcJ->SetBranchAddress("JpsiRap", &JpsiRap);
  mcJ->SetBranchAddress("JpsiEta", &JpsiEta);
  mcJ->SetBranchAddress("JpsiPhi", &JpsiPhi);
  mcJ->SetBranchAddress("Jpsict", &Jpsict);
  mcJ->SetBranchAddress("JpsictErr", &JpsictErr);
  mcJ->SetBranchAddress("JpsiVprob", &vProb);
  mcJ->SetBranchAddress("trigger", &trigger);

  int mEvt = mcJ->GetEntries();
  int perc = mEvt / 100;

  TFile *outfile = new TFile("MC_mass_cos.root", "recreate");
  TTree *newtree = new TTree("MC_cos", "");

  double cosa;
  TBranch *cos_tree = newtree->Branch("costh", &cosa);
  TBranch *pT_tree = newtree->Branch("JpsiPt", &JpsiPt);
  TBranch *y_tree = newtree->Branch("JpsiRap", &JpsiRap);
  TBranch *mass_tree = newtree->Branch("JpsiMass", &JpsiMass);
  
  for(int i = 0; i < mEvt; i++) {
    mcJ->GetEntry(i);
    if( muPPt > 5.6 && muNPt > 5.6 && 
	abs(muPEta) < 1.6 && abs(muNEta) < 1.6 &&
	vProb > 0.01 &&
	JpsiPt > 12 && JpsiPt < 70 &&
	abs(JpsiRap) < 1.2 &&
	abs(Jpsict/JpsictErr) < 2.5 &&
	trigger == 1)
      {
	TLorentzVector *p4_jpsi = new TLorentzVector;
	p4_jpsi->SetPtEtaPhiM(JpsiPt, JpsiEta, JpsiPhi, JpsiMass);
	TLorentzVector *p4_muP = new TLorentzVector;
	p4_muP->SetPtEtaPhiM(muPPt, muPEta, muPPhi, muPMass);
	      
	cosa = costh(p4_jpsi, p4_muP);

	newtree->Fill();
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with MC" << endl; 
  }

  outfile->Write();
  outfile->Close();
  
}
