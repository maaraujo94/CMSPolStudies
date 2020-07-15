// code that runs over data, applies all cuts but mass and saves a root file with pT, y, mass and costheta_HX - for mass window studies
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

void jpsi2012_massdata()
{
  TChain *dataJ = new TChain("jpsi_tuple");

  dataJ->Add("/eos/user/m/maaraujo/Jpsi2012/Data/r2012B_MuPk_jpsi_v8_flat_skim.root");
  dataJ->Add("/eos/user/m/maaraujo/Jpsi2012/Data/r2012C_MuPk_jpsi_v8_1_flat_skim.root");
  dataJ->Add("/eos/user/m/maaraujo/Jpsi2012/Data/r2012D_MuPk_jpsi_v8_1_flat_skim.root");

  Double_t muPPt, muNPt, muPEta, muNEta, muPPhi, muPMass, JpsiMass, JpsiPt, JpsiRap, JpsiEta, JpsiPhi, Jpsict, JpsictErr, vProb;
  Long64_t trigger;

  dataJ->SetBranchAddress("muPPt", &muPPt);
  dataJ->SetBranchAddress("muNPt", &muNPt);
  dataJ->SetBranchAddress("muPEta", &muPEta);
  dataJ->SetBranchAddress("muNEta", &muNEta);
  dataJ->SetBranchAddress("muPPhi", &muPPhi);
  dataJ->SetBranchAddress("muPMass", &muPMass);
  dataJ->SetBranchAddress("JpsiMass", &JpsiMass);
  dataJ->SetBranchAddress("JpsiPt", &JpsiPt);
  dataJ->SetBranchAddress("JpsiRap", &JpsiRap);
  dataJ->SetBranchAddress("JpsiEta", &JpsiEta);
  dataJ->SetBranchAddress("JpsiPhi", &JpsiPhi);
  dataJ->SetBranchAddress("Jpsict", &Jpsict);
  dataJ->SetBranchAddress("JpsictErr", &JpsictErr);
  dataJ->SetBranchAddress("JpsiVprob", &vProb);
  dataJ->SetBranchAddress("trigger", &trigger);

  int dEvt = dataJ->GetEntries();
  int perc = dEvt / 100;

  TFile *outfile = new TFile("data_mass_cos.root", "recreate");
  TTree *newtree = new TTree("data_cos", "");

  double cosa;
  TBranch *cos_tree = newtree->Branch("costh", &cosa);
  TBranch *pT_tree = newtree->Branch("JpsiPt", &JpsiPt);
  TBranch *y_tree = newtree->Branch("JpsiRap", &JpsiRap);
  TBranch *mass_tree = newtree->Branch("JpsiMass", &JpsiMass);
  
  for(int i = 0; i < dEvt; i++) {
    dataJ->GetEntry(i);
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
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with data" << endl; 
  }

  outfile->Write();
  outfile->Close();
}
