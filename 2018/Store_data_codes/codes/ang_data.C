// code that runs over data, applies all cuts and saves a root file with Jpsi pT, y, mass and lts + costheta_HX
// current cuts: based on trigger and distributions; no mass/lts cut

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

void ang_data()
{
  // open tree and read data
  TChain *dataJ = new TChain("jpsitree");

  dataJ->Add("/eos/user/m/maaraujo/JpsiRun2/Data/filtered-17-18-psi-19aug20.root");

  Double_t JpsiPt, Jpsict, JpsictErr, vProb, cosa, lts, mass, rap;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
  UInt_t trigger, run;
  
  dataJ->SetBranchAddress("muP_p4", &muP_p4);
  dataJ->SetBranchAddress("muM_p4", &muM_p4);
  dataJ->SetBranchAddress("mumu_p4", &mumu_p4);
  dataJ->SetBranchAddress("vProb", &vProb);
  dataJ->SetBranchAddress("trigger", &trigger);
  dataJ->SetBranchAddress("run", &run);
  dataJ->SetBranchAddress("ctpv", &Jpsict);
  dataJ->SetBranchAddress("ctpv_error", &JpsictErr);
  
  int dEvt = dataJ->GetEntries();
  int perc = dEvt / 100;

  TFile *outfile = new TFile("data18_cos.root", "recreate");
  TTree *newtree = new TTree("data_cos", "");

  newtree->Branch("costh", &cosa);
  newtree->Branch("JpsiPt", &JpsiPt);
  newtree->Branch("lts", &lts);
  newtree->Branch("JpsiMass", &mass);
  newtree->Branch("JpsiRap", &rap);

  for(int i = 0; i < dEvt; i++) {
    dataJ->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	vProb > 0.01 &&
	abs(mumu_p4->Rapidity()) < 1.2 &&
	(trigger&16) == 16 &&
	run > 313000)
      {
	JpsiPt = mumu_p4->Pt();
	rap = abs(mumu_p4->Rapidity());
	mass = mumu_p4->M();
	cosa = costh(mumu_p4, muP_p4);
	lts = abs(Jpsict/JpsictErr);
	  
	newtree->Fill();
	
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with data" << endl; 
  }
  
  outfile->Write();
  outfile->Close();
  
}
