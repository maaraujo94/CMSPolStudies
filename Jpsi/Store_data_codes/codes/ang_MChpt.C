// code that runs over MC (high pT), applies all cuts and saves a root file with Jpsi pT, y and mass + costheta_HX
// current cuts: based on trigger and distributions; no mass/ cut

const double gPI = TMath::Pi();
const double Mprot = 0.9382720;
const double sqrts = 13000.;

// macro to calculate angles
double *cos_B(TLorentzVector *B, TLorentzVector *psi, TLorentzVector *beam, TLorentzVector *targ)
{
  double y = B->Rapidity();
  static double ang[2];

  // B = B 4-vector in the laboratory

  TVector3 cm_to_B = -B->BoostVector();
  TVector3 B_to_cm = B->BoostVector();
  
  // calculate reference directions in the B rest frame

  TLorentzVector targ_B = *targ;  // beam, targ: 4-vectors of beams having positive (beam) and negative (targ) direction in the laboratory (use same convention as for the jpsi decay!)
  targ_B.Boost(cm_to_B);         // target (= negative beam) in the B rest frame
  TLorentzVector beam_B = *beam;
  beam_B.Boost(cm_to_B);         // (positive) beam in the B rest frame

  TVector3 beam_direction_B     = beam_B.Vect().Unit();
  TVector3 targ_direction_B     = targ_B.Vect().Unit();
  TVector3 B_direction          = B->Vect().Unit();

  TVector3 Yaxis = ( beam_direction_B.Cross( targ_direction_B ) ).Unit(); // use same convention as for the jpsi (beam x targ or targ x beam or particle_direction x beam? the sign depends on the convention)
  // symmetrize definition between rapidity>0 and rapidity <0:
  if ( y < 0 ) Yaxis = -Yaxis;

  // boost psi from the B rest frame into the proton-proton CM frame:

  TLorentzVector psi_B = *psi;   // psi = psi 4-vector in the lab
  psi_B.Boost(cm_to_B);


  TRotation rotation;

  TVector3 newYaxis = Yaxis;
  TVector3 newZaxis = B_direction;
  TVector3 newXaxis = newYaxis.Cross(newZaxis);

  rotation.SetToIdentity();
  rotation.RotateAxes(newXaxis,newYaxis,newZaxis);
  rotation.Invert(); // transforms coordinates from the "oldXYZ" frame to the "newXYZ" frame

  TVector3 psi_B_rotated = psi_B.Vect();

  psi_B_rotated.Transform(rotation);

  double TH = psi_B_rotated.Theta();

  double PHI = psi_B_rotated.Phi() * 180. / gPI;

  ang[0] = TH;
  ang[1] = PHI;

  return ang;
}

void ang_MChpt()
{
  Double_t mmPt, th, phi, lt, lterr, mass, rap, mPPt, mMPt, mPEta, mMEta, dR, mmPhi;
  double dPhi, dEta, dpT;
  double *ang;
  Float_t ct, ctErr;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;

  // beam and target vectors (always the same)
  double pbeam = sqrts/2.;
  double Ebeam = sqrt(pbeam*pbeam + Mprot*Mprot);
  
  TLorentzVector *beam = new TLorentzVector();
  TLorentzVector *targ = new TLorentzVector();
  beam->SetPxPyPzE( 0., 0., pbeam, Ebeam);
  targ->SetPxPyPzE( 0., 0., -pbeam, Ebeam);

  // 2017 tree
  TFile *fin7 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-mc-LOCAL17-highpt.root");
  TTree *tree7 = (TTree*)fin7->Get("jpsitree");
  
  tree7->SetBranchAddress("muP_p4", &muP_p4);
  tree7->SetBranchAddress("muM_p4", &muM_p4);
  tree7->SetBranchAddress("mumu_p4", &mumu_p4);
  tree7->SetBranchAddress("ctpv", &ct);
  tree7->SetBranchAddress("ctpv_error", &ctErr);

  int mEvt = tree7->GetEntries();
  int perc = mEvt / 100;

  TFile *fout7 = new TFile("MCh17_cos.root", "recreate");
  TTree *newtree7 = new TTree("MC_cos", "");

  newtree7->Branch("theta", &th);
  newtree7->Branch("phi", &phi);
  newtree7->Branch("dimPt", &mmPt);
  newtree7->Branch("lt", &lt);
  newtree7->Branch("lterr", &lterr);
  newtree7->Branch("Mass", &mass);
  newtree7->Branch("Rap", &rap);
  newtree7->Branch("dimPhi", &mmPhi);
  newtree7->Branch("muonPPt", &mPPt);
  newtree7->Branch("muonPEta", &mPEta);
  newtree7->Branch("muonMPt", &mMPt);
  newtree7->Branch("muonMEta", &mMEta);
  newtree7->Branch("DeltaR", &dR);
 
  for(int i = 0; i < mEvt; i++) {
    tree7->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4)
      {
	mPPt = muP_p4->Pt();
	mMPt = muM_p4->Pt();
	mPEta = muP_p4->Eta();
	mMEta = muM_p4->Eta();
	
	mmPt = mumu_p4->Pt();
	mmPhi = mumu_p4->Phi();
	rap = mumu_p4->Rapidity();
	mass = mumu_p4->M();
	lt = ct;
	lterr = ctErr;
	
	dPhi = muP_p4->Phi() - muM_p4->Phi();
	dEta = muP_p4->Eta() - muM_p4->Eta();
	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;
	
	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;

	ang = cos_B(mumu_p4, muP_p4, beam, targ);
	th = ang[0];
	phi = ang[1];

	newtree7->Fill();
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2017 MC" << endl; 
  }
  
  fout7->Write();
  fout7->Close();
  fin7->Close();

  // 2018 tree
  TFile *fin8 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-mc-LOCAL18-highpt.root");
  TTree *tree8 = (TTree*)fin8->Get("jpsitree");
  
  tree8->SetBranchAddress("muP_p4", &muP_p4);
  tree8->SetBranchAddress("muM_p4", &muM_p4);
  tree8->SetBranchAddress("mumu_p4", &mumu_p4);
  tree8->SetBranchAddress("ctpv", &ct);
  tree8->SetBranchAddress("ctpv_error", &ctErr);

  mEvt = tree8->GetEntries();
  perc = mEvt / 100;

  TFile *fout8 = new TFile("MCh18_cos.root", "recreate");
  TTree *newtree8 = new TTree("MC_cos", "");

  newtree8->Branch("theta", &th);
  newtree8->Branch("phi", &phi);
  newtree8->Branch("dimPt", &mmPt);
  newtree8->Branch("lt", &lt);
  newtree8->Branch("lterr", &lterr);
  newtree8->Branch("Mass", &mass);
  newtree8->Branch("Rap", &rap);
  newtree8->Branch("dimPhi", &mmPhi);
  newtree8->Branch("muonPPt", &mPPt);
  newtree8->Branch("muonPEta", &mPEta);
  newtree8->Branch("muonMPt", &mMPt);
  newtree8->Branch("muonMEta", &mMEta);
  newtree8->Branch("DeltaR", &dR);

  for(int i = 0; i < mEvt; i++) {
    tree8->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4  )
      {
	mPPt = muP_p4->Pt();
	mMPt = muM_p4->Pt();
	mPEta = muP_p4->Eta();
	mMEta = muM_p4->Eta();

	mmPt = mumu_p4->Pt();
	mmPhi = mumu_p4->Phi();
	rap = mumu_p4->Rapidity();
	mass = mumu_p4->M();
	lt = ct;
	lterr = ctErr;
	
	dPhi = muP_p4->Phi() - muM_p4->Phi();
	dEta = muP_p4->Eta() - muM_p4->Eta();
	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;
	
	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;

	ang = cos_B(mumu_p4, muP_p4, beam, targ);
	th = ang[0];
	phi = ang[1];

	newtree8->Fill();
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2018 MC" << endl; 
  }
  
  fout8->Write();
  fout8->Close();
  fin8->Close();
 
}
