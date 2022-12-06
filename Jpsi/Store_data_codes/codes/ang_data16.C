// code that runs over data, applies all cuts and saves a root file with Jpsi pT, y, mass and lts + costheta_HX
// current cuts: based on trigger and distributions; no mass/lts cut

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


void ang_data16()
{
  Double_t mmPt, ct, ctErr, vProb, th, phi, lt, lterr, mass, rap, mPPt, mMPt, mPEta, mMEta, dR;
  double dPhi, dEta, dpT;
  double *ang;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
  UInt_t trigger;

  // beam and target vectors (always the same)
  double pbeam = sqrts/2.;
  double Ebeam = sqrt(pbeam*pbeam + Mprot*Mprot);
  
  TLorentzVector *beam = new TLorentzVector();
  TLorentzVector *targ = new TLorentzVector();
  beam->SetPxPyPzE( 0., 0., pbeam, Ebeam);
  targ->SetPxPyPzE( 0., 0., -pbeam, Ebeam);

  // 2016 tree
  TFile *fin6 = new TFile("/eos/user/m/maaraujo/JpsiRun2/Data/filtered-all-psi-UL16.root");
  TTree *tree6 = (TTree*)fin6->Get("jpsitree");
  
  // setting branch address  
  tree6->SetBranchAddress("muP_p4", &muP_p4);
  tree6->SetBranchAddress("muM_p4", &muM_p4);
  tree6->SetBranchAddress("mumu_p4", &mumu_p4);
  tree6->SetBranchAddress("vProb", &vProb);
  tree6->SetBranchAddress("trigger", &trigger);
  tree6->SetBranchAddress("ctpv", &ct);
  tree6->SetBranchAddress("ctpv_error", &ctErr);
  
  int dEvt = tree6->GetEntries();
  int perc = dEvt / 100;

  TFile *fout6 = new TFile("data16_cos.root", "recreate");
  TTree *newtree6 = new TTree("data_cos", "");
  
  newtree6->Branch("theta", &th);
  newtree6->Branch("phi", &phi);
  newtree6->Branch("dimPt", &mmPt);
  newtree6->Branch("lt", &lt);
  newtree6->Branch("lterr", &lterr);
  newtree6->Branch("Mass", &mass);
  newtree6->Branch("Rap", &rap);
  newtree6->Branch("muonPPt", &mPPt);
  newtree6->Branch("muonPEta", &mPEta);
  newtree6->Branch("muonMPt", &mMPt);
  newtree6->Branch("muonMEta", &mMEta);
  newtree6->Branch("DeltaR", &dR);
  
  for(int i = 0; i < dEvt; i++) {
    tree6->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	vProb > 0.01 &&
	(trigger&2) == 2)
      {
	mPPt = muP_p4->Pt();
	mMPt = muM_p4->Pt();
	mPEta = muP_p4->Eta();
	mMEta = muM_p4->Eta();

	mmPt = mumu_p4->Pt();
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
	
	newtree6->Fill();
	
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2016 data" << endl; 
  }

  fout6->Write();
  fout6->Close();
  fin6->Close();
}
