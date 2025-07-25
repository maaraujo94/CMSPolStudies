// code that runs over data, applies all cuts and saves a root file with Jpsi pT, y, mass and lts + costheta_HX
// current cuts: based on trigger and distributions; no mass/lts cut

const double gPI = TMath::Pi();
const double Mprot = 0.9382720;
const double sqrts = 13000.;

// macro to calculate angles
double *cos_psi(TLorentzVector *psi, TLorentzVector *mu, TLorentzVector *beam, TLorentzVector *targ)
{
  double y = psi->Rapidity();
  static double ang[2];

  // psi = psi 4-vector in the laboratory

  TVector3 cm_to_psi = -psi->BoostVector();
  TVector3 psi_to_cm = psi->BoostVector();
  
  // calculate reference directions in the psi rest frame

  TLorentzVector targ_psi = *targ;  // beam, targ: 4-vectors of beams having positive (beam) and negative (targ) direction in the laboratory (use same convention as for the jpsi decay!)
  targ_psi.Boost(cm_to_psi);         // target (= negative beam) in the psi rest frame
  TLorentzVector beam_psi = *beam;
  beam_psi.Boost(cm_to_psi);         // (positive) beam in the psi rest frame

  TVector3 beam_direction_psi     = beam_psi.Vect().Unit();
  TVector3 targ_direction_psi     = targ_psi.Vect().Unit();
  TVector3 psi_direction          = psi->Vect().Unit();

  TVector3 Yaxis = ( beam_direction_psi.Cross( targ_direction_psi ) ).Unit(); // use same convention as for the jpsi (beam x targ or targ x beam or particle_direction x beam? the sign depends on the convention)
  // symmetrize definition between rapidity>0 and rapidity <0:
  if ( y < 0 ) Yaxis = -Yaxis;

  // boost mu from the psi rest frame into the proton-proton CM frame:

  TLorentzVector mu_psi = *mu;   // mu = mu 4-vector in the lab
  mu_psi.Boost(cm_to_psi);


  TRotation rotation;

  TVector3 newYaxis = Yaxis;
  TVector3 newZaxis = psi_direction;
  TVector3 newXaxis = newYaxis.Cross(newZaxis);

  rotation.SetToIdentity();
  rotation.RotateAxes(newXaxis,newYaxis,newZaxis);
  rotation.Invert(); // transforms coordinates from the "oldXYZ" frame to the "newXYZ" frame

  TVector3 mu_psi_rotated = mu_psi.Vect();

  mu_psi_rotated.Transform(rotation);

  double TH = mu_psi_rotated.Theta();

  double PHI = mu_psi_rotated.Phi() * 180. / gPI;

  ang[0] = TH;
  ang[1] = PHI;

  return ang;
}


void ang_data()
{
  Double_t mmPt, ct, ctErr, vProb, th, phi, lt, lterr, mass, rap, mPPt, mMPt, mPEta, mMEta, dR, mmPhi;
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

  // 2017 tree
  TFile *fin7 = new TFile("/eos/user/m/maaraujo/JpsiRun2/Data/filtered-all-psi-UL17.root");
  TTree *tree7 = (TTree*)fin7->Get("jpsitree");
  
  // setting branch address  
  tree7->SetBranchAddress("muP_p4", &muP_p4);
  tree7->SetBranchAddress("muM_p4", &muM_p4);
  tree7->SetBranchAddress("mumu_p4", &mumu_p4);
  tree7->SetBranchAddress("vProb", &vProb);
  tree7->SetBranchAddress("trigger", &trigger);
  tree7->SetBranchAddress("ctpv", &ct);
  tree7->SetBranchAddress("ctpv_error", &ctErr);
  
  int dEvt = tree7->GetEntries();
  int perc = dEvt / 100;

  TFile *fout7 = new TFile("data17_cos.root", "recreate");
  TTree *newtree7 = new TTree("data_cos", "");
  
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
  
  for(int i = 0; i < dEvt; i++) {
    tree7->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	vProb > 0.01 &&
	(trigger&16) == 16)
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
	
	ang = cos_psi(mumu_p4, muP_p4, beam, targ);
	th = ang[0];
	phi = ang[1];
	
	newtree7->Fill();
	
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2017 data" << endl; 
  }

  fout7->Write();
  fout7->Close();
  fin7->Close();

  // 2018 tree
  TFile *fin8 = new TFile("/eos/user/m/maaraujo/JpsiRun2/Data/filtered-all-psi-UL18.root");
  TTree *tree8 = (TTree*)fin8->Get("jpsitree");
  
  // setting branch address  
  tree8->SetBranchAddress("muP_p4", &muP_p4);
  tree8->SetBranchAddress("muM_p4", &muM_p4);
  tree8->SetBranchAddress("mumu_p4", &mumu_p4);
  tree8->SetBranchAddress("vProb", &vProb);
  tree8->SetBranchAddress("trigger", &trigger);
  tree8->SetBranchAddress("ctpv", &ct);
  tree8->SetBranchAddress("ctpv_error", &ctErr);
  
  dEvt = tree8->GetEntries();
  perc = dEvt / 100;

  TFile *fout8 = new TFile("data18_cos.root", "recreate");
  TTree *newtree8 = new TTree("data_cos", "");
  
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

  for(int i = 0; i < dEvt; i++) {
    tree8->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	vProb > 0.01 &&
	(trigger&16) == 16)
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

	ang = cos_psi(mumu_p4, muP_p4, beam, targ);
	th = ang[0];
	phi = ang[1];

	dPhi = muP_p4->Phi() - muM_p4->Phi();
	dEta = muP_p4->Eta() - muM_p4->Eta();
	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;
	
	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;

	newtree8->Fill();
	
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2018 data" << endl; 
  }
  
  fout8->Write();
  fout8->Close();
  fin8->Close();

}
