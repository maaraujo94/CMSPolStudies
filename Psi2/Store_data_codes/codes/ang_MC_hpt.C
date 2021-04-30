// code that runs over MC, applies all cuts and saves a root file with Jpsi pT, y and mass + costheta_HX
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

void ang_MC_hpt()
{
  TChain *tree = new TChain("psi2stree");

  tree->Add("/eos/user/m/maaraujo/Psi2SRun2/filtered-all-psi2s-mc-LOCAL18-highpt.root");

  Double_t mmPt, th, phi, mass, rap, lts;
  double *ang;
  Float_t ct, ctErr;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
  
  tree->SetBranchAddress("muP_p4", &muP_p4);
  tree->SetBranchAddress("muM_p4", &muM_p4);
  tree->SetBranchAddress("mumu_p4", &mumu_p4);
  tree->SetBranchAddress("ctpv", &ct);
  tree->SetBranchAddress("ctpv_error", &ctErr);
 
  int mEvt = tree->GetEntries();
  int perc = mEvt / 100;

  TFile *outfile = new TFile("MC18_hpt_cos.root", "recreate");
  TTree *newtree = new TTree("MC_cos", "");

  newtree->Branch("theta", &th);
  newtree->Branch("phi", &phi);
  newtree->Branch("dimPt", &mmPt);
  newtree->Branch("lts", &lts);
  newtree->Branch("Mass", &mass);
  newtree->Branch("Rap", &rap);

  // beam and target vectors (always the same)
  double pbeam = sqrts/2.;
  double Ebeam = sqrt(pbeam*pbeam + Mprot*Mprot);
  
  TLorentzVector *beam = new TLorentzVector();
  TLorentzVector *targ = new TLorentzVector();
  beam->SetPxPyPzE( 0., 0., pbeam, Ebeam);
  targ->SetPxPyPzE( 0., 0., -pbeam, Ebeam);

  for(int i = 0; i < mEvt; i++) {
    tree->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	abs(mumu_p4->Rapidity()) < 1.2  )
      {
	mmPt = mumu_p4->Pt();
	rap = abs(mumu_p4->Rapidity());
	mass = mumu_p4->M();
	lts = ct/ctErr;
	
	ang = cos_B(mumu_p4, muP_p4, beam, targ);
	th = ang[0];
	phi = ang[1];

	newtree->Fill();
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with MC" << endl; 
  }

  outfile->Write();
  outfile->Close();
  
}
