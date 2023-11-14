// code that runs over MC (trig and reco), applies cuts, and gets the costh dists

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

void ang_MC_reco()
{
  Double_t mmPt, th, lt, mass, rap, dR;
  double dPhi, dEta, dpT;
  double *ang;
  Float_t ct;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
  int mEvt, perc;
  
  // beam and target vectors (always the same)
  double pbeam = sqrts/2.;
  double Ebeam = sqrt(pbeam*pbeam + Mprot*Mprot);
  
  TLorentzVector *beam = new TLorentzVector();
  TLorentzVector *targ = new TLorentzVector();
  beam->SetPxPyPzE( 0., 0., pbeam, Ebeam);
  targ->SetPxPyPzE( 0., 0., -pbeam, Ebeam);
  
  // reco trees
  TFile *finNT1 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-notrigger_mc-LOCAL18.root");
  TTree *treeNT1 = (TTree*)finNT1->Get("jpsitree");
  TFile *finNT2 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-notrigger_mc-LOCAL18-highpt.root");
  TTree *treeNT2 = (TTree*)finNT2->Get("jpsitree");
  TFile *finNT3 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-notrigger_mc-LOCAL18-veryhighpt.root");
  TTree *treeNT3 = (TTree*)finNT3->Get("jpsitree");

  // low-pT
  mEvt = treeNT1->GetEntries();
  perc = mEvt / 100;

  treeNT1->SetBranchAddress("muP_p4", &muP_p4);
  treeNT1->SetBranchAddress("muM_p4", &muM_p4);
  treeNT1->SetBranchAddress("mumu_p4", &mumu_p4);
  treeNT1->SetBranchAddress("ctpv", &ct);

  TFile *fout1 = new TFile("MC18r_cos.root", "recreate");
  TTree *newtree1 = new TTree("MC_cos", "");

  newtree1->Branch("theta", &th);
  newtree1->Branch("dimPt", &mmPt);
  newtree1->Branch("lt", &lt);
  newtree1->Branch("Mass", &mass);
  newtree1->Branch("Rap", &rap);
  newtree1->Branch("DeltaR", &dR);

  for(int i = 0; i < mEvt; i++) {
    treeNT1->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4)
      {
	mmPt = mumu_p4->Pt();
	rap = mumu_p4->Rapidity();
	mass = mumu_p4->M();
	lt = ct;

	dPhi = muP_p4->Phi() - muM_p4->Phi();
	dEta = muP_p4->Eta() - muM_p4->Eta();
	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;

	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;

	ang = cos_B(mumu_p4, muP_p4, beam, targ);
	th = ang[0];

	newtree1->Fill();
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with low-pT reco MC" << endl; 
  }
  fout1->Write();
  fout1->Close();
  finNT1->Close();
  
  // mid-pT
  mEvt = treeNT2->GetEntries();
  perc = mEvt / 100;

  treeNT2->SetBranchAddress("muP_p4", &muP_p4);
  treeNT2->SetBranchAddress("muM_p4", &muM_p4);
  treeNT2->SetBranchAddress("mumu_p4", &mumu_p4);
  treeNT2->SetBranchAddress("ctpv", &ct);

  TFile *fout2 = new TFile("MCh18r_cos.root", "recreate");
  TTree *newtree2 = new TTree("MC_cos", "");

  newtree2->Branch("theta", &th);
  newtree2->Branch("dimPt", &mmPt);
  newtree2->Branch("lt", &lt);
  newtree2->Branch("Mass", &mass);
  newtree2->Branch("Rap", &rap);
  newtree2->Branch("DeltaR", &dR);

  for(int i = 0; i < mEvt; i++) {
    treeNT2->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4)
      {
	mmPt = mumu_p4->Pt();
	rap = mumu_p4->Rapidity();
	mass = mumu_p4->M();
	lt = ct;

	dPhi = muP_p4->Phi() - muM_p4->Phi();
	dEta = muP_p4->Eta() - muM_p4->Eta();
	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;

	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;

	ang = cos_B(mumu_p4, muP_p4, beam, targ);
	th = ang[0];

	newtree2->Fill();
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with mid-pT reco MC" << endl; 
  }
  fout2->Write();
  fout2->Close();
  finNT2->Close();

  // high-pT
  mEvt = treeNT3->GetEntries();
  perc = mEvt / 100;

  treeNT3->SetBranchAddress("muP_p4", &muP_p4);
  treeNT3->SetBranchAddress("muM_p4", &muM_p4);
  treeNT3->SetBranchAddress("mumu_p4", &mumu_p4);
  treeNT3->SetBranchAddress("ctpv", &ct);

  TFile *fout3 = new TFile("MCvh18r_cos.root", "recreate");
  TTree *newtree3 = new TTree("MC_cos", "");

  newtree3->Branch("theta", &th);
  newtree3->Branch("dimPt", &mmPt);
  newtree3->Branch("lt", &lt);
  newtree3->Branch("Mass", &mass);
  newtree3->Branch("Rap", &rap);
  newtree3->Branch("DeltaR", &dR);

  for(int i = 0; i < mEvt; i++) {
    treeNT3->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4)
      {
	mmPt = mumu_p4->Pt();
	rap = mumu_p4->Rapidity();
	mass = mumu_p4->M();
	lt = ct;

	dPhi = muP_p4->Phi() - muM_p4->Phi();
	dEta = muP_p4->Eta() - muM_p4->Eta();
	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;

	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;

	ang = cos_B(mumu_p4, muP_p4, beam, targ);
	th = ang[0];

	newtree3->Fill();
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with high-pT reco MC" << endl; 
  }
  fout3->Write();
  fout3->Close();
  finNT3->Close();
}
