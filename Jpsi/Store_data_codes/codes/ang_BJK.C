// code to run on lxplus to get angular dists in cut sample
// runs on B->J/psi K sample

const double gPI = TMath::Pi();
const double Mprot = 0.9382720;
const double sqrts = 13000.;

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
  //double costheta_HX = cos(theta_HX);

  return theta_HX;
}

void ang_BJK()
{
  // tree for data
  TChain *tree = new TChain("treeS");
  
  tree->Add("/eos/user/m/maaraujo/BtoPsiK/filtered-all-bupsik-psi25-DT18.root");

  // creating desired vars and setting branch address
  Double_t muPPt, muNPt, muPEta, muNEta, JpsiMass, JpsiPt, JpsiRap, Jpsict, JpsictErr, massB, BPt, lts;
  TLorentzVector *mumu_p4 = 0, *muP_p4 = 0, *muM_p4 = 0, *bu_p4 = 0;

  tree->SetBranchAddress("muP_p4", &muP_p4);
  tree->SetBranchAddress("muM_p4", &muM_p4);
  tree->SetBranchAddress("mumu_p4", &mumu_p4);
  tree->SetBranchAddress("ctpv", &Jpsict);
  tree->SetBranchAddress("ctpv_error", &JpsictErr);
  tree->SetBranchAddress("bu_p4", &bu_p4);

  // aux vars for reading tree
  int nEvt = tree->GetEntries();
  int perc = nEvt / 100;

  TFile *outfile = new TFile("btopk_cos.root", "recreate");
  TTree *newtree = new TTree("data_cos", "");

  double *ang_B, *ang_P;
  double B_th, B_phi, J_th, J_phi, J_th_comp;
  
  TBranch *TH_tree = newtree->Branch("TH", &B_th);
  TBranch *th_tree = newtree->Branch("th", &J_th);
  TBranch *PHI_tree = newtree->Branch("PHI", &B_phi);
  TBranch *phi_tree = newtree->Branch("phi", &J_phi);
  TBranch *lts_tree = newtree->Branch("lts", &lts);
  TBranch *JpT_tree = newtree->Branch("JpsiPt", &JpsiPt);
  TBranch *BpT_tree = newtree->Branch("BPt", &BPt);
  TBranch *BM_tree = newtree->Branch("BMass", &massB);
  TBranch *th_comp_tree = newtree->Branch("th_comp", &J_th_comp);

  // beam and target vectors (always the same)
  double pbeam = sqrts/2.;
  double Ebeam = sqrt(pbeam*pbeam + Mprot*Mprot);
  
  TLorentzVector *beam = new TLorentzVector();
  TLorentzVector *targ = new TLorentzVector();
  beam->SetPxPyPzE( 0., 0., pbeam, Ebeam);
  targ->SetPxPyPzE( 0., 0., -pbeam, Ebeam);
  
  // reading tree, filling histograms
  for(int i = 0; i < nEvt; i++) {
    tree->GetEntry(i);

    if(muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6
       //&& abs(muP_p4->Eta()) < 1.6 && abs(muM_p4->Eta()) < 1.6
       && abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4
       && mumu_p4->M() < 3.2 && mumu_p4->M() > 3.
       && abs(mumu_p4->Rapidity()) < 1.2) {
       //       && bu_p4->M() < 5.32 && bu_p4->M() > 5.24) {
    
      lts = abs(Jpsict/JpsictErr);
      JpsiPt = mumu_p4->Pt();
      BPt = bu_p4->Pt();
      massB = bu_p4->M();

      ang_B = cos_B(bu_p4, mumu_p4, beam, targ);
      B_th = ang_B[0];
      B_phi = ang_B[1];
      ang_P = cos_B(mumu_p4, muP_p4, beam, targ);
      J_th = ang_P[0];
      J_phi = ang_P[1];

      J_th_comp = costh(mumu_p4, muP_p4);
      
      newtree->Fill();
    }
    
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with data" << endl; 
  }

  outfile->Write();
  outfile->Close();
    
 
}
