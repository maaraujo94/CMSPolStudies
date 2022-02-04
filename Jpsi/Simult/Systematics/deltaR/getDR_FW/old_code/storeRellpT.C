// code that runs over MC, applies all cuts and saves a root file with the appropriate angular histos

const double gPI = TMath::Pi();

void storeRellpT()
{
  int mEvt, perc;
  Double_t dEta, dPhi, dR, dpT;
  Float_t ct;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;

  // output histos
  TH2D* r_RpT = new TH2D("r_RpT", "reco p_{T} > 50", 75, 0, 0.3, 75, 0, 120);
  TH2D* t_RpT = new TH2D("t_RpT", "trig p_{T} > 50", 75, 0, 0.3, 75, 0, 120);
  
  // reco trees
  TFile *finNT2 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-notrigger_mc-LOCAL18-highpt.root");
  TTree *treeNT2 = (TTree*)finNT2->Get("jpsitree");
  TFile *finNT3 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-notrigger_mc-LOCAL18-veryhighpt.root");
  TTree *treeNT3 = (TTree*)finNT3->Get("jpsitree");
  
  // mid-pT
  mEvt = treeNT2->GetEntries();
  perc = mEvt / 100;

  treeNT2->SetBranchAddress("muP_p4", &muP_p4);
  treeNT2->SetBranchAddress("muM_p4", &muM_p4);
  treeNT2->SetBranchAddress("mumu_p4", &mumu_p4);
  treeNT2->SetBranchAddress("ctpv", &ct);

  for(int i = 0; i < mEvt; i++) {
    treeNT2->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	mumu_p4->Pt() > 50 && mumu_p4->Pt() < 66 &&
	mumu_p4->M() > 2.92 && mumu_p4->M() < 3.28 &&
	abs(mumu_p4->Rapidity()) < 1.2 &&
	ct > -0.005 && ct < 0.05)
      {
	dEta = muP_p4->Eta() - muM_p4->Eta();
	dPhi = muP_p4->Phi() - muM_p4->Phi();
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;

	dR = sqrt(dEta*dEta+1.4*1.4*dPhi*dPhi);
	dpT = abs(muP_p4->Pt()-muM_p4->Pt());

	if(dPhi < 0)
	  r_RpT->Fill(dR, dpT);
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with mid-pT reco MC" << endl; 
  }

  // high-pT
  mEvt = treeNT3->GetEntries();
  perc = mEvt / 100;

  treeNT3->SetBranchAddress("muP_p4", &muP_p4);
  treeNT3->SetBranchAddress("muM_p4", &muM_p4);
  treeNT3->SetBranchAddress("mumu_p4", &mumu_p4);
  treeNT3->SetBranchAddress("ctpv", &ct);

  for(int i = 0; i < mEvt; i++) {
    treeNT3->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	mumu_p4->Pt() > 66 && mumu_p4->Pt() < 120 &&
	mumu_p4->M() > 2.92 && mumu_p4->M() < 3.28 &&
	abs(mumu_p4->Rapidity()) < 1.2 &&
	ct > -0.005 && ct < 0.05)
      {
	dEta = muP_p4->Eta() - muM_p4->Eta();
	dPhi = muP_p4->Phi() - muM_p4->Phi();
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;

	dR = sqrt(dEta*dEta+1.4*1.4*dPhi*dPhi);
	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	
	if(dPhi < 0)
	  r_RpT->Fill(dR, dpT);
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with high-pT reco MC" << endl; 
  }

  finNT2->Close();
  finNT3->Close();

  // reco trees
  TFile *finR2 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-mc-LOCAL18-highpt.root");
  TTree *treeR2 = (TTree*)finR2->Get("jpsitree");
  TFile *finR3 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-mc-LOCAL18-veryhighpt.root");
  TTree *treeR3 = (TTree*)finR3->Get("jpsitree");
  
  // mid-pT
  mEvt = treeR2->GetEntries();
  perc = mEvt / 100;

  treeR2->SetBranchAddress("muP_p4", &muP_p4);
  treeR2->SetBranchAddress("muM_p4", &muM_p4);
  treeR2->SetBranchAddress("mumu_p4", &mumu_p4);
  treeR2->SetBranchAddress("ctpv", &ct);

  for(int i = 0; i < mEvt; i++) {
    treeR2->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	mumu_p4->Pt() > 50 && mumu_p4->Pt() < 66 &&
	mumu_p4->M() > 2.92 && mumu_p4->M() < 3.28 &&
	abs(mumu_p4->Rapidity()) < 1.2 &&
	ct > -0.005 && ct < 0.05)
      {
	dEta = muP_p4->Eta() - muM_p4->Eta();
	dPhi = muP_p4->Phi() - muM_p4->Phi();
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;

	dR = sqrt(dEta*dEta+1.4*1.4*dPhi*dPhi);
	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	
	if(dPhi < 0)
	  t_RpT->Fill(dR, dpT);
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with mid-pT trig MC" << endl; 
  }

  // high-pT
  mEvt = treeR3->GetEntries();
  perc = mEvt / 100;

  treeR3->SetBranchAddress("muP_p4", &muP_p4);
  treeR3->SetBranchAddress("muM_p4", &muM_p4);
  treeR3->SetBranchAddress("mumu_p4", &mumu_p4);
  treeR3->SetBranchAddress("ctpv", &ct);

  for(int i = 0; i < mEvt; i++) {
    treeR3->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	mumu_p4->Pt() > 66 && mumu_p4->Pt() < 120 &&
	mumu_p4->M() > 2.92 && mumu_p4->M() < 3.28 &&
	abs(mumu_p4->Rapidity()) < 1.2 &&
	ct > -0.005 && ct < 0.05)
      {
	dEta = muP_p4->Eta() - muM_p4->Eta();
	dPhi = muP_p4->Phi() - muM_p4->Phi();
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;
	
	dR = sqrt(dEta*dEta+1.4*1.4*dPhi*dPhi);
	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	
	if(dPhi < 0)
	  t_RpT->Fill(dR, dpT);
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with high-pT trig MC" << endl; 
  }

  finR2->Close();
  finR3->Close();

  r_RpT->GetXaxis()->SetTitle("#deltaR^{ell}");
  r_RpT->GetYaxis()->SetTitle("|#deltap_{T}|");
  t_RpT->GetXaxis()->SetTitle("#deltaR^{ell}");
  t_RpT->GetYaxis()->SetTitle("|#deltap_{T}|");
  
  TFile *fout = new TFile("RellpTstore_14.root", "recreate");
  r_RpT->Write();
  t_RpT->Write();
  fout->Close();
  
}
