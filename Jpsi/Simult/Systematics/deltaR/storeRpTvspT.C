// code that runs over MC, applies all cuts and saves a root file with the appropriate angular histos

const double gPI = TMath::Pi();

void storeRpTvspT()
{
  int mEvt, perc;
  Double_t dEta, dPhi, dR, dpT, pT;
  Float_t ct;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;

  // output histos
  TH2D* r_RpTpT = new TH2D("r_RpTpT", "reco", 100, 0, 0.4, 95, 25, 120);
  TH2D* t_RpTpT = new TH2D("t_RpTpT", "trig", 100, 0, 0.4, 95, 25, 120);
  
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

  for(int i = 0; i < mEvt; i++) {
    treeNT1->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	mumu_p4->Pt() > 25 && mumu_p4->Pt() < 46 &&
	mumu_p4->M() > 2.92 && mumu_p4->M() < 3.28 &&
	abs(mumu_p4->Rapidity()) < 1.2 &&
	ct > -0.005 && ct < 0.05)
      {
	if(muP_p4->Pt() > muM_p4->Pt()) {
	  dPhi = muP_p4->Phi() - muM_p4->Phi();
	  dEta = muP_p4->Eta() - muM_p4->Eta();
	}
	else {
	  dPhi = muM_p4->Phi() - muP_p4->Phi();
	  dEta = muM_p4->Eta() - muP_p4->Eta();
	}
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;

	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;
	pT = mumu_p4->Pt();
	
	r_RpTpT->Fill(dR, pT);
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with low-pT reco MC" << endl; 
  }

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
	mumu_p4->Pt() > 46 && mumu_p4->Pt() < 66 &&
	mumu_p4->M() > 2.92 && mumu_p4->M() < 3.28 &&
	abs(mumu_p4->Rapidity()) < 1.2 &&
	ct > -0.005 && ct < 0.05)
      {
	if(muP_p4->Pt() > muM_p4->Pt()) {
	  dPhi = muP_p4->Phi() - muM_p4->Phi();
	  dEta = muP_p4->Eta() - muM_p4->Eta();
	}
	else {
	  dPhi = muM_p4->Phi() - muP_p4->Phi();
	  dEta = muM_p4->Eta() - muP_p4->Eta();
	}
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;

	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;
	pT = mumu_p4->Pt();
	
	r_RpTpT->Fill(dR, pT);
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
	if(muP_p4->Pt() > muM_p4->Pt()) {
	  dPhi = muP_p4->Phi() - muM_p4->Phi();
	  dEta = muP_p4->Eta() - muM_p4->Eta();
	}
	else {
	  dPhi = muM_p4->Phi() - muP_p4->Phi();
	  dEta = muM_p4->Eta() - muP_p4->Eta();
	}
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;

	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;
	pT = mumu_p4->Pt();
	
	r_RpTpT->Fill(dR, pT);
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with high-pT reco MC" << endl; 
  }

  finNT1->Close();
  finNT2->Close();
  finNT3->Close();

  // reco trees
  TFile *finR1 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-mc-LOCAL18.root");
  TTree *treeR1 = (TTree*)finR1->Get("jpsitree");
  TFile *finR2 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-mc-LOCAL18-highpt.root");
  TTree *treeR2 = (TTree*)finR2->Get("jpsitree");
  TFile *finR3 = new TFile("/eos/user/m/maaraujo/JpsiRun2/MC/filtered-all-psi-mc-LOCAL18-veryhighpt.root");
  TTree *treeR3 = (TTree*)finR3->Get("jpsitree");
  
  // low-pT
  mEvt = treeR1->GetEntries();
  perc = mEvt / 100;

  treeR1->SetBranchAddress("muP_p4", &muP_p4);
  treeR1->SetBranchAddress("muM_p4", &muM_p4);
  treeR1->SetBranchAddress("mumu_p4", &mumu_p4);
  treeR1->SetBranchAddress("ctpv", &ct);

  for(int i = 0; i < mEvt; i++) {
    treeR1->GetEntry(i);
    if( muP_p4->Pt() > 5.6 && muM_p4->Pt() > 5.6 && 
	abs(muP_p4->Eta()) < 1.4 && abs(muM_p4->Eta()) < 1.4 &&
	mumu_p4->Pt() > 25 && mumu_p4->Pt() < 46 &&
	mumu_p4->M() > 2.92 && mumu_p4->M() < 3.28 &&
	abs(mumu_p4->Rapidity()) < 1.2 &&
	ct > -0.005 && ct < 0.05)
      {
	if(muP_p4->Pt() > muM_p4->Pt()) {
	  dPhi = muP_p4->Phi() - muM_p4->Phi();
	  dEta = muP_p4->Eta() - muM_p4->Eta();
	}
	else {
	  dPhi = muM_p4->Phi() - muP_p4->Phi();
	  dEta = muM_p4->Eta() - muP_p4->Eta();
	}
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;

	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;
	pT = mumu_p4->Pt();
	
	t_RpTpT->Fill(dR, pT);
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with low-pT trig MC" << endl; 
  }

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
	mumu_p4->Pt() > 46 && mumu_p4->Pt() < 66 &&
	mumu_p4->M() > 2.92 && mumu_p4->M() < 3.28 &&
	abs(mumu_p4->Rapidity()) < 1.2 &&
	ct > -0.005 && ct < 0.05)
      {
	if(muP_p4->Pt() > muM_p4->Pt()) {
	  dPhi = muP_p4->Phi() - muM_p4->Phi();
	  dEta = muP_p4->Eta() - muM_p4->Eta();
	}
	else {
	  dPhi = muM_p4->Phi() - muP_p4->Phi();
	  dEta = muM_p4->Eta() - muP_p4->Eta();
	}
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;

	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;
	pT = mumu_p4->Pt();
	
	t_RpTpT->Fill(dR, pT);
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
	if(muP_p4->Pt() > muM_p4->Pt()) {
	  dPhi = muP_p4->Phi() - muM_p4->Phi();
	  dEta = muP_p4->Eta() - muM_p4->Eta();
	}
	else {
	  dPhi = muM_p4->Phi() - muP_p4->Phi();
	  dEta = muM_p4->Eta() - muP_p4->Eta();
	}
	if(dPhi > gPI) dPhi -= 2.*gPI;
	if(dPhi < -gPI) dPhi += 2.*gPI;
	
	dpT = abs(muP_p4->Pt()-muM_p4->Pt());
	dR = sqrt(dEta*dEta+dPhi*dPhi)+log(dpT)/45.;
	pT = mumu_p4->Pt();
	
	t_RpTpT->Fill(dR, pT);
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with high-pT trig MC" << endl; 
  }

  finR1->Close();
  finR2->Close();
  finR3->Close();

  r_RpTpT->GetXaxis()->SetTitle("#deltaR_{pt}");
  r_RpTpT->GetYaxis()->SetTitle("p_{T}");
  t_RpTpT->GetXaxis()->SetTitle("#deltaR_{pt}");
  t_RpTpT->GetYaxis()->SetTitle("p_{T}");
  
  TFile *fout = new TFile("RpTpTstore.root", "recreate");
  r_RpTpT->Write();
  t_RpTpT->Write();
  fout->Close();
  
}
