// code that runs over MC, applies all cuts and saves a root file with the appropriate angular histos

const double gPI = TMath::Pi();

double R_eff(double dEta, double dPhi)
{
  double R_val = 0;
  if(dEta < 0.06) R_val = dPhi-0.12;
  else if(dEta < 0.12) R_val = dPhi + 1.5*dEta - 0.21;
  else R_val = dPhi;

  return R_val;
}

void storeRd()
{
  int mEvt, perc;
  Double_t dEta, dPhi, dR, dpT;
  Float_t ct;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;

  // output histos
  TH2D** r_RpT = new TH2D*[3];
  TH2D** t_RpT = new TH2D*[3];
  for(int i = 0; i < 3; i++) {
    r_RpT[i] = new TH2D(Form("r_RpT_%d", i), Form("reco p_{T} > 50 (%d)", i+1), 150, -0.3, 0.3, 75, 0, 120);
    t_RpT[i] = new TH2D(Form("t_RpT_%d", i), Form("trig p_{T} > 50 (%d)", i+1), 150, -0.3, 0.3, 75, 0, 120);
  }
  
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

	if(abs(dEta) < 0.06)
	  r_RpT[0]->Fill(abs(dPhi)-0.12, dpT);
	else if(abs(dEta) < 0.12)
	  r_RpT[1]->Fill(abs(dPhi)+1.5*abs(dEta)-0.21, dpT);
	else
	  r_RpT[2]->Fill(abs(dPhi), dpT);
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

	if(abs(dEta) < 0.06)
	  r_RpT[0]->Fill(abs(dPhi)-0.12, dpT);
	else if(abs(dEta) < 0.12)
	  r_RpT[1]->Fill(abs(dPhi)+1.5*abs(dEta)-0.21, dpT);
	else
	  r_RpT[2]->Fill(abs(dPhi), dpT);
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

	if(abs(dEta) < 0.06)
	  t_RpT[0]->Fill(abs(dPhi)-0.12, dpT);
	else if(abs(dEta) < 0.12)
	  t_RpT[1]->Fill(abs(dPhi)+1.5*abs(dEta)-0.21, dpT);
	else
	  t_RpT[2]->Fill(abs(dPhi), dpT);
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

	if(abs(dEta) < 0.06)
	  t_RpT[0]->Fill(abs(dPhi)-0.12, dpT);
	else if(abs(dEta) < 0.12)
	  t_RpT[1]->Fill(abs(dPhi)+1.5*abs(dEta)-0.21, dpT);
	else
	  t_RpT[2]->Fill(abs(dPhi), dpT);
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with high-pT trig MC" << endl; 
  }

  finR2->Close();
  finR3->Close();

  for(int i = 0; i < 3; i++) {
    r_RpT[i]->GetXaxis()->SetTitle("#deltaR");
    r_RpT[i]->GetYaxis()->SetTitle("|#deltap_{T}|");
    t_RpT[i]->GetXaxis()->SetTitle("#deltaR");
    t_RpT[i]->GetYaxis()->SetTitle("|#deltap_{T}|");
  }
  
  TFile *fout = new TFile("Rdstore.root", "recreate");
  for(int i = 0; i < 3; i++) {
    r_RpT[i]->Write();
    t_RpT[i]->Write();
  }
  fout->Close();
 
}
