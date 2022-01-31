// code that runs over MC, applies all cuts and saves a root file with the appropriate angular histos

const double gPI = TMath::Pi();

void storePhi()
{
  int mEvt, perc;
  Double_t phiP, phiM, phi1, phi2;
  Float_t ct;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;

  // output histos
  TH1D* r_PhiD = new TH1D(Form("r_PhiD"), Form("reco p_{T} > 50"), 100, -0.25, 0.25);
  TH1D* t_PhiD = new TH1D(Form("t_PhiD"), Form("trig p_{T} > 50"), 100, -0.25, 0.25);
  TH1D* r_PhiPM = new TH1D(Form("r_PhiPM"), Form("reco p_{T} > 50"), 100, -0.25, 0.25);
  TH1D* t_PhiPM = new TH1D(Form("t_PhiPM"), Form("trig p_{T} > 50"), 100, -0.25, 0.25);
  
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
	phiP = muP_p4->Phi();
	phiM = muM_p4->Phi();

	if(muP_p4->Pt() > muM_p4->Pt()) {
	  phi1 = phiP;
	  phi2 = phiM;
	}
	else {
	  phi1 = phiM;
	  phi2 = phiP;
	}
	  
	r_PhiD->Fill(phi1-phi2);
	r_PhiPM->Fill(phiP-phiM);
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
	phiP = muP_p4->Phi();
	phiM = muM_p4->Phi();

	if(muP_p4->Pt() > muM_p4->Pt()) {
	  phi1 = phiP;
	  phi2 = phiM;
	}
	else {
	  phi1 = phiM;
	  phi2 = phiP;
	}
	  
	r_PhiD->Fill(phi1-phi2);
	r_PhiPM->Fill(phiP-phiM);
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
	phiP = muP_p4->Phi();
	phiM = muM_p4->Phi();

	if(muP_p4->Pt() > muM_p4->Pt()) {
	  phi1 = phiP;
	  phi2 = phiM;
	}
	else {
	  phi1 = phiM;
	  phi2 = phiP;
	}
	  
	t_PhiD->Fill(phi1-phi2);
	t_PhiPM->Fill(phiP-phiM);
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
	phiP = muP_p4->Phi();
	phiM = muM_p4->Phi();

	if(muP_p4->Pt() > muM_p4->Pt()) {
	  phi1 = phiP;
	  phi2 = phiM;
	}
	else {
	  phi1 = phiM;
	  phi2 = phiP;
	}
	  
	t_PhiD->Fill(phi1-phi2);
	t_PhiPM->Fill(phiP-phiM);
     }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with high-pT trig MC" << endl; 
  }

  finR2->Close();
  finR3->Close();

  r_PhiPM->GetXaxis()->SetTitle("#phi^{+}-#phi^{-}");
  r_PhiD->GetXaxis()->SetTitle("#phi^{h}-#phi^{l}");
  t_PhiPM->GetXaxis()->SetTitle("#phi^{+}-#phi^{-}");
  t_PhiD->GetXaxis()->SetTitle("#phi^{h}-#phi^{l}");

  TFile *fout = new TFile("Phistore_d.root", "recreate");
  r_PhiPM->Write();
  r_PhiD->Write();
  t_PhiPM->Write(); 
  t_PhiD->Write();
  fout->Close();
  
}
