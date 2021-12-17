// code that runs over MC, applies all cuts and saves a root file with the appropriate angular histos

const double gPI = TMath::Pi();

void storeEtaPt()
{
  int mEvt, perc;
  Double_t dEta, dPt;
  Float_t ct;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;

  // output histos
  TH2D* rC_EtaPt = new TH2D("rC_EtaPt", "reco p_{T} > 50", 75, -0.3, 0.3, 75, 0, 120);
  TH2D* tC_EtaPt = new TH2D("tC_EtaPt", "trig p_{T} > 50", 75, -0.3, 0.3, 75, 0, 120);
  TH2D* r_EtaPt = new TH2D(Form("r_EtaPt"), Form("reco p_{T} > 50"), 75, 0, 0.3, 75, 0, 120);;
  TH2D* t_EtaPt = new TH2D(Form("t_EtaPt"), Form("trig p_{T} > 50"), 75, 0, 0.3, 75, 0, 120);
  
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
	  dEta = muP_p4->Eta() - muM_p4->Eta();
	}
	else {
	  dEta = muM_p4->Eta() - muP_p4->Eta();
	}
	dPt = abs(muP_p4->Pt() - muM_p4->Pt());
	
	rC_EtaPt->Fill(dEta, dPt);
	r_EtaPt->Fill(abs(dEta), dPt);
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
	  dEta = muP_p4->Eta() - muM_p4->Eta();
	}
	else {
	  dEta = muM_p4->Eta() - muP_p4->Eta();
	}
	dPt = abs(muP_p4->Pt() - muM_p4->Pt());
	
	rC_EtaPt->Fill(dEta, dPt);
	r_EtaPt->Fill(abs(dEta), dPt);
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
	  dEta = muP_p4->Eta() - muM_p4->Eta();
	}
	else {
	  dEta = muM_p4->Eta() - muP_p4->Eta();
	}
	dPt = abs(muP_p4->Pt() - muM_p4->Pt());
	
	tC_EtaPt->Fill(dEta, dPt);
	t_EtaPt->Fill(abs(dEta), dPt);
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
	  dEta = muP_p4->Eta() - muM_p4->Eta();
	}
	else {
	  dEta = muM_p4->Eta() - muP_p4->Eta();
	}
	dPt = abs(muP_p4->Pt() - muM_p4->Pt());
	
	tC_EtaPt->Fill(dEta, dPt);
	t_EtaPt->Fill(abs(dEta), dPt);
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with high-pT trig MC" << endl; 
  }

  finR2->Close();
  finR3->Close();

  r_EtaPt->GetXaxis()->SetTitle("|#delta#eta|");
  r_EtaPt->GetYaxis()->SetTitle("|#deltap_{T}|");
  t_EtaPt->GetXaxis()->SetTitle("|#delta#eta|");
  t_EtaPt->GetYaxis()->SetTitle("|#deltap_{T}|");
  
  rC_EtaPt->GetXaxis()->SetTitle("#delta#eta");
  rC_EtaPt->GetYaxis()->SetTitle("|#deltap_{T}|");
  tC_EtaPt->GetXaxis()->SetTitle("#delta#eta");
  tC_EtaPt->GetYaxis()->SetTitle("|#deltap_{T}|");
  
  TFile *fout = new TFile("EtaPtstore.root", "recreate");
  rC_EtaPt->Write();
  tC_EtaPt->Write();
  r_EtaPt->Write();
  t_EtaPt->Write();
  fout->Close();
  
}
