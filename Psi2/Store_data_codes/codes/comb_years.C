// macro to run through both year datafiles and combine

void comb_years()
{
  Double_t mmPt, th, phi, lt, lterr, mass, rap, dR;
  Double_t mPPt, mMPt, mPEta, mMEta;
  int tEvt, perc;
  
  // PART 1 : the data
  TFile *fData = new TFile("../dataS_cos.root", "recreate");
  TTree *tData = new TTree("data_cos", "");
  
  tData->Branch("theta", &th);
  tData->Branch("phi", &phi);
  tData->Branch("dimPt", &mmPt);
  tData->Branch("lt", &lt);
  tData->Branch("lterr", &lterr);
  tData->Branch("Mass", &mass);
  tData->Branch("Rap", &rap);
  tData->Branch("muonPPt", &mPPt);
  tData->Branch("muonPEta", &mPEta);
  tData->Branch("muonMPt", &mMPt);
  tData->Branch("muonMEta", &mMEta);
  tData->Branch("DeltaR", &dR);

  // open files and read TTrees - 2017
  TFile *finD7 = new TFile("../data17_cos.root");
  TTree *treeD7 = (TTree*)finD7->Get("data_cos");

  tEvt = treeD7->GetEntries();
  perc = tEvt / 100;
  treeD7->SetBranchAddress("dimPt", &mmPt);
  treeD7->SetBranchAddress("theta", &th);
  treeD7->SetBranchAddress("phi", &phi);
  treeD7->SetBranchAddress("lt", &lt);
  treeD7->SetBranchAddress("lterr", &lterr);
  treeD7->SetBranchAddress("Mass", &mass);
  treeD7->SetBranchAddress("Rap", &rap);
  treeD7->SetBranchAddress("muonPPt", &mPPt);
  treeD7->SetBranchAddress("muonPEta", &mPEta);
  treeD7->SetBranchAddress("muonMPt", &mMPt);
  treeD7->SetBranchAddress("muonMEta", &mMEta);
  treeD7->SetBranchAddress("DeltaR", &dR);

  for(int i = 0; i < tEvt; i++) {
    treeD7->GetEntry(i);
    tData->Fill();
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2017 data" << endl; 
  }
  finD7->Close();

  // open files and read TTrees - 2018
  TFile *finD8 = new TFile("../data18_cos.root");
  TTree *treeD8 = (TTree*)finD8->Get("data_cos");

  tEvt = treeD8->GetEntries();
  perc = tEvt / 100;
  treeD8->SetBranchAddress("dimPt", &mmPt);
  treeD8->SetBranchAddress("theta", &th);
  treeD8->SetBranchAddress("phi", &phi);
  treeD8->SetBranchAddress("lt", &lt);
  treeD8->SetBranchAddress("lterr", &lterr);
  treeD8->SetBranchAddress("Mass", &mass);
  treeD8->SetBranchAddress("Rap", &rap);
  treeD8->SetBranchAddress("muonPPt", &mPPt);
  treeD8->SetBranchAddress("muonPEta", &mPEta);
  treeD8->SetBranchAddress("muonMPt", &mMPt);
  treeD8->SetBranchAddress("muonMEta", &mMEta);
  treeD8->SetBranchAddress("DeltaR", &dR);

  for(int i = 0; i < tEvt; i++) {
    treeD8->GetEntry(i);
    tData->Fill();
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2018 data" << endl; 
  }
  finD8->Close();

  fData->Write();
  fData->Close();

  // PART 2 : the low-pt MC
  TFile *fMC = new TFile("../MCS_cos.root", "recreate");
  TTree *tMC = new TTree("MC_cos", "");
  
  tMC->Branch("theta", &th);
  tMC->Branch("phi", &phi);
  tMC->Branch("dimPt", &mmPt);
  tMC->Branch("lt", &lt);
  tMC->Branch("lterr", &lterr);
  tMC->Branch("Mass", &mass);
  tMC->Branch("Rap", &rap);
  tMC->Branch("muonPPt", &mPPt);
  tMC->Branch("muonPEta", &mPEta);
  tMC->Branch("muonMPt", &mMPt);
  tMC->Branch("muonMEta", &mMEta);
  tMC->Branch("DeltaR", &dR);
  
  // open files and read TTrees - 2017
  TFile *finM7 = new TFile("../MC17_cos.root");
  TTree *treeM7 = (TTree*)finM7->Get("MC_cos");

  tEvt = treeM7->GetEntries();
  perc = tEvt / 100;
  treeM7->SetBranchAddress("dimPt", &mmPt);
  treeM7->SetBranchAddress("theta", &th);
  treeM7->SetBranchAddress("phi", &phi);
  treeM7->SetBranchAddress("lt", &lt);
  treeM7->SetBranchAddress("lterr", &lterr);
  treeM7->SetBranchAddress("Mass", &mass);
  treeM7->SetBranchAddress("Rap", &rap);
  treeM7->SetBranchAddress("muonPPt", &mPPt);
  treeM7->SetBranchAddress("muonPEta", &mPEta);
  treeM7->SetBranchAddress("muonMPt", &mMPt);
  treeM7->SetBranchAddress("muonMEta", &mMEta);
  treeM7->SetBranchAddress("DeltaR", &dR);

  for(int i = 0; i < tEvt; i++) {
    treeM7->GetEntry(i);
    tMC->Fill();
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2017 MC (low-pT)" << endl; 
  }
  finM7->Close();

  // open files and read TTrees - 2018
  TFile *finM8 = new TFile("../MC18_cos.root");
  TTree *treeM8 = (TTree*)finM8->Get("MC_cos");

  tEvt = treeM8->GetEntries();
  perc = tEvt / 100;
  treeM8->SetBranchAddress("dimPt", &mmPt);
  treeM8->SetBranchAddress("theta", &th);
  treeM8->SetBranchAddress("phi", &phi);
  treeM8->SetBranchAddress("lt", &lt);
  treeM8->SetBranchAddress("lterr", &lterr);
  treeM8->SetBranchAddress("Mass", &mass);
  treeM8->SetBranchAddress("Rap", &rap);
  treeM8->SetBranchAddress("muonPPt", &mPPt);
  treeM8->SetBranchAddress("muonPEta", &mPEta);
  treeM8->SetBranchAddress("muonMPt", &mMPt);
  treeM8->SetBranchAddress("muonMEta", &mMEta);
  treeM8->SetBranchAddress("DeltaR", &dR);

  for(int i = 0; i < tEvt; i++) {
    treeM8->GetEntry(i);
    tMC->Fill();
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2018 MC (low-pT)" << endl; 
  }
  finM8->Close();

  fMC->Write();
  fMC->Close();

  // PART 4 : the high-pt MC
  TFile *fMCvh = new TFile("../MCvhS_cos.root", "recreate");
  TTree *tMCvh = new TTree("MC_cos", "");
  
  tMCvh->Branch("theta", &th);
  tMCvh->Branch("phi", &phi);
  tMCvh->Branch("dimPt", &mmPt);
  tMCvh->Branch("lt", &lt);
  tMCvh->Branch("lterr", &lterr);
  tMCvh->Branch("Mass", &mass);
  tMCvh->Branch("Rap", &rap);
  tMCvh->Branch("muonPPt", &mPPt);
  tMCvh->Branch("muonPEta", &mPEta);
  tMCvh->Branch("muonMPt", &mMPt);
  tMCvh->Branch("muonMEta", &mMEta);
  tMCvh->Branch("DeltaR", &dR);

  // open files and read TTrees - 2017
  TFile *finMvh7 = new TFile("../MCvh17_cos.root");
  TTree *treeMvh7 = (TTree*)finMvh7->Get("MC_cos");

  tEvt = treeMvh7->GetEntries();
  perc = tEvt / 100;
  treeMvh7->SetBranchAddress("dimPt", &mmPt);
  treeMvh7->SetBranchAddress("theta", &th);
  treeMvh7->SetBranchAddress("phi", &phi);
  treeMvh7->SetBranchAddress("lt", &lt);
  treeMvh7->SetBranchAddress("lterr", &lterr);
  treeMvh7->SetBranchAddress("Mass", &mass);
  treeMvh7->SetBranchAddress("Rap", &rap);
  treeMvh7->SetBranchAddress("muonPPt", &mPPt);
  treeMvh7->SetBranchAddress("muonPEta", &mPEta);
  treeMvh7->SetBranchAddress("muonMPt", &mMPt);
  treeMvh7->SetBranchAddress("muonMEta", &mMEta);
  treeMvh7->SetBranchAddress("DeltaR", &dR);

  for(int i = 0; i < tEvt; i++) {
    treeMvh7->GetEntry(i);
    tMCvh->Fill();
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2017 MC (high-pT)" << endl; 
  }
  finMvh7->Close();

  // open files and read TTrees - 2018
  TFile *finMvh8 = new TFile("../MCvh18_cos.root");
  TTree *treeMvh8 = (TTree*)finMvh8->Get("MC_cos");

  tEvt = treeMvh8->GetEntries();
  perc = tEvt / 100;
  treeMvh8->SetBranchAddress("dimPt", &mmPt);
  treeMvh8->SetBranchAddress("theta", &th);
  treeMvh8->SetBranchAddress("phi", &phi);
  treeMvh8->SetBranchAddress("lt", &lt);
  treeMvh8->SetBranchAddress("lterr", &lterr);
  treeMvh8->SetBranchAddress("Mass", &mass);
  treeMvh8->SetBranchAddress("Rap", &rap);
  treeMvh8->SetBranchAddress("muonPPt", &mPPt);
  treeMvh8->SetBranchAddress("muonPEta", &mPEta);
  treeMvh8->SetBranchAddress("muonMPt", &mMPt);
  treeMvh8->SetBranchAddress("muonMEta", &mMEta);
  treeMvh8->SetBranchAddress("DeltaR", &dR);

  for(int i = 0; i < tEvt; i++) {
    treeMvh8->GetEntry(i);
    tMCvh->Fill();
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2018 MC (high-pT)" << endl; 
  }
  finMvh8->Close();

  fMCvh->Write();
  fMCvh->Close();
}
