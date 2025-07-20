#import "../ptbins.C"
#import "../ptcut.C"

// saving lifetime histos for all 3 mass regions
void bkgSave()
{
  // section for storing the lifetime histograms
  // prepare binning and histograms for plots
  int tbins = 55;
  double lowt = -0.05, hit = 0.5; // plotting in mm, not cm
  double m_min[] = {3.4, 3.57, 3.82};
  double m_max[] = {3.52, 3.81, 4.0};

  TH2D **ltHist = new TH2D*[3];
  string lbl[] = {"LSB", "SR", "RSB"};
  for(int i = 0; i < 3; i++) 
    ltHist[i] = new TH2D(Form("ltH_%s", lbl[i].c_str()), Form("Run 2 %s data c#tau", lbl[i].c_str()), tbins, lowt, hit, nPtBins, ptBins);
    
  // open and read the data tree
  TFile *fin = new TFile("../../Store_data_codes/dataS_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
    
  Double_t data_pt, data_lt, data_m, data_y;  
double mPPt, mMPt, mPEta, mMEta;

  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  treeD->SetBranchAddress("muonPEta", &mPEta);
  treeD->SetBranchAddress("muonMEta", &mMEta);
  treeD->SetBranchAddress("muonPPt", &mPPt);
  treeD->SetBranchAddress("muonMPt", &mMPt);

  // cycle over data , fill the lifetime histograms
  int dEvt = treeD->GetEntries();
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      // filling in 3 mass intervals
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_y) < 1.2) 
	if((abs(mPEta) > eta_lim || mPPt > pt_cut) && (abs(mMEta) > eta_lim || mMPt > pt_cut))
	  for(int j = 0; j < 3; j++) 
	    if(data_m < m_max[j] && data_m > m_min[j]) 
	      ltHist[j]->Fill(data_lt*10, data_pt); // filling with mm! Remember!
    }
  
  fin->Close();
  
  TFile *fout = new TFile("files/ltStore.root", "recreate");
  for(int i = 0; i < 3; i++) 
    ltHist[i]->Write();
  fout->Close();
  cout << "lifetime histograms all filled" << endl;
}

 
