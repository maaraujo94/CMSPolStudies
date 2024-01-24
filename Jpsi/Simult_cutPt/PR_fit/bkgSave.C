#import "../ptcut.C"

#import "../ptbins.C"

void bkgSave()
{
  // section for storing the lifetime histograms
  // prepare binning and histograms for plots
  int tbins = 55;
  double lowt = -0.05, hit = 0.5; // plotting in mm, not cm
  TH2D *ltHist = new TH2D("ltH", "Run 2 data c#tau", tbins, lowt, hit, nPtBins, ptBins);
    
  // open and read the data tree
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/dataS_cos.root");
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
  
  // cycle over data , fill the lifetime histogram
  int dEvt = treeD->GetEntries();
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
if((abs(mPEta) > eta_lim || mPPt > pt_cut) && (abs(mMEta) > eta_lim || mMPt > pt_cut))
      // filling flat mass SR (3.0 - 3.2 GeV)
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && data_m < 3.2 && data_m > 3.0 && abs(data_y) < 1.2) {
	ltHist->Fill(data_lt*10, data_pt); // filling with mm! Remember!
      }
    }
  fin->Close();
    
  TFile *fout = new TFile("files/ltStore.root", "recreate");
  ltHist->Write();
  fout->Close();
  cout << "lifetime histograms all filled" << endl;
}

 
