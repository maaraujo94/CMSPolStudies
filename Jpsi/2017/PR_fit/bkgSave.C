#import "../ptbins.C"

void bkgSave()
{
  // section for storing the lifetime histograms
  // prepare binning and histograms for plots
  int tbins = 55;
  double lowt = -0.05, hit = 0.5; // plotting in mm, not cm
  TH2D *ltHist = new TH2D("ltH", "2017 data c#tau", tbins, lowt, hit, nPtBins, ptBins);
    
  // open and read the data tree
  TFile *fin = new TFile("../../Store_data_codes/data17_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
    
  Double_t data_pt, data_lt, data_m, data_y;  
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  // cycle over data , fill the lifetime histogram
  int dEvt = treeD->GetEntries();
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
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

 
