#import "../ptbins.C"
#import "mbins.C"

// saving lifetime histos for all 3 mass regions
void bkgSave_N()
{
  // section for storing the lifetime histograms
  // prepare binning and histograms for plots
  int tbins = 55;
  double lowt = -0.05, hit = 0.5; // plotting in mm, not cm

  TH2D **ltHist = new TH2D*[nmBins];
  for(int i = 0; i < nmBins; i++) 
    ltHist[i] = new TH2D(Form("ltH_%s", lbl[i].c_str()), Form("2018 %s data c#tau", lbl[i].c_str()), tbins, lowt, hit, nPtBins, ptBins);
    
  // open and read the data tree
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/data18_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
    
  Double_t data_pt, data_lt, data_m, data_y;  
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  // cycle over data , fill the lifetime histograms
  int dEvt = treeD->GetEntries();
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      // filling in 8 mass intervals
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_y) < 1.2) 
	for(int j = 0; j < nmBins; j++) 
	  if(data_m < m_max[j] && data_m > m_min[j]) 
	    ltHist[j]->Fill(data_lt*10, data_pt); // filling with mm! Remember!
    }

  fin->Close();
    
  TFile *fout = new TFile("files/ltStore_N.root", "recreate");
  for(int i = 0; i < nmBins; i++) 
    ltHist[i]->Write();
  fout->Close();
  cout << "lifetime histograms all filled" << endl;
}
