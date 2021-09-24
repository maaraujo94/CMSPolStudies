void ltSave()
{
  // section for storing the lifetime histograms
  // prepare binning and histograms for plots
  const int nPtBins = 1;
  double ptBins[nPtBins+1];
  ptBins[0] = 25;
  ptBins[1] = 120;
 
  int nReg = 3;
  string regN[] = {"SR", "LSB", "RSB"};
    
  int tbins = 110;
  double lowt = -0.05, hit = 0.5; // plotting in mm, not cm
  TH1D **ltHist = new TH1D*[nReg];
  for(int ip = 0; ip < nReg; ip++) {
    ltHist[ip] = new TH1D(Form("ltH%s", regN[ip].c_str()), Form("2018 data %s c#tau (%.0f < p_{T} < %.0f GeV)", regN[ip].c_str(), ptBins[0], ptBins[nPtBins]), tbins, lowt, hit);
  }
  
  // open and read the data tree
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/data18_cos.root");
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
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins]) {
	// filling flat mass SR (3.0 - 3.2 GeV)
	if(data_m < 3.2 && data_m > 3.0 && abs(data_y) < 1.2) 
	  ltHist[0]->Fill(data_lt*10); // filling with mm! Remember!
	
	// filling flat mass LSB (2.92 - 2.95 GeV)
	if(data_m < 2.95 && data_m > 2.92 && abs(data_y) < 1.2) 
	  ltHist[1]->Fill(data_lt*10); // filling with mm! Remember!
	  
	// filling flat mass RSB (3.21 - 3.28 GeV)
	if(data_m < 3.28 && data_m > 3.21 && abs(data_y) < 1.2)
	  ltHist[2]->Fill(data_lt*10); // filling with mm! Remember!
	
      }
    }
  fin->Close();
    
  TFile *fout = new TFile("files/ltStore.root", "recreate");
  for(int ip = 0; ip < nReg; ip++) {
    ltHist[ip]->Write();
  }
  fout->Close();
}
