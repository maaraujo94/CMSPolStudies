#import "../ptbins.C"
#import "mbins.C"
#import "../etacut.C"

void bkgSave()
{
  // section for storing the lifetime histograms
  // prepare binning and histograms for plots
  int tbins = 55;
  double lowt = -0.05, hit = 0.5; // plotting in mm, not cm

  TH1D **ltHist = new TH1D*[nmBins];
  for(int i = 0; i < nmBins; i++) 
    ltHist[i] = new TH1D(Form("ltH_%s", lbl[i].c_str()), Form("Run 2 %s data c#tau", lbl[i].c_str()), tbins, lowt, hit);
  
  // open and read the data tree
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/dataS_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
    
  Double_t data_pt, data_lt, data_m, data_y;  
  double mPEta, mMEta;

  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  treeD->SetBranchAddress("muonPEta", &mPEta);
  treeD->SetBranchAddress("muonMEta", &mMEta);

  // cycle over data , fill the lifetime histogram
  int dEvt = treeD->GetEntries();
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      // filling in the SB mass intervals
      if((abs(mPEta) < eta_lo || abs(mPEta) > eta_hi) && (abs(mMEta) < eta_lo || abs(mMEta) > eta_hi))
	if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_y) < 1.2) 
	  for(int j = 0; j < nmBins; j++) 
	    if(data_m < m_max[j] && data_m > m_min[j]) 
	      ltHist[j]->Fill(data_lt*10); // filling with mm! Remember!
    }
  fin->Close();
  
  TFile *fout = new TFile("files/ltStore.root", "recreate");
  for(int i = 0; i < nmBins; i++) 
    ltHist[i]->Write();
 fout->Close();
  cout << "lifetime histograms all filled" << endl;

}
