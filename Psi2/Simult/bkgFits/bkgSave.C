#import "../ptbins.C"

// macro to save data mass distributions
void bkgSave()
{
  // section for storing the mass histograms
  // prepare histograms for plots - fine
  int mbins = 65;
  double lowm = 3.35, him = 4.0;
  TH2D *h_d2d = new TH2D("mH", "Run 2 PR data M(#mu#mu)", mbins, lowm, him, nPtBins, ptBins);
  TH2D *hNP_d2d = new TH2D("mH_NP", "Run 2 NP data M(#mu#mu)", mbins, lowm, him, nPtBins, ptBins);

  // filling all the histos at once    
  // open and read the data tree
  TFile *fin1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/dataS_cos.root");
  TTree *tree1 = (TTree*)fin1->Get("data_cos");
    
  // data
  Double_t data_pt, data_lt, data_m, data_y;  
    
  tree1->SetBranchAddress("dimPt", &data_pt);
  tree1->SetBranchAddress("Rap", &data_y);
  tree1->SetBranchAddress("Mass", &data_m);
  tree1->SetBranchAddress("lt", &data_lt);
    
  // cycle over data , fill the lifetime histogram
  int dEvt = tree1->GetEntries();
  for(int i = 0; i < dEvt; i++)
    {
      tree1->GetEntry(i);
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_lt) < 0.005 && abs(data_y) < 1.2) {
	h_d2d->Fill(data_m, data_pt);
      }
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && data_lt > 0.01 && data_lt < 0.05 && abs(data_y) < 1.2) {
	hNP_d2d->Fill(data_m, data_pt);
      }
    }
  fin1->Close();

  TFile *fout = new TFile("files/mStore.root", "recreate");
  h_d2d->Write();
  hNP_d2d->Write();
  fout->Close();
  cout << "mass histograms all filled" << endl;
}
