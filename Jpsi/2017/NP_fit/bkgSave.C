int DO_MASS = 1;

void bkgSave()
{
  // section for storing the mass histograms
  if(DO_MASS == 1) {

    // prepare binning and histograms for plots
    const int nPtBins = 9;
    double ptBins[nPtBins+1];
    for(int i=0; i<5; i++) ptBins[i] = 5.*i+25.;
    for(int i=0; i<4; i++) ptBins[i+5] = 50+10.*i;
    ptBins[9] = 120;
    for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
    cout << endl;

    // prepare mass histograms
    int mbins = 40;
    double lowm = 2.9, him = 3.3;
    TH2D *h_d2d = new TH2D("mH", "2017 data M(#mu#mu)", mbins, lowm, him, nPtBins, ptBins);
        
    // filling all the histos at once    
    // open and read the data tree
    TFile *fin1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/data17_cos.root");
    TTree *tree1 = (TTree*)fin1->Get("data_cos");
    
    // data
    Double_t data_pt, data_lt, data_m, data_y;  
    
    tree1->SetBranchAddress("dimPt", &data_pt);
    tree1->SetBranchAddress("Rap", &data_y);
    tree1->SetBranchAddress("Mass", &data_m);
    tree1->SetBranchAddress("lt", &data_lt);
    
    // cycle over data , fill the mass histogram
    int dEvt = tree1->GetEntries();
    for(int i = 0; i < dEvt; i++)
      {
	tree1->GetEntry(i);
	if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && data_lt > 0.01 && data_lt < 0.05 && abs(data_y) < 1.2) {
	  h_d2d->Fill(data_m, data_pt);
	}
      }
    fin1->Close();

    TFile *fout = new TFile("files/mStore.root", "recreate");
    h_d2d->Write();	
    fout->Close();
    cout << "mass histograms all filled" << endl;
  }
  
}
