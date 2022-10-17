#import "../rcut.C"

int DO_LT = 1;
int DO_MASS = 1;

void bkgSave()
{
  // section for storing the lifetime histograms
  if(DO_LT == 1) {
    // prepare binning and histograms for plots
    const int nPtBins = 19;
    double ptBins[nPtBins+1];
    for(int i = 0; i < 10; i++) ptBins[i] = 25 + 2.5*i;
    for(int i = 0; i < 6; i++) ptBins[i+10] = 50 + 5.*i;
    for(int i = 0; i < 2; i++) ptBins[i+16] = 80 + 10.*i;
    for(int i = 0; i < 2; i++) ptBins[i+18] = 100 + 20.*i;
    for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
    cout << endl;
    
    int tbins = 110;
    double lowt = -0.05, hit = 0.5; // plotting in mm, not cm
    TH2D *ltHist = new TH2D("ltH", "Run 2 data c#tau", tbins, lowt, hit, nPtBins, ptBins);
    
    // open and read the data tree
    TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/dataS_cos.root");
    TTree *treeD = (TTree*)fin->Get("data_cos");
    
    Double_t data_pt, data_lt, data_m, data_y;  
Double_t dR;
    treeD->SetBranchAddress("dimPt", &data_pt);
    treeD->SetBranchAddress("Rap", &data_y);
    treeD->SetBranchAddress("Mass", &data_m);
    treeD->SetBranchAddress("lt", &data_lt);
    treeD->SetBranchAddress("DeltaR", &dR);
  
    // cycle over data , fill the lifetime histogram
    int dEvt = treeD->GetEntries();
    for(int i = 0; i < dEvt; i++)
      {
	treeD->GetEntry(i);
	if(dR > r_cut)
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
    TH2D *h_d2d = new TH2D("mH", "Run 2 data M(#mu#mu)", mbins, lowm, him, nPtBins, ptBins);
        
    // filling all the histos at once    
    // open and read the data tree
    TFile *fin1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/dataS_cos.root");
    TTree *tree1 = (TTree*)fin1->Get("data_cos");
    
    // data
    Double_t data_pt, data_lt, data_m, data_y;  
    Double_t dR;
    
    tree1->SetBranchAddress("dimPt", &data_pt);
    tree1->SetBranchAddress("Rap", &data_y);
    tree1->SetBranchAddress("Mass", &data_m);
    tree1->SetBranchAddress("lt", &data_lt);
    tree1->SetBranchAddress("DeltaR", &dR);
    
    // cycle over data , fill the lifetime histogram
    int dEvt = tree1->GetEntries();
    for(int i = 0; i < dEvt; i++)
      {
	tree1->GetEntry(i);
	if(dR > r_cut)
	  if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_lt) < 0.005 && abs(data_y) < 1.2) {
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
