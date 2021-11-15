int DO_LT = 1;
int DO_MASS = 1;

void bkgSave()
{
  // section for storing the lifetime histograms
  if(DO_LT == 1) {
    // prepare binning and histograms for plots
    const int nPtBins = 7;
    double ptBins[nPtBins+1];
    for(int i=0; i<3; i++) ptBins[i] = 7.*i+25.;
    for(int i=0; i<4; i++) ptBins[i+3] = 46.+10.*i;
    ptBins[7] = 120;
    for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
    cout << endl;

    double m_min = 3.57, m_max = 3.81;
    
    int tbins = 110;
    double lowt = -0.05, hit = 0.5; // plotting in mm, not cm
    TH1D **ltHist = new TH1D*[nPtBins];
    for(int ip = 0; ip < nPtBins; ip++) {
      ltHist[ip] = new TH1D(Form("ltH%.0f", ptBins[ip]), Form("2017 data c#tau (%.0f < p_{T} < %.0f GeV)", ptBins[ip], ptBins[ip+1]), tbins, lowt, hit);
    }
    
    // open and read the data tree
    TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/data17_cos.root");
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
	// filling flat mass SR (3.57 - 3.81 GeV)
	if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && data_m < m_max && data_m > m_min && abs(data_y) < 1.2) {
	  for(int i_p = 0; i_p < nPtBins; i_p++)
	    if(data_pt > ptBins[i_p] && data_pt < ptBins[i_p+1])
	      ltHist[i_p]->Fill(data_lt*10); // filling with mm! Remember!
	}
      }
    fin->Close();
    
    TFile *fout = new TFile("files/ltStore.root", "recreate");
    for(int ip = 0; ip < nPtBins; ip++) {
      ltHist[ip]->Write();
    }
    fout->Close();
    cout << "lifetime histograms all filled" << endl;
  }

  // section for storing the mass histograms
  if(DO_MASS == 1) {

    // prepare binning and histograms for plots
    const int nPtBins = 7;
    double ptBins[nPtBins+1];
    for(int i=0; i<3; i++) ptBins[i] = 7.*i+25.;
    for(int i=0; i<4; i++) ptBins[i+3] = 46.+10.*i;
    ptBins[7] = 120;
    for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
    cout << endl;

    // prepare mass histograms
    int mbins = 65;
    double lowm = 3.35, him = 4.0;
    TH1D **h_d1d = new TH1D*[nPtBins];
    for(int ip = 0; ip < nPtBins; ip++)
      h_d1d[ip] = new TH1D(Form("mH%.0f", ptBins[ip]), Form("2017 data M(#mu#mu) (%.0f < p_{T} < %.0f GeV)",  ptBins[ip], ptBins[ip+1]), mbins, lowm, him);
        
    // filling all the histos at once    
    // open and read the data tree
    TFile *fin1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/data17_cos.root");
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
	  for(int i_p = 0; i_p < nPtBins; i_p++)
	    if(data_pt > ptBins[i_p] && data_pt < ptBins[i_p+1])
	      h_d1d[i_p]->Fill(data_m);
	}
      }
    fin1->Close();

    TFile *fout = new TFile("files/mStore.root", "recreate");
    for(int ip = 0; ip < nPtBins; ip++) {
      h_d1d[ip]->Write();	
    }
    fout->Close();
    cout << "mass histograms all filled" << endl;
  }
  
}
