#import "ltPerPt.C"
#import "plotLtPars.C"

int DO_FILL = 0;

void ltBkg()
{
  // prepare binning and histograms for plots 
  const int nPtBins = 7;
  double ptBins[nPtBins+1];
  for(int i=0; i<3; i++) ptBins[i] = 7.*i+25.;
  for(int i=0; i<4; i++) ptBins[i+3] = 46.+10.*i;
  ptBins[7] = 100;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  if(DO_FILL == 1) {
    int tbins = 120;
    double lowt = -0.1, hit = 0.5;
    TH1D **ltHist = new TH1D*[nPtBins];
    for(int ip = 0; ip < nPtBins; ip++) {
      ltHist[ip] = new TH1D(Form("ltH%.0f", ptBins[ip]), Form("2018 data c#tau (%.0f < p_{T} < %.0f GeV)", ptBins[ip], ptBins[ip+1]), tbins, lowt, hit);
    }
    
    // filling all the histos at once
    // start by getting the parameters for the sigma(y; pt) function
    ifstream ifile;
    string data;
    int pt_bins = 7;
    double pt_min[pt_bins], pt_max[pt_bins], c1[pt_bins], c2[pt_bins], m[pt_bins], aux, sig;
    ifile.open("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/2018/plots_massFit/fit_sig.txt");
    getline(ifile, data);
    getline(ifile, data);
    for(int i = 0; i < 7; i++) {
      ifile >> pt_min[i] >> pt_max[i] >> c1[i] >> aux >> c2[i] >> aux >> m[i] >> aux;
    }
    ifile.close();
    
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
	// filling flat mass SR (3.0 - 3.2 GeV)
	if(data_pt > 25 && data_pt < 100 && data_m < 3.2 && data_m > 3.0) {
	  for(int i_p = 0; i_p < nPtBins; i_p++)
	    if(data_pt > ptBins[i_p] && data_pt < ptBins[i_p+1])
	      ltHist[i_p]->Fill(data_lt*10);
	}
      }
    fin->Close();
    
    TFile *fout = new TFile("files/ltStore.root", "recreate");
    for(int i = 0; i < 3; i++) {
      for(int ip = 0; ip < nPtBins; ip++) {
	ltHist[ip]->Write();
      }
    }
    fout->Close();
    cout << "histograms all filled" << endl;
  }
  
  else {
    ofstream ftable;
    ftable.open("files/lt_fit.txt");
    ftable << "pt_min\t pt_max\t N_PR\t eN_PR\t N_NP\t eN_NP\t f\t ef\t mu\t emu\t sigma1\t esigma1\t sigma2\t esigma2\t lambda\t elambda\t chi2\t NDF\t f_NP\n";
    ftable.close();
    
    for(int i = 0; i < nPtBins; i++) {
      ltPerPt(ptBins[i], ptBins[i+1]);
    }
    plotLtPars();
  }
}
