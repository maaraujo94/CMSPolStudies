#import "../../../ptbins.C"

void bkgSave()
{
  // setup the histos
  TH1D *h_LSB = new TH1D("lsbH", "Run 2 Data (LSB)", nPtBins, ptBins);
  TH1D *h_RSB = new TH1D("rsbH", "Run 2 Data (RSB)", nPtBins, ptBins);
  TH1D *h_tot = new TH1D("totH", "Run 2 Data (SR)", nPtBins, ptBins);

  TH1D *hNP_LSB = new TH1D("lsbH_NP", "Run 2 Data (NP LSB)", nPtBins, ptBins);
  TH1D *hNP_RSB = new TH1D("rsbH_NP", "Run 2 Data (NP RSB)", nPtBins, ptBins);
  TH1D *hNP_tot = new TH1D("totH_NP", "Run 2 Data (NP SR)", nPtBins, ptBins);

  // open file, get data
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/dataS_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
    
  int dEvt = treeD->GetEntries();
  
  // definitions to store data events
  Double_t data_pt, data_lt, data_m;

  // using same mass limit as mass fit (excludes first 3 mass bins)
  double m_min[] = {3.43, 3.57, 3.82};
  double m_max[] = {3.52, 3.81, 4.0};

  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
    
  // cycle over data, fill the histograms
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      // PR
      if(abs(data_lt) < 0.005) {
	// SR
	if(data_m < m_max[1] && data_m > m_min[1]) {
	  h_tot->Fill(data_pt);
	}
	// LSB
	else if(data_m < m_max[0] && data_m > m_min[0]) {
	  h_LSB->Fill(data_pt);
	}
	// RSB
	else if(data_m < m_max[2] && data_m > m_min[2]) {
	  h_RSB->Fill(data_pt);
	}
      }
      // NP
      if(data_lt > 0.01 && data_lt < 0.08) {
	// SR
	if(data_m < m_max[1] && data_m > m_min[1]) {
	  hNP_tot->Fill(data_pt);
	}
	// LSB
	else if(data_m < m_max[0] && data_m > m_min[0]) {
	  hNP_LSB->Fill(data_pt);
	}
	// RSB
	else if(data_m < m_max[2] && data_m > m_min[2]) {
	  hNP_RSB->Fill(data_pt);
	}
      }
    }
    
  TFile *fout = new TFile("files/yieldHistos.root", "recreate");
  h_LSB->Write();
  h_RSB->Write();
  h_tot->Write();
  hNP_LSB->Write();
  hNP_RSB->Write();
  hNP_tot->Write();
  fout->Close();
}
