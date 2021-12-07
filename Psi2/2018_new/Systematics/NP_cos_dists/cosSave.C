void cosSave()
{
  // prepare binning and histograms for plots
  const int nPtBins = 7;
  double ptBins[nPtBins+1];
  for(int i=0; i<3; i++) ptBins[i] = 7.*i+25.;
  for(int i=0; i<4; i++) ptBins[i+3] = 46.+10.*i;
  ptBins[7] = 120;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // histograms for data (peak and NP) and MC
  TH2D *NP1Hist_ab = new TH2D("NP1H_ab", "2018 Data (NP 1)", 20, 0, 1., nPtBins, ptBins);
  TH2D *NP2Hist_ab = new TH2D("NP2H_ab", "2018 Data (NP 2)", 20, 0, 1., nPtBins, ptBins);
  TH2D *mcHist_ab = new TH2D("mcH_ab", "2018 MC", 20, 0, 1., nPtBins, ptBins);

  // open files and read TTrees
  TFile *fin = new TFile("../../../Store_data_codes/data18_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  
  int dEvt = treeD->GetEntries();
  
  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m, data_y;

  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  // cycle over data, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && data_m > 3.57 && data_m < 3.81 && abs(data_y) < 1.2) {
	if(data_lt > 0.01 && data_lt < 0.02 ) {
	  NP1Hist_ab->Fill(abs(cos(data_th)), data_pt);
	}
	else if(data_lt > 0.03 && data_lt < 0.05 ) {
	  NP2Hist_ab->Fill(abs(cos(data_th)), data_pt);
	}
      }
    }

  NP1Hist_ab->SetStats(0);
  NP1Hist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  NP1Hist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  NP2Hist_ab->SetStats(0);
  NP2Hist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  NP2Hist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");


  TFile *outfile = new TFile("cosStore.root", "recreate");
  NP1Hist_ab->Write();
  NP2Hist_ab->Write();
  outfile->Close();
  
}
