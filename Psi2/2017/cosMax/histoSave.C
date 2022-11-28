// code to get the 2d fine-binned data/mc ratio hist for |costh|max
// saves ratio (|costh|)

void histoSave()
{
  // fine pT binning - to be rebinned after background subtraction
  const int nPtBins = 20;
  double ptBins[nPtBins+1];
  for(int i = 0; i < 12; i++) ptBins[i] = 20 + 2.5*i;
  for(int i = 0; i < 6; i++) ptBins[i+12] = 50 + 5.*i;
  for(int i = 0; i < 3; i++) ptBins[i+18] = 80 + 10.*i;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  double m_min = 3.57, m_max = 3.81;

  // histograms for data and MC
  TH2D *dataHist_ab = new TH2D("dataH_ab", "2017 Data (PR)", 20, 0, 1., nPtBins, ptBins);
  TH2D *mcHist_ab = new TH2D("mcH_ab", "2017 MC", 20, 0, 1., nPtBins, ptBins);

  // open files and read TTrees
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/data17_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/MCm17_cos.root");
  TTree *treeM2 = (TTree*)fin2->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m2Evt = treeM2->GetEntries();
  
  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m;
  Double_t mc_th, mc_pt, mc_lt, mc_m;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  treeM2->SetBranchAddress("theta", &mc_th);
  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("lt", &mc_lt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_lt) < 0.005 && data_m > m_min && data_m < m_max) {      
	dataHist_ab->Fill(abs(cos(data_th)), data_pt);
      }
    }
  
  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
      if(mc_pt > ptBins[0] && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && mc_m > m_min && mc_m < m_max) {
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);	
      }
    }
  

  dataHist_ab->SetStats(0);
  dataHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  dataHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  mcHist_ab->SetStats(0);
  mcHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  mcHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  TH2D *ratioHist_ab = new TH2D("ratioH_ab", "2017 Data/MC", 20, 0, 1., nPtBins, ptBins);
  ratioHist_ab = (TH2D*)dataHist_ab->Clone("ratioH_ab");
  ratioHist_ab->Sumw2();
  ratioHist_ab->Divide(mcHist_ab);
  ratioHist_ab->SetTitle("2017 Data/MC");
    
  TFile *outfile = new TFile("histoStore.root", "recreate");
  dataHist_ab->Write();
  mcHist_ab->Write();
  ratioHist_ab->Write();
  outfile->Close();

}
