// code to get the 2d fine-binned data/mc ratio hist for |costh|max
// saves ratio (|costh|)

void histoSave()
{
  // fine pT binning - to be rebinned after background subtraction
  const int nPtBins = 49;
  double ptBins[nPtBins+1];
  for(int i = 0; i < 15; i++) ptBins[i] = 25+i;
  for(int i = 0; i < 30; i++) ptBins[i+15] = 40 + 2.*i;
  for(int i = 0; i < 5; i++) ptBins[i+45] = 100+5.*i;
  // new limits for the MC
  ptBins[18] = 45;
  ptBins[19] = 47.5;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // histograms for data and MC
  TH2D *dataHist_ab = new TH2D("dataH_ab", "Run 2 Data (PR)", 20, 0, 1., nPtBins, ptBins);
  TH2D *mcHist_ab = new TH2D("mcH_ab", "Run 2 MC", 20, 0, 1., nPtBins, ptBins);

  // open files and read TTrees
  TFile *fin = new TFile("../../Store_data_codes/dataS_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  
  TFile *fin1 = new TFile("../../Store_data_codes/MCS_cos.root");
  TTree *treeM1 = (TTree*)fin1->Get("MC_cos");
  
  TFile *fin2 = new TFile("../../Store_data_codes/MCmS_cos.root");
  TTree *treeM2 = (TTree*)fin2->Get("MC_cos");
  
  TFile *fin3 = new TFile("../../Store_data_codes/MChS_cos.root");
  TTree *treeM3 = (TTree*)fin3->Get("MC_cos");
  
  TFile *fin4 = new TFile("../../Store_data_codes/MCvhS_cos.root");
  TTree *treeM4 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries(); 
  int m4Evt = treeM4->GetEntries();
 
  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m;
  Double_t mc_th, mc_pt, mc_lt, mc_m;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  treeM1->SetBranchAddress("theta", &mc_th);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("lt", &mc_lt);

  treeM2->SetBranchAddress("theta", &mc_th);
  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("lt", &mc_lt);

  treeM3->SetBranchAddress("theta", &mc_th);
  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("lt", &mc_lt);

  treeM4->SetBranchAddress("theta", &mc_th);
  treeM4->SetBranchAddress("dimPt", &mc_pt);
  treeM4->SetBranchAddress("Mass", &mc_m);
  treeM4->SetBranchAddress("lt", &mc_lt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_lt) < 0.005 && data_m > 3.0 && data_m < 3.2) {      
	dataHist_ab->Fill(abs(cos(data_th)), data_pt);
      }
    }
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if(mc_pt > ptBins[0] && mc_pt < 45 && abs(mc_lt) < 0.005 && mc_m > 3.0 && mc_m < 3.2) {
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
      }
    }

  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
      if(mc_pt > 45 && mc_pt < 50 && abs(mc_lt) < 0.005 && mc_m > 3.0 && mc_m < 3.2) {
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);	
      }
    }

  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);
      if(mc_pt > 50 && mc_pt < 70 && abs(mc_lt) < 0.005 && mc_m > 3.0 && mc_m < 3.2) {
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);	
      }
    }
  
  for(int i = 0; i < m4Evt; i++)
    {
      treeM4->GetEntry(i);
      if(mc_pt > 70 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && mc_m > 3.0 && mc_m < 3.2) {
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
      }
    }
  
  dataHist_ab->SetStats(0);
  dataHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  dataHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  mcHist_ab->SetStats(0);
  mcHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  mcHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  TH2D *ratioHist_ab = new TH2D("ratioH_ab", "Run 2 Data/MC", 20, 0, 1., nPtBins, ptBins);
  ratioHist_ab = (TH2D*)dataHist_ab->Clone("ratioH_ab");
  ratioHist_ab->Sumw2();
  ratioHist_ab->Divide(mcHist_ab);
  ratioHist_ab->SetTitle("Run 2 Data/MC");
    
  TFile *outfile = new TFile("histoStore.root", "recreate");
  dataHist_ab->Write();
  mcHist_ab->Write();
  ratioHist_ab->Write();
  outfile->Close();

}
