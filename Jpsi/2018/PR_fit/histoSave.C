// code to get the 2d fine-binned data/mc ratio hist
// plots data, mc and ratio
// saves ratio (normal and |costh|), number of entries in each histo bin

void histoSave()
{
  // fine pT binning - to be rebinned after background subtraction
  const int nPtBins = 38;
  double ptBins[nPtBins+1];
  for(int i = 0; i < 12; i++) ptBins[i] = 25+i;
  for(int i = 0; i < 6; i++) ptBins[i+12] = 37 + 1.5*i;
  for(int i = 0; i < 12; i++) ptBins[i+18] = 46 + 2.5*i;
  for(int i = 0; i < 6; i++) ptBins[i+30] = 76 + 4*i;
  for(int i = 0; i < 3; i++) ptBins[i+36] = 100 + 10*i;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // histograms for data and MC
  TH2D *dataHist = new TH2D("dataH", "2018 Data (PR)", 40, -1., 1., nPtBins, ptBins);
  TH2D *mcHist = new TH2D("mcH", "2018 MC", 40, -1., 1., nPtBins, ptBins);
  
  TH2D *dataHist_ab = new TH2D("dataH_ab", "2018 Data (PR)", 20, 0, 1., nPtBins, ptBins);
  TH2D *mcHist_ab = new TH2D("mcH_ab", "2018 MC", 20, 0, 1., nPtBins, ptBins);

  // open files and read TTrees
  TFile *fin = new TFile("../../Store_data_codes/data18_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/MC18_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../../Store_data_codes/MCh18_cos.root");
  TTree *treeM2 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("../../Store_data_codes/MCvh18_cos.root");
  TTree *treeM3 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries();
  
  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m, data_y;
  Double_t mc_th, mc_pt, mc_lt, mc_m, mc_y;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  treeM1->SetBranchAddress("theta", &mc_th);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("lt", &mc_lt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_lt) < 0.01 && abs(data_y) < 1.2 && data_m > 3.0 && data_m < 3.2) {
      
	dataHist->Fill(cos(data_th), data_pt);
	dataHist_ab->Fill(abs(cos(data_th)), data_pt);
       
      }
    }
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if(mc_pt > ptBins[0] && mc_pt < 46 && abs(mc_lt) < 0.01 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {

	mcHist->Fill(cos(mc_th), mc_pt);
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
      }
    }

  treeM2->SetBranchAddress("theta", &mc_th);
  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Rap", &mc_y);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("lt", &mc_lt);

  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
      if(mc_pt > 46 && mc_pt < 66 && abs(mc_lt) < 0.01 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {

	mcHist->Fill(cos(mc_th), mc_pt);
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
	
      }
    }

  treeM3->SetBranchAddress("theta", &mc_th);
  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Rap", &mc_y);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("lt", &mc_lt);
  
  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);
      if(mc_pt > 66 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.01 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {

	mcHist->Fill(cos(mc_th), mc_pt);
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
	
      }
    }
  
  dataHist_ab->SetStats(0);
  dataHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  dataHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  mcHist_ab->SetStats(0);
  mcHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  mcHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  TFile *outfile = new TFile("files/histoStore.root", "recreate");
  dataHist->Write();
  dataHist_ab->Write();
  mcHist->Write();
  mcHist_ab->Write();
  outfile->Close();

  cout << dataHist->GetEntries() << " data events and " << mcHist->GetEntries() << " MC events" << endl;

  // split between pT ranges - must set by hand
  cout << dataHist->GetYaxis()->GetBinLowEdge(1) << " " << dataHist->GetYaxis()->GetBinUpEdge(18) << endl;
  cout << dataHist->GetYaxis()->GetBinLowEdge(19) << " " << dataHist->GetYaxis()->GetBinUpEdge(26) << endl;
  cout << dataHist->GetYaxis()->GetBinLowEdge(27) << " " << dataHist->GetYaxis()->GetBinUpEdge(38) << endl;
  
  cout << dataHist->Integral(1, 18) << " data events and " << mcHist->Integral(1, 18) << " MC events in pT range 1" << endl;
  cout << dataHist->Integral(19, 26) << " data events and " << mcHist->Integral(19, 26) << " MC events in pT range 2" << endl;
  cout << dataHist->Integral(27, 38) << " data events and " << mcHist->Integral(27, 38) << " MC events in pT range 3" << endl;  

}
