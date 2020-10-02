void nrEvts()
{
  // open files and read TTrees
  TFile *fin1 = new TFile("../../Store_data_codes/2018/data18_cos.root");
  TTree *treeD = (TTree*)fin1->Get("data_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/2018/MC18_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../../Store_data_codes/2018/MC18_hpt_cos.root");
  TTree *treeM2 = (TTree*)fin3->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();

  cout << treeD->GetEntries("lts>2.5 && JpsiPt < 100") << " data events after cuts and " << treeM1->GetEntries("JpsiPt<46") + treeM2->GetEntries("JpsiPt<100") << " MC events after cuts" << endl;

  // prepare binning and histograms for plots
  const int nPtBins = 45;//68;
  double ptBins[nPtBins+1];
  for(int i=0 ; i<15; i++) ptBins[i] = (i+25.);
  for(int i=0; i<31; i++) ptBins[i+15] = (40.+2.*i);
  //for(int i=0; i<14; i++) ptBins[i+45] = (100+2.5*i);
  //for(int i=0; i<5; i++) ptBins[i+59] = (135. + 5.*i);
  //for(int i=0; i<5; i++) ptBins[i+64] = (160. + 10.*i);
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;
  
  TH1D *dataHist = new TH1D("dataH", "Data (NP)", nPtBins, ptBins);
  TH1D *mcHist = new TH1D("mcH", "MC", nPtBins, ptBins);
  TH1D *mcHist_test = new TH1D("mcH_t", "MC test", nPtBins, ptBins);

  // definitions to store data and MC events
  Double_t data_pt, lts;
  Double_t mc1_pt, mc2_pt;

  treeD->SetBranchAddress("JpsiPt", &data_pt);
  treeD->SetBranchAddress("lts", &lts);
  
  treeM1->SetBranchAddress("JpsiPt", &mc1_pt);
  treeM2->SetBranchAddress("JpsiPt", &mc2_pt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(lts > 2.5) {
	dataHist->Fill(data_pt);
      }
    }

  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if(mc1_pt < 46) {
	mcHist->Fill(mc1_pt);
      }
    }
 
  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
      mcHist->Fill(mc2_pt);
    }
  
  TFile *outfile = new TFile("files/nrEvts.root", "recreate");
  dataHist->Write();
  mcHist->Write();
  outfile->Close();

}
