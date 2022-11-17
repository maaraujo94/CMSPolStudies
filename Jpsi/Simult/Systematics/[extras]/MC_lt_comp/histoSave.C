// code to get the 2d data/mc ratio hist (pr, np, peak and sidebands)
// saves data, mc, ratio

void histoSave()
{
  const int nPtBins = 18;
  double ptBins[nPtBins+1] = {25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 55, 60, 65, 70, 75, 80, 90, 100};

// histograms for MC (PR and NP)
  TH2D *PRHist  = new TH2D("PRH", "MC (PR)",        20, 0, 1., nPtBins, ptBins);
  TH2D *NPHist  = new TH2D("NPH", "MC (NP)",        20, 0, 1., nPtBins, ptBins);
  TH2D *dHist   = new TH2D("dH",  "Data (NP Peak)", 20, 0, 1., nPtBins, ptBins);

  // open files and read TTrees
  TFile *finD = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/data18_cos.root");
  TTree *treeD = (TTree*)finD->Get("data_cos");

  TFile *finNP = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCNP18_cos.root");
  TTree *treeNP = (TTree*)finNP->Get("MC_cos");

  TFile *fin1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MC18_cos.root");
  TTree *treeM1 = (TTree*)fin1->Get("MC_cos");
  TFile *fin2 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCm18_cos.root");
  TTree *treeM2 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCh18_cos.root");
  TTree *treeM3 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCvh18_cos.root");
  TTree *treeM4 = (TTree*)fin4->Get("MC_cos");

  int dEvt = treeD->GetEntries();
  int NPEvt = treeNP->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries();
  int m4Evt = treeM4->GetEntries(); 
 
  // definitions to store data and MC events
  Double_t mc_th, mc_pt, mc_lt, mc_m, mc_y;
  
  treeD->SetBranchAddress("theta", &mc_th);
  treeD->SetBranchAddress("dimPt", &mc_pt);
  treeD->SetBranchAddress("Rap", &mc_y);
  treeD->SetBranchAddress("Mass", &mc_m);
  treeD->SetBranchAddress("lt", &mc_lt);

  treeNP->SetBranchAddress("theta", &mc_th);
  treeNP->SetBranchAddress("dimPt", &mc_pt);
  treeNP->SetBranchAddress("Rap", &mc_y);
  treeNP->SetBranchAddress("Mass", &mc_m);
  treeNP->SetBranchAddress("lt", &mc_lt);
  
  treeM1->SetBranchAddress("theta", &mc_th);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("lt", &mc_lt);

  treeM2->SetBranchAddress("theta", &mc_th);
  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Rap", &mc_y);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("lt", &mc_lt);

  treeM3->SetBranchAddress("theta", &mc_th);
  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Rap", &mc_y);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("lt", &mc_lt);

  treeM4->SetBranchAddress("theta", &mc_th);
  treeM4->SetBranchAddress("dimPt", &mc_pt);
  treeM4->SetBranchAddress("Rap", &mc_y);
  treeM4->SetBranchAddress("Mass", &mc_m);
  treeM4->SetBranchAddress("lt", &mc_lt);

  // cycle over NP data
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      // pt and y conditions
      if(mc_pt > ptBins[0] && mc_pt < ptBins[nPtBins] && mc_lt > 0.01 && mc_lt < 0.05 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	dHist->Fill(abs(cos(mc_th)), mc_pt);
      }
    }                                                                          

  // cycle over MC, fill the costh histogram acc to binning
  for(int i = 0; i < NPEvt; i++)
    {
      treeNP->GetEntry(i);
      // pt and y conditions
      if(mc_pt > ptBins[0] && mc_pt < ptBins[nPtBins] && mc_lt > 0.01 && mc_lt < 0.05 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	NPHist->Fill(abs(cos(mc_th)), mc_pt);
      }
    }
  
  // MC sample 1
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if(mc_pt > ptBins[0] && mc_pt < 45 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	PRHist->Fill(abs(cos(mc_th)), mc_pt);
      }
    }

  // MC sample 2
  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
      if(mc_pt > 45 && mc_pt < 50 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	PRHist->Fill(abs(cos(mc_th)), mc_pt);
      }
    }

  // MC sample 3
  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);
      if(mc_pt > 50 && mc_pt < 70 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	PRHist->Fill(abs(cos(mc_th)), mc_pt);
      }
    }

  // MC sample 4
  for(int i = 0; i < m4Evt; i++)
    {
      treeM4->GetEntry(i);
      if(mc_pt > 70 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	PRHist->Fill(abs(cos(mc_th)), mc_pt);
      }
    }

  // get ratio
  TH2D *rHist  = new TH2D();
  TH2D *rMCHist  = new TH2D();

  rHist = (TH2D*)dHist->Clone("rH");
  rHist->Sumw2();
  rHist->Divide(NPHist);
  rHist->SetTitle("Data/MC (NP)");
  
  rMCHist = (TH2D*)PRHist->Clone("rMCH");
  rMCHist->Sumw2();
  rMCHist->Divide(NPHist);
  rMCHist->SetTitle("PR/NP (MC)");

  // store data, MC and data/MC ratios
  TFile *outfile = new TFile("mcComp.root", "recreate");
  dHist->Write();
  PRHist->Write();
  NPHist->Write();
  rHist->Write();
  rMCHist->Write();
  outfile->Close();
  
}
