// code to get the 2d data/mc ratio hist (pr, np, peak and sidebands)
// saves data, mc, ratio

#import "../../ptbins.C"

void histoSave()
{
  // histograms for "data" (reweighed MC) and unpol MC
  TH2D *MCpHist = new TH2D("MCpH", "MC (+0.4 #lambda)", 20, 0, 1., nPtBins, ptBins);
  TH2D *MCmHist = new TH2D("MCmH", "MC (-0.1 #lambda)", 20, 0, 1., nPtBins, ptBins);
  TH2D *MCHist  = new TH2D("MCH",  "MC",                20, 0, 1., nPtBins, ptBins);

  // open files and read TTrees
  TFile *fin1 = new TFile("../../../Store_data_codes/MCS_cos.root");
  TTree *treeM1 = (TTree*)fin1->Get("MC_cos");
  TFile *fin2 = new TFile("../../../Store_data_codes/MCmS_cos.root");
  TTree *treeM2 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../../../Store_data_codes/MChS_cos.root");
  TTree *treeM3 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("../../../Store_data_codes/MCvhS_cos.root");
  TTree *treeM4 = (TTree*)fin4->Get("MC_cos");
  
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries();
  int m4Evt = treeM4->GetEntries(); 
 
  // definitions to store data and MC events
  Double_t mc_th, mc_pt, mc_lt, mc_m, mc_y;
  double w_p, w_m;
  
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

  // cycle over MC, filling unpol + reweighed histos
  // MC sample 1
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if(mc_pt > ptBins[0] && mc_pt < 45 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	w_p = 1+0.4*cos(mc_th)*cos(mc_th);
	w_m = 1-0.1*cos(mc_th)*cos(mc_th);
	
	MCHist->Fill(abs(cos(mc_th)), mc_pt);
	MCpHist->Fill(abs(cos(mc_th)), mc_pt, w_p);
	MCmHist->Fill(abs(cos(mc_th)), mc_pt, w_m);
      }
    }

  // MC sample 2
  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
      if(mc_pt > 45 && mc_pt < 50 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	w_p = 1+0.4*cos(mc_th)*cos(mc_th);
	w_m = 1-0.1*cos(mc_th)*cos(mc_th);
	
	MCHist->Fill(abs(cos(mc_th)), mc_pt);
	MCpHist->Fill(abs(cos(mc_th)), mc_pt, w_p);
	MCmHist->Fill(abs(cos(mc_th)), mc_pt, w_m);
      }
    }

  // MC sample 3
  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);
      if(mc_pt > 50 && mc_pt < 70 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	w_p = 1+0.4*cos(mc_th)*cos(mc_th);
	w_m = 1-0.1*cos(mc_th)*cos(mc_th);
	
	MCHist->Fill(abs(cos(mc_th)), mc_pt);
	MCpHist->Fill(abs(cos(mc_th)), mc_pt, w_p);
	MCmHist->Fill(abs(cos(mc_th)), mc_pt, w_m);
      }
    }

  // MC sample 4
  for(int i = 0; i < m4Evt; i++)
    {
      treeM4->GetEntry(i);
      if(mc_pt > 70 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	w_p = 1+0.4*cos(mc_th)*cos(mc_th);
	w_m = 1-0.1*cos(mc_th)*cos(mc_th);
	
	MCHist->Fill(abs(cos(mc_th)), mc_pt);
	MCpHist->Fill(abs(cos(mc_th)), mc_pt, w_p);
	MCmHist->Fill(abs(cos(mc_th)), mc_pt, w_m);
      }
    }

  // get ratios
  TH2D *rMCpHist = new TH2D();
  TH2D *rMCmHist = new TH2D();

  rMCpHist = (TH2D*)MCpHist->Clone("rMCpH");
  rMCpHist->Sumw2();
  rMCpHist->Divide(MCHist);
  rMCpHist->SetTitle("Data/MC (+0.4 #lambda)");
  
  rMCmHist = (TH2D*)MCmHist->Clone("rMCmH");
  rMCmHist->Sumw2();
  rMCmHist->Divide(MCHist);
  rMCmHist->SetTitle("Data/MC (-0.1 #lambda)");

  // store data, MC and data/MC ratios
  TFile *outfile = new TFile("files/histoStore.root", "recreate");
  MCHist->Write();
  MCpHist->Write();
  MCmHist->Write();
  rMCpHist->Write();
  rMCmHist->Write();
  outfile->Close();
}
