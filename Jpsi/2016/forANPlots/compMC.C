// macro to compare the costheta dists of diff MC samples in same pT range

#import "../ptbins.C"

void compMC()
{
  // histograms defined in pairs for comparison
  TH2D **h_MC = new TH2D*[2];
  h_MC[0] = new TH2D("mcH0", "MC", 20, 0, 1., nPtBins, ptBins);
  h_MC[1] = new TH2D("mcH1", "pre MC", 20, 0, 1., nPtBins, ptBins);

  // open files and read TTrees
  TFile *fin0 = new TFile("../../Store_data_codes/MC16_cos.root");
  TTree *treeM0 = (TTree*)fin0->Get("MC_cos");
  TFile *fin1 = new TFile("../../Store_data_codes/MC16p_cos.root");
  TTree *treeM1 = (TTree*)fin1->Get("MC_cos");
  
  int m0Evt = treeM0->GetEntries();
  int m1Evt = treeM1->GetEntries();
  
  // definitions to store MC events
  Double_t mc_th, mc_pt, mc_lt, mc_m, mc_y;
  
  treeM0->SetBranchAddress("theta", &mc_th);
  treeM0->SetBranchAddress("dimPt", &mc_pt);
  treeM0->SetBranchAddress("Rap", &mc_y);
  treeM0->SetBranchAddress("Mass", &mc_m);
  treeM0->SetBranchAddress("lt", &mc_lt);

  treeM1->SetBranchAddress("theta", &mc_th);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("lt", &mc_lt);

  // start cycling through trees
  for(int i = 0; i < m0Evt; i++)
    {
      treeM0->GetEntry(i);

      if(mc_pt > ptBins[0] && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) 
	h_MC[0]->Fill(abs(cos(mc_th)), mc_pt);
    }

  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);

      if(mc_pt > ptBins[0] && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) 
	h_MC[1]->Fill(abs(cos(mc_th)), mc_pt);
    }

  fin0->Close();
  fin1->Close();

  // get the ratio between histos
  TH2D *hr_MC = new TH2D();
  hr_MC = (TH2D*)h_MC[0]->Clone("MCr");
  hr_MC->Sumw2();
  hr_MC->Divide(h_MC[1]);
  hr_MC->SetTitle(Form("MC comparison"));
  
  TFile *outfile = new TFile("files/MCcompStore.root", "recreate");
  h_MC[0]->Write();
  h_MC[1]->Write();
  hr_MC->Write();
  outfile->Close();
}
