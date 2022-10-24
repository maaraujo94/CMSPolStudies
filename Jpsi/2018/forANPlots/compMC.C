// macro to compare the costheta dists of diff MC samples in same pT range

void compMC()
{
  // pT limits - all comparisons
  int nB[] = {3, 3, 6};
  int nBL = 3;
  double ptBinsL[] = {40, 42, 44, 46};
  int nBM = 3;
  double ptBinsM[] = {46, 48, 50, 52};
  int nBH = 6;
  double ptBinsH[] = {66, 72, 78, 84, 90, 100, 120};

  // histograms defined in pairs for comparison
  TH2D **h_MC = new TH2D*[6];
  h_MC[0] = new TH2D("mcH0", "Low pT MC", 20, 0, 1., nBL, ptBinsL);
  h_MC[1] = new TH2D("mcH1", "Mid pT MC", 20, 0, 1., nBL, ptBinsL);
  h_MC[2] = new TH2D("mcH2", "Mid pT MC", 20, 0, 1., nBM, ptBinsM);
  h_MC[3] = new TH2D("mcH3", "High pT MC", 20, 0, 1., nBM, ptBinsM);
  h_MC[4] = new TH2D("mcH4", "High pT MC", 20, 0, 1., nBH, ptBinsH);
  h_MC[5] = new TH2D("mcH5", "Very High pT MC", 20, 0, 1., nBH, ptBinsH);

  // open files and read TTrees
  TFile *fin0 = new TFile("../../Store_data_codes/MC18_cos.root");
  TTree *treeM0 = (TTree*)fin0->Get("MC_cos");
  TFile *fin1 = new TFile("../../Store_data_codes/MCm18_cos.root");
  TTree *treeM1 = (TTree*)fin1->Get("MC_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/MCh18_cos.root");
  TTree *treeM2 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../../Store_data_codes/MCvh18_cos.root");
  TTree *treeM3 = (TTree*)fin3->Get("MC_cos");
  
  int m0Evt = treeM0->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries();
  
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

  // start cycling through trees
  for(int i = 0; i < m0Evt; i++)
    {
      treeM0->GetEntry(i);

      if(mc_pt > 40 && mc_pt < 46 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) 
	h_MC[0]->Fill(abs(cos(mc_th)), mc_pt);
    }

  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);

      if(mc_pt > 40 && mc_pt < 46 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) 
	h_MC[1]->Fill(abs(cos(mc_th)), mc_pt);

      else if(mc_pt > 46 && mc_pt < 52 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) 
	h_MC[2]->Fill(abs(cos(mc_th)), mc_pt);
    }

  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);

      if(mc_pt > 46 && mc_pt < 52 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) 
	h_MC[3]->Fill(abs(cos(mc_th)), mc_pt);

      else if(mc_pt > 66 && mc_pt < 120 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) 
	h_MC[4]->Fill(abs(cos(mc_th)), mc_pt);
    }

  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);

      if(mc_pt > 66 && mc_pt < 120 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) 
	h_MC[5]->Fill(abs(cos(mc_th)), mc_pt);
    }

  fin0->Close();
  fin1->Close();
  fin2->Close();
  fin3->Close();

  // get the ratio between histos
  TH2D **hr_MC = new TH2D*[3];
  string lbl_MC[3] = {"low", "mid", "high"};
  for(int i = 0; i < 3; i++) {
    hr_MC[i] = new TH2D();
    hr_MC[i] = (TH2D*)h_MC[i*2]->Clone(Form("ratio_%d", i));
    hr_MC[i]->Sumw2();
    hr_MC[i]->Divide(h_MC[i*2+1]);
    hr_MC[i]->SetTitle(Form("%s p_{T} comparison", lbl_MC[i].c_str()));
  }

  TFile *outfile = new TFile("files/MCcompStore.root", "recreate");
  for(int i = 0; i < 3; i++) {
    h_MC[2*i]->Write();
    h_MC[2*i+1]->Write();
    hr_MC[i]->Write();
  }
  outfile->Close();
}
