// macro to plot the mass, rapidity and pT dists of data and MC
// also plots the lifetime of data

void store_ANdists_MC()
{
  // PART 1 : creating the histograms

  // two pT dists: PRSR data + the MC
  TH1D **h_pT = new TH1D*[4]; 
  for(int i = 0; i < 4; i++)
    h_pT[i] = new TH1D(Form("h_pT%d", i), Form("p_{T} distributions"), 100, 0, 200);

  // 2 y dists: (PRSR data + MC) over 1 pT region
  TH1D **h_y = new TH1D*[2]; 
  for(int i = 0; i < 2; i++)
    h_y[i] = new TH1D(Form("h_y%d", i), Form("y distributions"), 60, -1.5, 1.5);

  // 2 M dists: (PR data + MC) over 1 pT regions
  TH1D **h_m = new TH1D*[2]; 
  for(int i = 0; i < 2; i++)
    h_m[i] = new TH1D(Form("h_m%d", i), Form("M distributions"), 130, 3.35, 4.);

  // 1 lifetime dists: data in 1 pT regions
  // in micron!
  TH1D **h_lt = new TH1D*[1];
  for(int i = 0; i < 1; i++) 
    h_lt[i] = new TH1D(Form("h_lt%d", i), "c#tau distribution", 105, -150, 900);

  // 3 costh dists: PR data + NP data + MC over 1 pT regions
  TH1D **h_cos = new TH1D*[3]; 
  for(int i = 0; i < 3; i++)
    h_cos[i] = new TH1D(Form("h_cos%d", i), Form("cos#theta distributions"), 20, 0, 1.);

  double m_min[] = {3.4, 3.57, 3.82};
  double m_max[] = {3.52, 3.81, 4.0};
  
  // PART 2 : open files and read TTrees
  TFile *fin = new TFile("../../Store_data_codes/data17_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/MCm17_cos.root");
  TTree *treeM2 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../../Store_data_codes/MC17_cos.root");
  TTree *treeM3 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("../../Store_data_codes/MCvh17_cos.root");
  TTree *treeM4 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries();
  int m4Evt = treeM4->GetEntries();
  
  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m, data_y;
  Double_t mc_th, mc_pt, mc_lt, mc_m, mc_y;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
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

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);

      // fill the pT histo
      if(data_m > m_min[1] && data_m < m_max[1] && abs(data_lt) < 0.005 && abs(data_y) < 1.2) {
	h_pT[0]->Fill(data_pt);
      }

      // fill the y histo
      if(data_m > m_min[1] && data_m < m_max[1] && abs(data_lt) < 0.005) {
	if(data_pt > 20 && data_pt < 100) 
	  h_y[0]->Fill(data_y);
      }
      
      // fill the M histo
      if(abs(data_lt) < 0.005 && abs(data_y) < 1.2) {
	if(data_pt > 20 && data_pt < 100)
	  h_m[0]->Fill(data_m);
      }
      
      // fill the lt histo - in MICRON
      if(data_m > m_min[1] && data_m < m_max[1] && abs(data_y) < 1.2) {
	if(data_pt > 20 && data_pt < 100) 
	  h_lt[0]->Fill(data_lt*1e4);
      }

      // fill the costh histo - PRSR
      if(data_m > m_min[1] && data_m < m_max[1] && abs(data_y) < 1.2 && abs(data_lt) < 0.005) {
	if(data_pt > 25 && data_pt < 100) 
	  h_cos[0]->Fill(abs(cos(data_th)));
      }
      
      // fill the costh histo - NPSR
      if(data_m > m_min[1] && data_m < m_max[1] && abs(data_y) < 1.2 && data_lt > 0.01 && data_lt < 0.08) {
	if(data_pt > 20 && data_pt < 100) 
	  h_cos[1]->Fill(abs(cos(data_th)));
      }
      
    }
  
  cout << "data filled" << endl;
  
  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);

      // fill the pT histo
      if(mc_m > m_min[1] && mc_m < m_max[1] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	h_pT[1]->Fill(mc_pt);
      }

      // fill the y histo
      if(mc_m > m_min[1] && mc_m < m_max[1] && abs(mc_lt) < 0.005) {
	if(mc_pt > 20 && mc_pt < 100) 
	  h_y[1]->Fill(mc_y);
      }

      // fill the m histo
      if(abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	if(mc_pt > 20 && mc_pt < 100) 
	  h_m[1]->Fill(mc_m);
      }

      // fill the costh histo
      if(mc_m > m_min[1] && mc_m < m_max[1] && abs(mc_y) < 1.2 && abs(mc_lt) < 0.005) {
	if(mc_pt > 20 && mc_pt < 100) {
	  h_cos[2]->Fill(abs(cos(mc_th)));
	}
      }

    }

  cout << "MC filled" << endl;

  // just the pt of the other MC
  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);

      // fill the pT histo
      if(mc_m > m_min[1] && mc_m < m_max[1] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	h_pT[2]->Fill(mc_pt);
      }

    }

  cout << "MC filled" << endl;

  for(int i = 0; i < m4Evt; i++)
    {
      treeM4->GetEntry(i);

      // fill the pT histo
      if(mc_m > m_min[1] && mc_m < m_max[1] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	h_pT[3]->Fill(mc_pt);
      }

    }

  cout << "MC filled" << endl;

  
  fin->Close();
  fin2->Close();
  fin3->Close();
  fin4->Close();
  
  TFile *fout = new TFile("files/store_ANdists.root", "recreate");

  // store the pT dists
  string lbl_pt[] = {"Data", "MC", "lowPtMC", "highPtMC"};
  for(int i = 0; i < 4; i++) {
    h_pT[i]->Write(Form("h_pT_%s", lbl_pt[i].c_str()));
  }

  // store the y dists
  string lbl_y[] = {"Data", "MC"};
  for(int i = 0; i < 2; i++) {
    h_y[i]->Write(Form("h_y_%s", lbl_y[i].c_str()));
  }

  // store the M dists
  string lbl_m[] = {"Data", "MC"};
  for(int i = 0; i < 2; i++) {
    h_m[i]->Write(Form("h_m_%s", lbl_m[i].c_str()));
  }

  // store the lifetime dists
  string lbl_lt[] = {"Data"};
  for(int i = 0; i < 1; i++) {
    h_lt[i]->Write(Form("h_lt_%s", lbl_lt[i].c_str()));
  }

  // store the cos dists
  string lbl_cos[] = {"PR", "NP", "MC"};
  for(int i = 0; i < 3; i++) {
    h_cos[i]->Write(Form("h_cos_%s", lbl_cos[i].c_str()));
  }

  cout << "all histos stored" << endl;
  
  fout->Close();
}
