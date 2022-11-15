// macro to plot the mass, rapidity and pT dists of data and MC
// also plots the lifetime of data

void store_ANdists()
{
  // PART 1 : creating the histograms

  // 2 pT dists: PRSR data + the MC
  TH1D **h_pT = new TH1D*[2]; 
  for(int i = 0; i < 2; i++)
    h_pT[i] = new TH1D(Form("h_pT%d", i), Form("p_{T} distributions"), 100, 0, 200);

  // 2 y dists: (PRSR data + MC) over 1 pT region
  TH1D **h_y = new TH1D*[2]; 
  for(int i = 0; i < 2; i++)
    h_y[i] = new TH1D(Form("h_y%d", i), Form("y distributions"), 60, -1.5, 1.5);

  // 2 M dists: (PR data + MC) over 1 pT region
  TH1D **h_m = new TH1D*[2]; 
  for(int i = 0; i < 2; i++)
    h_m[i] = new TH1D(Form("h_m%d", i), Form("M distributions"), 80, 3.35, 4.0);

  // 1 lifetime dist: data
  // in micron!
  TH1D **h_lt = new TH1D*[1];
  for(int i = 0; i < 1; i++) 
    h_lt[i] = new TH1D(Form("h_lt%d", i), "c#tau distribution", 75, -150, 600);

  // 3 costh dists: PR data + NP data + MC over 1 pT region
  TH1D **h_cos = new TH1D*[3]; 
  for(int i = 0; i < 3; i++)
    h_cos[i] = new TH1D(Form("h_cos%d", i), Form("cos#theta distributions"), 20, 0, 1.);
  
  // PART 2 : open files and read TTrees
  TFile *fin = new TFile("../../Store_data_codes/data18_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2a = new TFile("../../Store_data_codes/MCm18_cos.root");
  TTree *treeM1a = (TTree*)fin2a->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1aEvt = treeM1a->GetEntries();
  
  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m, data_y;
  Double_t mc_th, mc_pt, mc_lt, mc_m, mc_y;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  treeM1a->SetBranchAddress("theta", &mc_th);
  treeM1a->SetBranchAddress("dimPt", &mc_pt);
  treeM1a->SetBranchAddress("Rap", &mc_y);
  treeM1a->SetBranchAddress("Mass", &mc_m);
  treeM1a->SetBranchAddress("lt", &mc_lt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);

      // fill the pT histo
      if(data_m > 3.57 && data_m < 3.81 && abs(data_lt) < 0.005 && abs(data_y) < 1.2) {
	h_pT[0]->Fill(data_pt);
      }

      // fill the y histo
      if(data_m > 3.57 && data_m < 3.81 && abs(data_lt) < 0.005) {
	if(data_pt > 20 && data_pt < 100) 
	  h_y[0]->Fill(data_y);
      }
      
      // fill the M histo
      if(abs(data_lt) < 0.005 && abs(data_y) < 1.2) {
	if(data_pt > 20 && data_pt < 100)
	  h_m[0]->Fill(data_m);
      }
      
      // fill the lt histo - in MICRON
      if(data_m > 3.57 && data_m < 3.81 && abs(data_y) < 1.2) {
	if(data_pt > 25 && data_pt < 100) 
	  h_lt[0]->Fill(data_lt*1e4);
      }

      // fill the costh histo - PRSR
      if(data_m > 3.57 && data_m < 3.81 && abs(data_y) < 1.2 && abs(data_lt) < 0.005) {
	if(data_pt > 25 && data_pt < 120) 
	  h_cos[0]->Fill(abs(cos(data_th)));
      }
      // fill the costh histo - NPSR
      if(data_m > 3.57 && data_m < 3.81 && abs(data_y) < 1.2 && data_lt > 0.01 && data_lt < 0.05) {
	if(data_pt > 25 && data_pt < 100) 
	  h_cos[1]->Fill(abs(cos(data_th)));
      }

    }

  cout << "data filled" << endl;
  
  for(int i = 0; i < m1aEvt; i++)
    {
      treeM1a->GetEntry(i);

      // fill the pT histo
      if(mc_m > 3.57 && mc_m < 3.81 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	h_pT[1]->Fill(mc_pt);
      }
      
      // fill the y histo
      if(mc_m > 3.57 && mc_m < 3.81 && abs(mc_lt) < 0.005) {
	if(mc_pt > 20 && mc_pt < 100) 
	  h_y[1]->Fill(mc_y);
      }
      
      // fill the m histo
      if(abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	if(mc_pt > 20 && mc_pt < 100) 
	  h_m[1]->Fill(mc_m);
      }

      // fill the costh histo
      if(mc_m > 3.57 && mc_m < 3.81 && abs(mc_y) < 1.2 && abs(mc_lt) < 0.005) {
	if(mc_pt > 20 && mc_pt < 100) {
	  h_cos[2]->Fill(abs(cos(mc_th)));
	}
      }

    }

  cout << "MC filled" << endl;

  fin->Close();
  fin2a->Close();

  TFile *fout = new TFile("files/store_ANdists.root", "recreate");

  // store the pT dists
  string lbl_pt[] = {"Data", "MC"};
  for(int i = 0; i < 2; i++) {
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
  string lbl_lt[] = {"full"};
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
