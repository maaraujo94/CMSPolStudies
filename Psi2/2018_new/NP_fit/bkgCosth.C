// code to get the 2d sideband/mc ratio hist

void bkgCosth()
{
  // open files and read TTrees
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/data18_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/MC18_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin4 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/MCvh18_cos.root");
  TTree *treeM3 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m3Evt = treeM3->GetEntries();

  // prepare binning and histograms for plots
  const int nPtBins = 7;
  double ptBins[nPtBins+1];
  for(int i=0; i<3; i++) ptBins[i] = 7.*i+25.;
  for(int i=0; i<4; i++) ptBins[i+3] = 46.+10.*i;
  ptBins[7] = 120;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // define aux vals for mass limits
  double m_min[] = {3.4, 3.57, 3.82};
  double m_max[] = {3.52, 3.81, 4.0};

  // data needs LSB/RSB NP (2)
  // MC only needs SR PR (1)
  string lbl[2] = {"LSB", "RSB"};
  TH2D **dataHist = new TH2D*[2];
  TH2D **dataHist_ab = new TH2D*[2];
  for(int i = 0; i < 2; i++) {
    dataHist[i] = new TH2D(Form("dataH%d", i), Form("2018 Data (%s)", lbl[i].c_str()), 40, -1., 1., nPtBins, ptBins);
    dataHist_ab[i] = new TH2D(Form("dataH%d_ab", i), Form("2018 Data (%s)", lbl[i].c_str()), 20, 0, 1., nPtBins, ptBins);
  }
  
  TH2D *mcHist = new TH2D("mcH", "2018 MC", 40, -1., 1., nPtBins, ptBins);
  TH2D *mcHist_ab = new TH2D("mcH_ab", "2018 MC", 20, 0, 1., nPtBins, ptBins);
  
  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m, data_y;
  Double_t mc_th, mc_pt, mc_lt, mc_m, mc_y;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("lt", &data_lt);
  
  treeM1->SetBranchAddress("theta", &mc_th);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("lt", &mc_lt);
  
  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_y) < 1.2 && data_lt > 0.01 && data_lt < 0.05) {
	// LSB
	if(data_m < m_max[0] && data_m > m_min[0]) {
	  dataHist[0]->Fill(cos(data_th), data_pt);
	  dataHist_ab[0]->Fill(abs(cos(data_th)), data_pt);
	}
	// RSB
	else if(data_m < m_max[2] && data_m > m_min[2]) {
	  dataHist[1]->Fill(cos(data_th), data_pt);
	  dataHist_ab[1]->Fill(abs(cos(data_th)), data_pt);
	}
      }
    }
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if(mc_pt > ptBins[0] && mc_pt < 46 && abs(mc_lt) < 0.005 && mc_m < m_max[1] && mc_m > m_min[1] && abs(mc_y) < 1.2) {
	mcHist->Fill(cos(mc_th), mc_pt);
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
      }
    }

  treeM3->SetBranchAddress("theta", &mc_th);
  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("Rap", &mc_y);
  treeM3->SetBranchAddress("lt", &mc_lt);

  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);
      if(mc_pt > 46 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && mc_m < m_max[1] && mc_m > m_min[1] && abs(mc_y) < 1.2) {
	mcHist->Fill(cos(mc_th), mc_pt);
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
      }
    }
  
  // store all costh histograms ( no need to plot )
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.11);
  c->SetLogz();

  for(int i = 0; i < 2; i++) {
    dataHist_ab[i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    dataHist_ab[i]->GetYaxis()->SetTitle("p_{T} (GeV)");
  }
  mcHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  mcHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");
    
  // get the ratio histogram for each bin
  TH2D **ratioHist = new TH2D*[2];
  TH2D **ratioHist_ab = new TH2D*[2];

  TFile *outfile = new TFile("files/bkgHist.root", "recreate");
  
  for(int i = 0; i < 2; i++) {
    ratioHist_ab[i] = new TH2D(Form("ratioH%d_ab", i), Form("2018 Data/MC (%s)", lbl[i].c_str()), 20, 0, 1., nPtBins, ptBins);
    
    ratioHist_ab[i] = (TH2D*)dataHist_ab[i]->Clone(Form("ratioH%d_ab", i));
    ratioHist_ab[i]->Sumw2();
    ratioHist_ab[i]->Divide(mcHist_ab);
    ratioHist_ab[i]->SetTitle(Form("2018 Data/MC (%s)", lbl[i].c_str()));

    dataHist_ab[i]->Write();
    ratioHist_ab[i]->Write();
  }

  outfile->Close();
  int nBinsX = dataHist_ab[0]->GetNbinsX();

}
