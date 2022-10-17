#import "../rcut.C"

// code to get the 2d sideband/mc ratio hist

void bkgCosth()
{
  // prepare binning for plots
  TH2D *h_ref = new TH2D();  
  TFile *fref = new TFile("files/histoStore.root");
  fref->GetObject("ratioH_ab", h_ref);
  h_ref->SetDirectory(0);
  fref->Close();

  const int nBinsX = h_ref->GetNbinsX();
  double minX = h_ref->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_ref->GetXaxis()->GetBinUpEdge(nBinsX);
  const int nPtBins = h_ref->GetNbinsY();
  const double *ptBins = h_ref->GetYaxis()->GetXbins()->GetArray();

  // open files and read TTrees
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/dataS_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCS_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MChS_cos.root");
  TTree *treeM2 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCvhS_cos.root");
  TTree *treeM3 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries();

  // data needs LSB/RSB PR (2)
  // MC only needs SR PR (1)
  string lbl[2] = {"LSB", "RSB"};
  TH2D **dataHist_ab = new TH2D*[2];
  for(int i = 0; i < 2; i++) {
    dataHist_ab[i] = new TH2D(Form("dataH%d_ab", i), Form("Full Data (%s)", lbl[i].c_str()), nBinsX, minX, maxX, nPtBins, ptBins);
  }
  
  TH2D *mcHist_ab = new TH2D("mcH_ab", "Full MC", nBinsX, minX, maxX, nPtBins, ptBins);

  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m, data_y;
  Double_t dR;
  Double_t mc_th, mc_pt, mc_lt, mc_m, mc_y;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("lt", &data_lt);
  treeD->SetBranchAddress("DeltaR", &dR);
  
  treeM1->SetBranchAddress("theta", &mc_th);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("lt", &mc_lt);
  treeM1->SetBranchAddress("DeltaR", &dR);
  
  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(dR > r_cut)
	if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_y) < 1.2 && abs(data_lt) < 0.005) {
	  // LSB
	  if(data_m < 2.95 && data_m > 2.92) {
	    dataHist_ab[0]->Fill(abs(cos(data_th)), data_pt);
	  }
	  // RSB
	  else if(data_m < 3.28 && data_m > 3.21) {
	    dataHist_ab[1]->Fill(abs(cos(data_th)), data_pt);
	  }
	}
    }
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if(dR > r_cut)
	if(mc_pt > ptBins[0] && mc_pt < 47.5 && abs(mc_lt) < 0.005 && mc_m < 3.2 && mc_m > 3.0 && abs(mc_y) < 1.2) {
	  mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
	}
    }

  treeM2->SetBranchAddress("theta", &mc_th);
  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("Rap", &mc_y);
  treeM2->SetBranchAddress("lt", &mc_lt);
  treeM2->SetBranchAddress("DeltaR", &dR);

  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
      if(dR > r_cut)
	if(mc_pt > 47.5 && mc_pt < 70 && abs(mc_lt) < 0.005 && mc_m < 3.2 && mc_m > 3.0 && abs(mc_y) < 1.2) {
	  mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
	}
    }
  

  treeM3->SetBranchAddress("theta", &mc_th);
  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("Rap", &mc_y);
  treeM3->SetBranchAddress("lt", &mc_lt);
  treeM3->SetBranchAddress("DeltaR", &dR);

  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);
      if(dR > r_cut)
	if(mc_pt > 70 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && mc_m < 3.2 && mc_m > 3.0 && abs(mc_y) < 1.2) {
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
  TH2D **ratioHist_ab = new TH2D*[2];

  TFile *outfile = new TFile("files/bkgHist.root", "recreate");
  
  for(int i = 0; i < 2; i++) {
    ratioHist_ab[i] = new TH2D(Form("ratioH%d_ab", i), Form("Full Data/MC (%s)", lbl[i].c_str()), nBinsX, minX, maxX, nPtBins, ptBins);
    
    ratioHist_ab[i] = (TH2D*)dataHist_ab[i]->Clone(Form("ratioH%d_ab", i));
    ratioHist_ab[i]->Sumw2();
    ratioHist_ab[i]->Divide(mcHist_ab);
    ratioHist_ab[i]->SetTitle(Form("Full Data/MC (%s)", lbl[i].c_str()));

    dataHist_ab[i]->Write();
    ratioHist_ab[i]->Write();
  }

  outfile->Close();
 
}
