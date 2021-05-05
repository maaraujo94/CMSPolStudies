// code to get the 2d fine-binned data/mc ratio hist
// plots data, mc and ratio
// saves ratio (normal and |costh|), number of entries in each histo bin

void ratioSave()
{
  // open files and read TTrees
  TFile *fin = new TFile("../Store_data_codes/data18_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("../Store_data_codes/MC18_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../Store_data_codes/MCh18_cos.root");
  TTree *treeM2 = (TTree*)fin3->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();
  
  // prepare binning and histograms for plots
  const int nPtBins = 9;
  double ptBins[10] = {25, 28, 32, 36, 40, 46, 54, 66, 80, 100};
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // define the mass windows for the signal and the sidebands
  double Mq = 3.685, sigW = 0.12; // changed psi mass to 3.685 for easier bins
  double bkgW = sigW;
  double m_min[] = {3.4, Mq-sigW, Mq+bkgW};
  double m_max[] = {Mq-bkgW, Mq+sigW, 4.0};
  string wname[3] = {"L", "S", "R"};
  string wtit[3] = {"LSB", "Peak", "RSB"};

  TH2D **dataHist = new TH2D*[3];
  TH2D **mcHist = new TH2D*[3];
  TH2D **dataHist_ab = new TH2D*[3];
  TH2D **mcHist_ab = new TH2D*[3];

  for(int i_w = 0; i_w < 3; i_w++) {
    dataHist[i_w] = new TH2D(Form("dataH_%s", wname[i_w].c_str()), Form("Data (PR) %s", wtit[i_w].c_str()), 40, -1., 1., nPtBins, ptBins);
    mcHist[i_w] = new TH2D(Form("mcH_%s", wname[i_w].c_str()), Form("MC %s", wtit[i_w].c_str()), 40, -1., 1., nPtBins, ptBins);
    
    dataHist_ab[i_w] = new TH2D(Form("dataH_ab_%s", wname[i_w].c_str()), Form("Data (PR) %s", wtit[i_w].c_str()), 20, 0, 1., nPtBins, ptBins);
    mcHist_ab[i_w] = new TH2D(Form("mcH_ab_%s", wname[i_w].c_str()), Form("MC %s", wtit[i_w].c_str()), 20, 0, 1., nPtBins, ptBins);
  }
  
  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lts, data_m, data_y;
  Double_t mc_th, mc_pt, mc_lts, mc_m, mc_y;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lts", &data_lts);
  
  treeM1->SetBranchAddress("theta", &mc_th);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("lts", &mc_lts);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_lts) < 2.5) {
	for(int i_w = 0; i_w < 3; i_w++) {
	  if(data_m < m_max[i_w] && data_m > m_min[i_w]) {
	    dataHist[i_w]->Fill(cos(data_th), data_pt);
	    dataHist_ab[i_w]->Fill(abs(cos(data_th)), data_pt);
	  }
	}
      }
    }
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);

      if(mc_pt > ptBins[0] && mc_pt < 46 && abs(mc_lts) < 2.5) {
	for(int i_w = 0; i_w < 3; i_w++) {
	  if(mc_m < m_max[i_w] && mc_m > m_min[i_w]) {
	    mcHist[i_w]->Fill(cos(mc_th), mc_pt);
	    mcHist_ab[i_w]->Fill(abs(cos(mc_th)), mc_pt);
	  }
	}
      }
    }

  treeM2->SetBranchAddress("theta", &mc_th);
  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Rap", &mc_y);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("lts", &mc_lts);

  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);

      if(mc_pt > 46 && mc_pt < ptBins[nPtBins] && abs(mc_lts) < 2.5) {
	for(int i_w = 0; i_w < 3; i_w++) {
	  if(mc_m < m_max[i_w] && mc_m > m_min[i_w]) {
	    mcHist[i_w]->Fill(cos(mc_th), mc_pt);
	    mcHist_ab[i_w]->Fill(abs(cos(mc_th)), mc_pt);
	  }
	}
      }
    }

  // plot all signal costh histograms
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.11);
  c->SetLogz();

  dataHist_ab[1]->SetStats(0);
  dataHist_ab[1]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  dataHist_ab[1]->GetYaxis()->SetTitle("p_{T} (GeV)");
  dataHist_ab[1]->Draw("COLZ");
  c->SaveAs("plots/data_2d_abs.pdf");
  c->Clear();

  mcHist_ab[1]->SetStats(0);
  mcHist_ab[1]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  mcHist_ab[1]->GetYaxis()->SetTitle("p_{T} (GeV)");
  mcHist_ab[1]->Draw("COLZ");
  c->SaveAs("plots/mc_2d_abs.pdf");
  c->Clear();
  
  // get the ratio histogram for each bin
  TH2D **ratioHist = new TH2D*[3];
  TH2D **ratioHist_ab = new TH2D*[3];
  TFile *outfile = new TFile("files/ratioHist.root", "recreate");

  for(int i_w = 0; i_w < 3; i_w++) {

    ratioHist[i_w] = new TH2D(Form("ratioH_%s", wname[i_w].c_str()), Form("ratioH %s", wtit[i_w].c_str()), 40, -1., 1., nPtBins, ptBins);

    c->SetLogz(0);
    ratioHist[i_w] = (TH2D*)dataHist[i_w]->Clone(Form("ratioHist_%s", wname[i_w].c_str()));
    ratioHist[i_w]->Sumw2();
    ratioHist[i_w]->Divide(mcHist[1]);
    
    ratioHist_ab[i_w] = new TH2D(Form("ratioH_ab_%s", wname[i_w].c_str()), Form("ratioH abs %s", wtit[i_w].c_str()), 20, 0, 1., nPtBins, ptBins);
  
    c->SetLogz(0);
    ratioHist_ab[i_w] = (TH2D*)dataHist_ab[i_w]->Clone(Form("ratioHist_ab_%s", wname[i_w].c_str()));
    ratioHist_ab[i_w]->Sumw2();
    ratioHist_ab[i_w]->Divide(mcHist_ab[1]);
    ratioHist_ab[i_w]->SetStats(0);
    ratioHist_ab[i_w]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    ratioHist_ab[i_w]->GetYaxis()->SetTitle("p_{T} (GeV)");
    ratioHist_ab[i_w]->SetTitle(Form("%s/MC", wtit[i_w].c_str()));
    ratioHist_ab[i_w]->Draw("COLZ");
    c->SaveAs(Form("plots/ratio%s_2d_abs.pdf", wname[i_w].c_str()));
    c->Clear();

    dataHist[i_w]->Write();
    dataHist_ab[i_w]->Write();
    ratioHist[i_w]->Write();
    ratioHist_ab[i_w]->Write();
    if(i_w == 1) {
      mcHist[i_w]->Write();
      mcHist_ab[i_w]->Write();
    }
  }
  outfile->Close();

  cout << dataHist[1]->GetEntries() << " data events and " << mcHist[1]->GetEntries() << " MC events" << endl;

}
