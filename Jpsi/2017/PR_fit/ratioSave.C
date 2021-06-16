// code to get the 2d fine-binned data/mc ratio hist
// plots data, mc and ratio
// saves ratio (normal and |costh|), number of entries in each histo bin

// macro to determine sigma_m(y)
double sigma_y(double y, double c1, double c2, double m)
{
  // define fit limits
  double y1 = 0.25, y2 = 0.5, y3 = 0.7;
  
  // convert to easy system
  double m1 = (c2-c1)/(y2-y1);
  double b1 = c1-m1*y1;
  double m2 = m;
  double b2 = c2-m2*y3;

  // return the correct function depending on input y
  if(y < 0 || y > 1.2) {
    cout << "incorrect y input" << endl;
    return 0;
  }
  
  else if(y <= y1) return c1;
  else if(y <= y2) return b1+m1*y;
  else if(y <= y3) return c2;
  else return b2 + m2*y; 
}

void ratioSave()
{
  // start by getting the parameters for the sigma(y; pt) function
  ifstream ifile;
  string data;
  int pt_bins = 7;
  double pt_min[pt_bins], pt_max[pt_bins], c1[pt_bins], c2[pt_bins], m[pt_bins], aux, sig;
  ifile.open("../plots_massFit/fit_sig.txt");
  getline(ifile, data);
  getline(ifile, data);
  for(int i = 0; i < 7; i++) {
    ifile >> pt_min[i] >> pt_max[i] >> c1[i] >> aux >> c2[i] >> aux >> m[i] >> aux;
  }
  ifile.close();
  
  // open files and read TTrees
  TFile *fin = new TFile("../../Store_data_codes/data17_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/MC17_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../../Store_data_codes/MCh17_cos.root");
  TTree *treeM2 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("../../Store_data_codes/MCvh17_cos.root");
  TTree *treeM3 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries();

  // prepare binning and histograms for plots 
  const int nPtBins = 45;
  double ptBins[nPtBins+1];
  for(int i=0; i<15; i++) ptBins[i] = i+25.;
  for(int i=0; i<31; i++) ptBins[i+15] = 40.+2.*i;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;
 
  TH2D *dataHist = new TH2D("dataH", "2017 Data (PR)", 40, -1., 1., nPtBins, ptBins);
  TH2D *mcHist = new TH2D("mcH", "2017 MC", 40, -1., 1., nPtBins, ptBins);

  TH2D *dataHist_ab = new TH2D("dataH_ab", "2017 Data (PR)", 20, 0, 1., nPtBins, ptBins);
  TH2D *mcHist_ab = new TH2D("mcH_ab", "2017 MC", 20, 0, 1., nPtBins, ptBins);
  
  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_lte, data_m, data_y;
  Double_t mc_th, mc_pt, mc_lt, mc_lte, mc_m, mc_y;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  treeD->SetBranchAddress("lterr", &data_lte);
  
  treeM1->SetBranchAddress("theta", &mc_th);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("lt", &mc_lt);
  treeM1->SetBranchAddress("lterr", &mc_lte);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_lt/data_lte) < 2.5) {
	for(int j = 0; j < pt_bins; j++)
	  if(data_pt < pt_max[j] && data_pt > pt_min[j]) {
	    sig = sigma_y(data_y, c1[j], c2[j], m[j]);
	  }

	if(data_m < 3.097 + 2.5*sig && data_m > 3.097 - 2.5*sig) {
	  dataHist->Fill(cos(data_th), data_pt);
	  dataHist_ab->Fill(abs(cos(data_th)), data_pt);
	}
      }
    }
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if(mc_pt > ptBins[0] && mc_pt < 46 && abs(mc_lt/mc_lte) < 2.5) {
	for(int j = 0; j < pt_bins; j++)
	  if(mc_pt < pt_max[j] && mc_pt > pt_min[j])
	    sig = sigma_y(mc_y, c1[j], c2[j], m[j]);
	
	if(mc_m < 3.097 + 2.5*sig && mc_m > 3.097 - 2.5*sig) {
	  mcHist->Fill(cos(mc_th), mc_pt);
	  mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
	}
      }
    }

  treeM2->SetBranchAddress("theta", &mc_th);
  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Rap", &mc_y);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("lt", &mc_lt);
  treeM2->SetBranchAddress("lterr", &mc_lte);

  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
      if(mc_pt > 46 && mc_pt < 66 && abs(mc_lt/mc_lte) < 2.5) {
	for(int j = 0; j < pt_bins; j++)
	  if(mc_pt < pt_max[j] && mc_pt > pt_min[j])
	    sig = sigma_y(mc_y, c1[j], c2[j], m[j]);
	
	if(mc_m < 3.097 + 2.5*sig && mc_m > 3.097 - 2.5*sig) {
	  mcHist->Fill(cos(mc_th), mc_pt);
	  mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
	}
      }
    }

  treeM3->SetBranchAddress("theta", &mc_th);
  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Rap", &mc_y);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("lt", &mc_lt);
  treeM3->SetBranchAddress("lterr", &mc_lte);

  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);
      if(mc_pt > 66 && mc_pt < ptBins[nPtBins] && abs(mc_lt/mc_lte) < 2.5) {
	for(int j = 0; j < pt_bins; j++)
	  if(mc_pt < pt_max[j] && mc_pt > pt_min[j])
	    sig = sigma_y(mc_y, c1[j], c2[j], m[j]);
	
	if(mc_m < 3.097 + 2.5*sig && mc_m > 3.097 - 2.5*sig) {
	  mcHist->Fill(cos(mc_th), mc_pt);
	  mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
	}
      }
    }

  
  // plot all costh histograms
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.11);
  c->SetLogz();
  
  dataHist_ab->SetStats(0);
  dataHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  dataHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");
  dataHist_ab->Draw("COLZ");
  c->SaveAs("plots/data_2d_abs.pdf");
  c->Clear();

  mcHist_ab->SetStats(0);
  mcHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  mcHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");
  mcHist_ab->Draw("COLZ");
  c->SaveAs("plots/mc_2d_abs.pdf");
  c->Clear();
    
  // get the ratio histogram for each bin
  TH2D *ratioHist = new TH2D("ratioH", "ratioH", 40, -1., 1., nPtBins, ptBins);

  c->SetLogz(0);
  ratioHist = (TH2D*)dataHist->Clone(Form("ratioHist"));
  ratioHist->Sumw2();
  ratioHist->Divide(mcHist);

  TH2D *ratioHist_ab = new TH2D("ratioH_ab", "ratioH abs", 20, 0, 1., nPtBins, ptBins);

  c->SetLogz(0);
  ratioHist_ab = (TH2D*)dataHist_ab->Clone(Form("ratioHist_ab"));
  ratioHist_ab->Sumw2();
  ratioHist_ab->Divide(mcHist_ab);
  ratioHist_ab->SetTitle("2017 PR/MC");
  ratioHist_ab->Draw("COLZ");
  c->SaveAs("plots/ratio_2d_abs.pdf");
  c->Clear();

  TFile *outfile = new TFile("files/ratioHist.root", "recreate");
  dataHist_ab->Write();
  mcHist_ab->Write();
  ratioHist->Write();
  ratioHist_ab->Write();
  outfile->Close();

  cout << dataHist->GetEntries() << " data events and " << mcHist->GetEntries() << " MC events" << endl;

  cout << dataHist->Integral(1, 18) << " data events and " << mcHist->Integral(1, 18) << " MC events in pT range 1" << endl;
  cout << dataHist->Integral(19, 28) << " data events and " << mcHist->Integral(19, 28) << " MC events in pT range 2" << endl;
  cout << dataHist->Integral(29, 45) << " data events and " << mcHist->Integral(29, 45) << " MC events in pT range 3" << endl;  

}
