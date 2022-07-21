// plot the pT dists of peak, NP and SB data
void plotYields()
{
  const int nPtBins = 7;
  double ptBins[nPtBins+1];
  for(int i=0; i<3; i++) ptBins[i] = 7.*i+25.;
  for(int i=0; i<4; i++) ptBins[i+3] = 46.+10.*i;
  ptBins[7] = 120;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  double m_min[] = {3.4, 3.57, 3.82};
  double m_max[] = {3.52, 3.81, 4.0};

  // define the pT histos
  TH1D *h_tot = new TH1D("totH", "2017 Data", nPtBins, ptBins);
  TH1D *h_NP  = new TH1D("npH",  "2017 Data (NP)",  nPtBins, ptBins);
  TH1D *h_LSB = new TH1D("lsbH", "2017 Data (LSB)", nPtBins, ptBins);
  TH1D *h_RSB = new TH1D("rsbH", "2017 Data (RSB)", nPtBins, ptBins);
  
  // open file, get data
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/data17_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");

  int dEvt = treeD->GetEntries();
  
  // definitions to store data events
  Double_t data_pt, data_lt, data_m;
  
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  // cycle over data, fill the histograms
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      // PR SR
      if(abs(data_lt) < 0.005  && data_m < m_max[1] && data_m > m_min[1]) {
	h_tot->Fill(data_pt);
      }
      // NP SR
      else if(data_lt > 0.01 && data_lt < 0.05 && data_m < m_max[1] && data_m > m_min[1]) {
	h_NP->Fill(data_pt);
      }
      // LSB
      else if(abs(data_lt) < 0.005 && data_m < m_max[0] && data_m > m_min[0]) {
	h_LSB->Fill(data_pt);
      }
      // RSB
      else if(abs(data_lt) < 0.005 && data_m < m_max[2] && data_m > m_min[2]) {
	h_RSB->Fill(data_pt);
      }
    }

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);
  c->SetLogy();

  h_tot->SetLineColor(kBlack);
  h_tot->SetMarkerColor(kBlack);
  h_tot->SetStats(0);
  h_tot->SetMinimum(1e2);
  h_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  h_tot->Draw("error");
  h_NP->SetLineColor(kRed);
  h_NP->SetMarkerColor(kRed);
  h_NP->Draw("error same");
  h_LSB->SetLineColor(kBlue);
  h_LSB->SetMarkerColor(kBlue);
  h_LSB->Draw("error same");
  h_RSB->SetLineColor(kGreen);
  h_RSB->SetMarkerColor(kGreen);
  h_RSB->Draw("error same");

  TLegend *leg = new TLegend(0.7, 0.65, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_tot, "Total", "l");
  leg->AddEntry(h_NP, "NP", "l");
  leg->AddEntry(h_LSB, "LSB", "l");
  leg->AddEntry(h_RSB, "RSB", "l");
  leg->Draw();
  
  c->SaveAs("plots/ptcomp.pdf");
  c->Clear();
  c->SetLogy(0);
  
  TH1F *fr1 = c->DrawFrame(20, 0., 125, 1.);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("Yield ratio (a.u.)");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("Yield (NP) / Yield (Peak)");

  TH1D* f_NP = (TH1D*)h_NP->Clone("f_NP");
  f_NP->Sumw2();
  f_NP->Divide(h_tot);
  f_NP->SetStats(0);
  f_NP->SetLineColor(kRed);
  f_NP->SetMarkerColor(kRed);
  f_NP->Draw("error same");

  c->SaveAs("plots/ptcomp_fNP.pdf");
  c->Clear();

  TH1F *fr2 = c->DrawFrame(20, 0., 125, 1.);
  fr2->SetXTitle("p_{T} (GeV)");
  fr2->SetYTitle("Yield ratio (a.u.)");
  fr2->GetYaxis()->SetTitleOffset(1.5);
  fr2->GetYaxis()->SetLabelOffset(0.01);
  fr2->SetTitle("Yield (bkg) / Yield (Peak)");

  TH1D *f_SB = new TH1D("f_SB", "2017 Data f_{SB}", nPtBins, ptBins);
  f_SB->Sumw2();
  f_SB->Add(h_LSB, h_RSB, 1, 1);
  f_SB->Divide(h_tot);
  f_SB->SetStats(0);
  f_SB->SetLineColor(kBlue);
  f_SB->SetMarkerColor(kBlue);
  f_SB->Draw("error same");
  
  c->SaveAs("plots/ptcomp_fSB.pdf");
  c->Clear();
  c->Destructor();
}
