// plot the pT dists of NP peak and SB data
void plotYields_NP()
{
  const int nPtBins = 19;
  double ptBins[nPtBins+1];
  for(int i = 0; i < 10; i++) ptBins[i] = 25 + 2.5*i;
  for(int i = 0; i < 6; i++) ptBins[i+10] = 50 + 5.*i;
  for(int i = 0; i < 2; i++) ptBins[i+16] = 80 + 10.*i;
  for(int i = 0; i < 2; i++) ptBins[i+18] = 100 + 20.*i;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // define the pT histos
  TH1D *h_NP  = new TH1D("npH",  "Full NP",  nPtBins, ptBins);
  TH1D *h_LSB = new TH1D("lsbH", "Full NP (LSB)", nPtBins, ptBins);
  TH1D *h_RSB = new TH1D("rsbH", "Full NP (RSB)", nPtBins, ptBins);
  
  // open file, get data
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/dataS_cos.root");
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
      // NP SR
      if(data_lt > 0.01 && data_lt < 0.05 && data_m < 3.2 && data_m > 3.0) {
	h_NP->Fill(data_pt);
      }
      // LSB
      else if(data_lt > 0.01 && data_lt < 0.05 && data_m < 2.95 && data_m > 2.92) {
	h_LSB->Fill(data_pt);
      }
      // RSB
      else if(data_lt > 0.01 && data_lt < 0.05 && data_m < 3.28 && data_m > 3.21) {
	h_RSB->Fill(data_pt);
      }
    }
  
  //compare with fit results
  /*  TFile *fin_sb = new TFile("../bkgFits/files/mfit_NP2.root");
  TGraphErrors *fit_fsb = (TGraphErrors*)fin_sb->Get("fit_fBG");
  fin_sb->Close();
  TFile *fin_sbf = new TFile("../bkgFits/files/mfit_NP2free.root");
  TGraphErrors *fit_fsbf = (TGraphErrors*)fin_sbf->Get("fit_fBG");
  fin_sbf->Close();*/
  TFile *fin_sbc = new TFile("../bkgFits/files/bkgFrac_NP.root");
  TH1D *fit_fsbc = (TH1D*)fin_sbc->Get("fbkg_unc");
  fit_fsbc->SetDirectory(0);
  fin_sbc->Close();

  // now the plotting
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);
  c->SetLogy();

  h_NP->SetLineColor(kRed);
  h_NP->SetMarkerColor(kRed);
  h_NP->SetStats(0);
  h_NP->SetMinimum(1e2);
  h_NP->GetXaxis()->SetTitle("p_{T} (GeV)");
  h_NP->Draw("error");
  h_LSB->SetLineColor(kBlue);
  h_LSB->SetMarkerColor(kBlue);
  h_LSB->Draw("error same");
  h_RSB->SetLineColor(kGreen);
  h_RSB->SetMarkerColor(kGreen);
  h_RSB->Draw("error same");

  TLegend *leg = new TLegend(0.7, 0.65, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_NP, "NP", "l");
  leg->AddEntry(h_LSB, "LSB", "l");
  leg->AddEntry(h_RSB, "RSB", "l");
  leg->Draw();
  
  c->SaveAs("plots/ptcomp_NP.pdf");
  c->Clear();
  c->SetLogy(0);

  // compare f_bkg
  TH1F *fr2 = c->DrawFrame(20, 0., 125, 0.1);
  fr2->SetXTitle("p_{T} (GeV)");
  fr2->SetYTitle("Yield ratio (a.u.)");
  fr2->GetYaxis()->SetTitleOffset(1.5);
  fr2->GetYaxis()->SetLabelOffset(0.01);
  fr2->SetTitle("Yield (bkg) / Yield (Peak)");

  TH1D *f_SB = new TH1D("f_SB", "Full NP f_{SB}", nPtBins, ptBins);
  f_SB->Sumw2();
  f_SB->Add(h_LSB, h_RSB, 1, 1);
  f_SB->Divide(h_NP);
  f_SB->SetStats(0);
  f_SB->SetLineColor(kBlue);
  f_SB->SetMarkerColor(kBlue);
  f_SB->Draw("error same");
  
  c->SaveAs("plots/ptcomp_NP_fSB.pdf");
  c->Clear();

  // same but scaled and comparing to f_bkg from fit
  TH1F *fr2c = c->DrawFrame(20, 0., 125, 0.02);
  fr2c->SetXTitle("p_{T} (GeV)");
  fr2c->SetYTitle("Yield ratio (a.u.)");
  fr2c->GetYaxis()->SetTitleOffset(1.3);
  fr2c->GetYaxis()->SetLabelOffset(0.01);
  fr2c->SetTitle("Yield (bkg) / Yield (Peak)");

  f_SB->Scale(1./f_SB->Integral("width"));
  f_SB->SetLineColor(kBlue);
  f_SB->SetMarkerColor(kBlue);
  f_SB->Draw("error same");

  fit_fsbc->Scale(1./fit_fsbc->Integral("width"));
  fit_fsbc->SetLineColor(kBlack);
  fit_fsbc->SetMarkerColor(kBlack);
  fit_fsbc->Draw("error same");
  
  c->SaveAs("plots/fitcomp_NP_fSB.pdf");
  c->Clear();

  TFile *fout = new TFile("files/yieldComp_NP.root", "recreate");
  f_SB->Write();
  fit_fsbc->Write();
  fout->Close();
  
  c->Destructor();
}
