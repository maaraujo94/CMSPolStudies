// plot the pT dists of peak, NP and SB data
void plotYields()
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
  TH1D *h_tot = new TH1D("totH", "Full Data", nPtBins, ptBins);
  TH1D *h_NP  = new TH1D("npH",  "Full Data (NP)",  nPtBins, ptBins);
  TH1D *h_LSB = new TH1D("lsbH", "Full Data (LSB)", nPtBins, ptBins);
  TH1D *h_RSB = new TH1D("rsbH", "Full Data (RSB)", nPtBins, ptBins);
  
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
      // PR SR
      if(abs(data_lt) < 0.005  && data_m < 3.2 && data_m > 3.0) {
	h_tot->Fill(data_pt);
      }
      // NP SR
      else if(data_lt > 0.01 && data_lt < 0.05 && data_m < 3.2 && data_m > 3.0) {
	h_NP->Fill(data_pt);
      }
      // LSB
      else if(abs(data_lt) < 0.005 && data_m < 2.95 && data_m > 2.92) {
	h_LSB->Fill(data_pt);
      }
      // RSB
      else if(abs(data_lt) < 0.005 && data_m < 3.28 && data_m > 3.21) {
	h_RSB->Fill(data_pt);
      }
    }
  
  //compare with fit results
  TFile *fin_sb = new TFile("../bkgFits/files/bkgFrac.root");
  TH1D *fit_fsb = (TH1D*)fin_sb->Get("fbkg_unc");
  fit_fsb->SetDirectory(0);
  fin_sb->Close();
  TFile *fin_np = new TFile("../PR_fit/files/NPFrac.root");
  TH1D *fit_fnp = (TH1D*)fin_np->Get("fnp_unc");
  fit_fnp->SetDirectory(0);
  fin_np->Close();

  // now the plotting
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

  // comparing fNP
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

  // same but scaled and comparing to f_NP from fit
  TH1F *fr1c = c->DrawFrame(20, 0., 125, 0.02);
  fr1c->SetXTitle("p_{T} (GeV)");
  fr1c->SetYTitle("Yield ratio (a.u.)");
  fr1c->GetYaxis()->SetTitleOffset(1.5);
  fr1c->GetYaxis()->SetLabelOffset(0.01);
  fr1c->SetTitle("Yield (NP) / Yield (Peak)");

  f_NP->Scale(1./f_NP->Integral("width"));
  f_NP->SetLineColor(kBlue);
  f_NP->SetMarkerColor(kBlue);
  f_NP->Draw("error same");

  fit_fnp->Scale(1./fit_fnp->Integral("width"));
  fit_fnp->SetLineColor(kBlack);
  fit_fnp->SetMarkerColor(kBlack);
  fit_fnp->Draw("error same");
    
  c->SaveAs("plots/fitcomp_fNP.pdf");
  c->Clear();

  // now f_bkg
  TH1F *fr2 = c->DrawFrame(20, 0., 125, 0.1);
  fr2->SetXTitle("p_{T} (GeV)");
  fr2->SetYTitle("Yield ratio (a.u.)");
  fr2->GetYaxis()->SetTitleOffset(1.5);
  fr2->GetYaxis()->SetLabelOffset(0.01);
  fr2->SetTitle("Yield (bkg) / Yield (Peak)");

  TH1D *f_SB = new TH1D("f_SB", "Full Data f_{SB}", nPtBins, ptBins);
  f_SB->Sumw2();
  f_SB->Add(h_LSB, h_RSB, 1, 1);
  f_SB->Divide(h_tot);
  f_SB->SetStats(0);
  f_SB->SetLineColor(kBlue);
  f_SB->SetMarkerColor(kBlue);
  f_SB->Draw("error same");
  
  c->SaveAs("plots/ptcomp_fSB.pdf");
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

  fit_fsb->Scale(1./fit_fsb->Integral("width"));
  fit_fsb->SetLineColor(kBlack);
  fit_fsb->SetMarkerColor(kBlack);
  fit_fsb->Draw("error same");
    
  c->SaveAs("plots/fitcomp_fSB.pdf");
  c->Clear();

  TFile *fout = new TFile("files/yieldComp.root", "recreate");
  f_NP->Write();
  f_SB->Write();
  fit_fnp->Write();
  fit_fsb->Write();
  fout->Close();
  
  c->Destructor();
}
