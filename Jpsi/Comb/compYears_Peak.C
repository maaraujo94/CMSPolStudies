// code to plot the fit results

void compYears_Peak()
{
  // get the histo limits
  TFile *fIn = new TFile("../2017/PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get lambda values for each bin
  const int n_lbl = 1;
  string lbl[] = {"Data"};
  TFile *fInd7 = new TFile("../2017/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth7 = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lth7[i_t] = (TGraphErrors*)fInd7->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd7->Close();
  TFile *fInd8 = new TFile("../2018/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth8 = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lth8[i_t] = (TGraphErrors*)fInd8->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd8->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("peak J/#psi #lambda_{#theta}");

  graph_lth7[0]->SetLineColor(kBlue);
  graph_lth7[0]->SetMarkerSize(.25);
  graph_lth7[0]->SetMarkerStyle(20);
  graph_lth7[0]->SetMarkerColor(kBlue);
  graph_lth7[0]->SetMarkerColor(kBlue);
  graph_lth7[0]->Draw("p same");
  graph_lth8[0]->SetLineColor(kBlack);
  graph_lth8[0]->SetMarkerStyle(20);
  graph_lth8[0]->SetMarkerSize(.25);
  graph_lth8[0]->SetMarkerColor(kBlack);
  graph_lth8[0]->Draw("p same");
 
  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg = new TLegend(0.75, 0.75, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(graph_lth7[0], "2017", "pl");
  leg->AddEntry(graph_lth8[0], "2018", "pl");
  leg->Draw();
  
  c->SaveAs("par_lth_years_Peak.pdf");
  c->Destructor();
  
  fIn->Close();
}
