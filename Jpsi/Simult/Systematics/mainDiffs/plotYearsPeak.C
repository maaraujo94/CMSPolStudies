// code to compare the deviations for 2017 vs 2018
// using Peak lth bc cuts should work the same way

void plotYearsPeak()
{
  // PART 1 - regular fine-binned results for all studies
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  rHist->SetDirectory(0);
  fIn->Close();
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();

  // get the fit results - 8 sets
  TGraphErrors **graph_lth = new TGraphErrors*[2];
  // 1 - get 2017 and 2018 results
  TFile *fInd17 = new TFile("../../../2017/PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fInd17->Get(Form("graph_lambda_Data"));
  fInd17->Close();
  TFile *fInd18 = new TFile("../../../2018/PR_fit/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fInd18->Get(Form("graph_lambda_Data"));
  fInd18->Close();

  // get the differences
  double diffY[nBinspT], errY[nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    // 2017 - 2018 case can be calculated w their uncertainties
    double val7 = graph_lth[0]->GetY()[i];
    double val8 = graph_lth[1]->GetY()[i];
    diffY[i] = (val7 - val8);
    double err7 = graph_lth[0]->GetEY()[i];
    double err8 = graph_lth[1]->GetEY()[i];
    errY[i] = sqrt(pow(err7, 2) + pow(err8, 2));
  }

  TGraphErrors *g_lthY = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diffY, graph_lth[0]->GetEX(), errY);
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  double da_lim = 0.3;
  
  // draw 2017-2018 with the calculated uncertainty
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#Delta#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("Peak #Delta#lambda_{#theta} (2017-2018)");
  
  g_lthY->SetLineColor(kBlack);
  g_lthY->SetMarkerColor(kBlack);
  g_lthY->SetMarkerStyle(20);
  g_lthY->SetMarkerSize(.75);
  g_lthY->Draw("p same");

  TF1 *f1 = new TF1("f1", "[0]", 25, 120);
  f1->SetParameter(0, 0.02);
  g_lthY->Fit(f1, "0");

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  
  c->SaveAs("plots/lthPeak_Y.pdf");
  c->Clear();
  
  c->Destructor();

}
