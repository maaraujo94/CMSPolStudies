// code to compare the deviations for changed fNP

void plot_fNP()
{
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - 2 sets
  TGraphErrors **graph_lth = new TGraphErrors*[2];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  fIndB->Close();
  // 2 - get lambda values for changed fNP
  TFile *fIndL = new TFile("../fnp_est/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fIndL->Get(Form("graph_lambda_J"));
  fIndL->Close();

  // get the differences
  double diff[1][nBinspT], za[nBinspT], err[1][nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    diff[0][i] = (graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i]);
   
    double unc1 = graph_lth[0]->GetEY()[i];
    double unc2 = graph_lth[1]->GetEY()[i];
    err[0][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    za[i] = 0;
  }
  TGraphErrors *g_lthL = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[0], graph_lth[0]->GetEX(), err[0]);//za);

  TGraphErrors *g_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_lth[0]->GetEY());
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  
  double da_lim = 0.2;
  
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#Delta#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("prompt #Delta#lambda_{#theta} (f_{bkg})");
  
  g_lthL->SetLineColor(kBlue);
  g_lthL->SetMarkerColor(kBlue);
  g_lthL->SetMarkerStyle(20);
  g_lthL->SetMarkerSize(.75);
  g_lthL->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->Draw("ce3");

  TLegend *leg = new TLegend(0.77, 0.7, 0.97, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(g_lthL, "est f_{NP}", "pl");
  //leg->Draw();

  c->SaveAs("lth_absDiff_fNP.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();


}
