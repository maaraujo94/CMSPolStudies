// code to compare the deviations for several scenarios
// 1) 0.8*fBG- Run2
// 2) 1.2*fBG - Run2

void plot_SB()
{
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - 3 sets
  TGraphErrors **graph_lth = new TGraphErrors*[3];
  TGraphErrors **graph_lth_NP = new TGraphErrors*[3];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  graph_lth_NP[0] = (TGraphErrors*)fIndB->Get("graph_lambda_NP");
  fIndB->Close();
  // 2 - get lambda values for full LSB
  TFile *fIndL = new TFile("../LSB/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fIndL->Get(Form("graph_lambda_J"));
  graph_lth_NP[1] = (TGraphErrors*)fIndL->Get(Form("graph_lambda_NP"));
  fIndL->Close();
  // 3 - get lambda values for full RSB
  TFile *fIndH = new TFile("../RSB/files/finalFitRes.root");
  graph_lth[2] = (TGraphErrors*)fIndH->Get(Form("graph_lambda_J"));
  graph_lth_NP[2] = (TGraphErrors*)fIndH->Get(Form("graph_lambda_NP"));
  fIndH->Close();

  // get the differences
  double diff[2][nBinspT], za[nBinspT], err[2][nBinspT];
  double diff_NP[2][nBinspT], err_NP[2][nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    diff[0][i] = (graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[1][i] = (graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i]);
 
    double unc1 = graph_lth[0]->GetEY()[i];
    double unc2 = graph_lth[1]->GetEY()[i];
    err[0][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lth[2]->GetEY()[i];
    err[1][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    za[i] = 0;

    diff_NP[0][i] = (graph_lth_NP[1]->GetY()[i] - graph_lth_NP[0]->GetY()[i]);
    diff_NP[1][i] = (graph_lth_NP[2]->GetY()[i] - graph_lth_NP[0]->GetY()[i]);
 
    unc1 = graph_lth_NP[0]->GetEY()[i];
    unc2 = graph_lth_NP[1]->GetEY()[i];
    err_NP[0][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lth_NP[2]->GetEY()[i];
    err_NP[1][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));
  }
  TGraphErrors *g_lthL = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[0], graph_lth[0]->GetEX(), err[0]);//za);
  TGraphErrors *g_lthH = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[1], graph_lth[0]->GetEX(), err[1]);//za);

  TGraphErrors *g_lthL_NP = new TGraphErrors(nBinspT, graph_lth_NP[0]->GetX(), diff_NP[0], graph_lth_NP[0]->GetEX(), err_NP[0]);// za);
  TGraphErrors *g_lthH_NP = new TGraphErrors(nBinspT, graph_lth_NP[0]->GetX(), diff_NP[1], graph_lth_NP[0]->GetEX(), err_NP[1]);// za);

  TGraphErrors *g_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_lth[0]->GetEY());
    TGraphErrors *g_unc_NP = new TGraphErrors(nBinspT, graph_lth_NP[0]->GetX(), za, graph_lth_NP[0]->GetEX(), graph_lth_NP[0]->GetEY());

  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  
  double da_lim = 0.2;
  
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#Delta#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("prompt #Delta#lambda_{#theta} (SB)");
  
  g_lthL->SetLineColor(kBlue);
  g_lthL->SetMarkerColor(kBlue);
  g_lthL->SetMarkerStyle(20);
  g_lthL->SetMarkerSize(.75);
  g_lthL->Draw("p same");

  g_lthH->SetLineColor(kRed);
  g_lthH->SetMarkerColor(kRed);
  g_lthH->SetMarkerStyle(20);
  g_lthH->SetMarkerSize(.75);
  g_lthH->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->Draw("ce3");

  TLegend *leg = new TLegend(0.77, 0.7, 0.97, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(g_lthL, "LSB", "pl");
  leg->AddEntry(g_lthH, "RSB", "pl");
  leg->Draw();
  
  c->SaveAs("lth_absDiff_SB.pdf");
  c->Clear();

  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#Delta#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("non-prompt #Delta#lambda_{#theta} (SB)");
  
  g_lthL_NP->SetLineColor(kBlue);
  g_lthL_NP->SetMarkerColor(kBlue);
  g_lthL_NP->SetMarkerStyle(20);
  g_lthL_NP->SetMarkerSize(.75);
  g_lthL_NP->Draw("p same");

  g_lthH_NP->SetLineColor(kRed);
  g_lthH_NP->SetMarkerColor(kRed);
  g_lthH_NP->SetMarkerStyle(20);
  g_lthH_NP->SetMarkerSize(.75);
  g_lthH_NP->Draw("p same");

  g_unc_NP->SetLineColor(kBlack);
  g_unc_NP->SetFillColorAlpha(kBlack, 0.1);
  g_unc_NP->Draw("ce3");

  leg->Draw();
  
  c->SaveAs("lthNP_absDiff_SB.pdf");
  c->Clear();
  c->Destructor();

  fIn->Close();


}
