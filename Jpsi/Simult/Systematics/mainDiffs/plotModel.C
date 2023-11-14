// code to compare the deviations for several scenarios
// 1) linear mass bkg - Run2
// 2) SB-only mass bkg - Run 2

void plotModel()
{
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - 3 sets
  TGraphErrors **graph_lth = new TGraphErrors*[3];
  TGraphErrors **graph_lthNP = new TGraphErrors*[2];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  graph_lthNP[0] = (TGraphErrors*)fIndB->Get("graph_lambda_NP"); 
  fIndB->Close();
  // 1 - get results for SB-only fit
  TFile *fIndSB = new TFile("../SB_fit/PR_fit/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fIndSB->Get("graph_lambda_J");
  graph_lthNP[1] = (TGraphErrors*)fIndSB->Get("graph_lambda_NP"); 
  fIndSB->Close();
  // 2 - get lambda values for linear mass bkg
  TFile *fIndmb = new TFile("../Mass_bkg_alt/files/finalFitRes.root");
  graph_lth[2] = (TGraphErrors*)fIndmb->Get(Form("graph_lambda_J"));
  fIndmb->Close();

  // get the differences
  double diff[2][nBinspT], diffNP[1][nBinspT], za[nBinspT];
  double err[2][nBinspT], errNP[1][nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    diff[0][i] = (graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[1][i] = (graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i]);
 
    double unc1 = graph_lth[0]->GetEY()[i];
    double unc2 = graph_lth[1]->GetEY()[i];
    err[0][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lth[2]->GetEY()[i];
    err[1][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    diffNP[0][i] = (graph_lthNP[1]->GetY()[i] - graph_lthNP[0]->GetY()[i]);
 
    unc1 = graph_lthNP[0]->GetEY()[i];
    unc2 = graph_lthNP[1]->GetEY()[i];
    errNP[0][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    za[i] = 0;
  }
  TGraphErrors *g_lthSB = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[0], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthNPSB = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), diffNP[0], graph_lthNP[0]->GetEX(), za);
  TGraphErrors *g_lthMB = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[1], graph_lth[0]->GetEX(), za);
  
  TGraphErrors *g_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_lth[0]->GetEY());
  TGraphErrors *g_uncNP = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), za, graph_lthNP[0]->GetEX(), graph_lthNP[0]->GetEY());
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  
  double d_lim = 60;

  // FIRST - draw the lin mass bkg check
  double da_lim = 0.3;
  
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#Delta#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("prompt #Delta#lambda_{#theta} (linear mass bkg)");
  
  g_lthMB->SetLineColor(kRed);
  g_lthMB->SetMarkerColor(kRed);
  g_lthMB->SetMarkerStyle(20);
  g_lthMB->SetMarkerSize(.75);
  g_lthMB->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->Draw("ce3");

  c->SaveAs("plots/lth_absDiff_M.pdf");
  c->Clear();

  // SECOND - draw the SB-only fit lth (PR)
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#Delta#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("prompt #Delta#lambda_{#theta} (SB-only mass fit)");
  
  g_lthSB->SetLineColor(kRed);
  g_lthSB->SetMarkerColor(kRed);
  g_lthSB->SetMarkerStyle(20);
  g_lthSB->SetMarkerSize(.75);
  g_lthSB->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->Draw("ce3");

  c->SaveAs("plots/lth_absDiff_SB.pdf");
  c->Clear();
  
  // THIRD - draw the SB-only fit lth (NP)
  TH1F *fl3 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl3->SetXTitle("p_{T} (GeV)");
  fl3->SetYTitle("#Delta#lambda_{#theta}");
  fl3->GetYaxis()->SetTitleOffset(1.3);
  fl3->GetYaxis()->SetLabelOffset(0.01);
  fl3->SetTitle("non-prompt #Delta#lambda_{#theta} (SB-only mass fit)");
  
  g_lthNPSB->SetLineColor(kRed);
  g_lthNPSB->SetMarkerColor(kRed);
  g_lthNPSB->SetMarkerStyle(20);
  g_lthNPSB->SetMarkerSize(.75);
  g_lthNPSB->Draw("p same");

  g_uncNP->SetLineColor(kBlack);
  g_uncNP->SetFillColorAlpha(kBlack, 0.1);
  g_uncNP->Draw("ce3");

  c->SaveAs("plots/lthNP_absDiff_SB.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();


}
