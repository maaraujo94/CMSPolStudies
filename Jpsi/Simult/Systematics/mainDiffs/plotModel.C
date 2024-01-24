// code to compare the deviations for several scenarios
// 1) linear mass bkg - Run2
// 2) genDist uses only one SB
// 3) pT-indept mu2
// 4) fNP_psi estimate variation

void plotModel()
{
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - 6 sets
  TGraphErrors **graph_lth = new TGraphErrors*[6];
  TGraphErrors **graph_lthNP = new TGraphErrors*[5];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  graph_lthNP[0] = (TGraphErrors*)fIndB->Get("graph_lambda_NP"); 
  fIndB->Close();
  // 1 - get lambda values for linear mass bkg - not currently doing, haven't converged mass fit yet

  // 2 - get results for (L/R)SB-only bkg models
  TFile *fIndLSB = new TFile("../LSB/files/finalFitRes.root");
  graph_lth[2] = (TGraphErrors*)fIndLSB->Get("graph_lambda_J");
  graph_lthNP[2] = (TGraphErrors*)fIndLSB->Get("graph_lambda_NP"); 
  fIndLSB->Close();
  TFile *fIndRSB = new TFile("../RSB/files/finalFitRes.root");
  graph_lth[3] = (TGraphErrors*)fIndRSB->Get("graph_lambda_J");
  graph_lthNP[3] = (TGraphErrors*)fIndRSB->Get("graph_lambda_NP"); 
  fIndRSB->Close();
  // 3 - get results for NP mass fit with pT-indep mu2
  TFile *fIndmu2 = new TFile("../mu2_indep/files/finalFitRes.root");
  graph_lth[4] = (TGraphErrors*)fIndmu2->Get("graph_lambda_J");
  graph_lthNP[4] = (TGraphErrors*)fIndmu2->Get("graph_lambda_NP"); 
  fIndmu2->Close();
  // 4 - get results for variation of fNP_psi estimate
  TFile *fIndfnp = new TFile("../fnp_est/files/finalFitRes.root");
  graph_lth[5] = (TGraphErrors*)fIndfnp->Get("graph_lambda_J");
  fIndfnp->Close();
  
  // get the differences
  double diff[5][nBinspT], diffNP[4][nBinspT], za[nBinspT];
  double err[5][nBinspT], errNP[5][nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    diff[1][i] = (graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[2][i] = (graph_lth[3]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[3][i] = (graph_lth[4]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[4][i] = (graph_lth[5]->GetY()[i] - graph_lth[0]->GetY()[i]);
 
    double unc1 = graph_lth[0]->GetEY()[i];
    double unc2 = graph_lth[2]->GetEY()[i];
    err[1][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lth[3]->GetEY()[i];
    err[2][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lth[4]->GetEY()[i];
    err[3][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lth[5]->GetEY()[i];
    err[4][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    diffNP[1][i] = (graph_lthNP[2]->GetY()[i] - graph_lthNP[0]->GetY()[i]);
    diffNP[2][i] = (graph_lthNP[3]->GetY()[i] - graph_lthNP[0]->GetY()[i]);
    diffNP[3][i] = (graph_lthNP[4]->GetY()[i] - graph_lthNP[0]->GetY()[i]);
 
    unc1 = graph_lthNP[0]->GetEY()[i];
    unc2 = graph_lthNP[2]->GetEY()[i];
    errNP[1][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lthNP[3]->GetEY()[i];
    errNP[2][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lthNP[4]->GetEY()[i];
    errNP[3][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    za[i] = 0;
  }
  TGraphErrors *g_lthLSB = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[1], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthRSB = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[2], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthmu2 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[3], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthfnp = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[4], graph_lth[0]->GetEX(), za);
  
  TGraphErrors *g_lthNPLSB = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), diffNP[1], graph_lthNP[0]->GetEX(), za);
  TGraphErrors *g_lthNPRSB = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), diffNP[2], graph_lthNP[0]->GetEX(), za);
  TGraphErrors *g_lthNPmu2 = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), diffNP[3], graph_lthNP[0]->GetEX(), za);
  
  TGraphErrors *g_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_lth[0]->GetEY());
  TGraphErrors *g_uncNP = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), za, graph_lthNP[0]->GetEX(), graph_lthNP[0]->GetEY());
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  
  double d_lim = 60;

  // FIRST - draw the lin mass bkg check
  double da_lim = 0.3;
  
  /*  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
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
  c->Clear();*/

  // SECOND - draw the L/RSB-only gen Dist
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#Delta#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("prompt #Delta#lambda_{#theta} (bkg gen with only L/RSB)");
  
  g_lthLSB->SetLineColor(kRed);
  g_lthLSB->SetMarkerColor(kRed);
  g_lthLSB->SetMarkerStyle(20);
  g_lthLSB->SetMarkerSize(.75);
  g_lthLSB->Draw("p same");

  g_lthRSB->SetLineColor(kBlue);
  g_lthRSB->SetMarkerColor(kBlue);
  g_lthRSB->SetMarkerStyle(20);
  g_lthRSB->SetMarkerSize(.75);
  g_lthRSB->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->Draw("ce3");

  c->SaveAs("plots/lth_absDiff_SB.pdf");
  c->Clear();

  // THIRD - draw the mu pt-indep case
  TH1F *fl3 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl3->SetXTitle("p_{T} (GeV)");
  fl3->SetYTitle("#Delta#lambda_{#theta}");
  fl3->GetYaxis()->SetTitleOffset(1.3);
  fl3->GetYaxis()->SetLabelOffset(0.01);
  fl3->SetTitle("prompt #Delta#lambda_{#theta} (#mu_{2} p_{T}-independent)");
  
  g_lthmu2->SetLineColor(kRed);
  g_lthmu2->SetMarkerColor(kRed);
  g_lthmu2->SetMarkerStyle(20);
  g_lthmu2->SetMarkerSize(.75);
  g_lthmu2->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->Draw("ce3");

  c->SaveAs("plots/lth_absDiff_mu2.pdf");
  c->Clear();

  // FOURTH - draw the fnp_psi est case
  TH1F *fl4 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl4->SetXTitle("p_{T} (GeV)");
  fl4->SetYTitle("#Delta#lambda_{#theta}");
  fl4->GetYaxis()->SetTitleOffset(1.3);
  fl4->GetYaxis()->SetLabelOffset(0.01);
  fl4->SetTitle("prompt #Delta#lambda_{#theta} (f_{NP}^{#psi} correction)");
  
  g_lthfnp->SetLineColor(kRed);
  g_lthfnp->SetMarkerColor(kRed);
  g_lthfnp->SetMarkerStyle(20);
  g_lthfnp->SetMarkerSize(.75);
  g_lthfnp->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->Draw("ce3");

  c->SaveAs("plots/lth_absDiff_fnp.pdf");
  c->Clear();

  // now the same for the NP lth
  // SECOND - draw the L/RSB-only gen Dist
  TH1F *flNP2 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  flNP2->SetXTitle("p_{T} (GeV)");
  flNP2->SetYTitle("#Delta#lambda_{#theta}");
  flNP2->GetYaxis()->SetTitleOffset(1.3);
  flNP2->GetYaxis()->SetLabelOffset(0.01);
  flNP2->SetTitle("non-prompt #Delta#lambda_{#theta} (bkg gen with only L/RSB)");
  
  g_lthNPLSB->SetLineColor(kRed);
  g_lthNPLSB->SetMarkerColor(kRed);
  g_lthNPLSB->SetMarkerStyle(20);
  g_lthNPLSB->SetMarkerSize(.75);
  g_lthNPLSB->Draw("p same");

  g_lthNPRSB->SetLineColor(kBlue);
  g_lthNPRSB->SetMarkerColor(kBlue);
  g_lthNPRSB->SetMarkerStyle(20);
  g_lthNPRSB->SetMarkerSize(.75);
  g_lthNPRSB->Draw("p same");

  g_uncNP->SetLineColor(kBlack);
  g_uncNP->SetFillColorAlpha(kBlack, 0.1);
  g_uncNP->Draw("ce3");

  c->SaveAs("plots/lthNP_absDiff_SB.pdf");
  c->Clear();

  // THIRD - draw the mu pt-indep case
  TH1F *flNP3 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  flNP3->SetXTitle("p_{T} (GeV)");
  flNP3->SetYTitle("#Delta#lambda_{#theta}");
  flNP3->GetYaxis()->SetTitleOffset(1.3);
  flNP3->GetYaxis()->SetLabelOffset(0.01);
  flNP3->SetTitle("non-prompt #Delta#lambda_{#theta} (#mu_{2} p_{T}-independent)");
  
  g_lthNPmu2->SetLineColor(kRed);
  g_lthNPmu2->SetMarkerColor(kRed);
  g_lthNPmu2->SetMarkerStyle(20);
  g_lthNPmu2->SetMarkerSize(.75);
  g_lthNPmu2->Draw("p same");

  g_uncNP->SetLineColor(kBlack);
  g_uncNP->SetFillColorAlpha(kBlack, 0.1);
  g_uncNP->Draw("ce3");

  c->SaveAs("plots/lthNP_absDiff_mu2.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();


}
