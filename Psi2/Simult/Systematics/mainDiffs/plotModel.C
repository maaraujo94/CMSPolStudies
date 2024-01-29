// code to compare the deviations for change in model
// linear mass bkg - Run2
// fSB through evt counting - Run 2
// lifetime fit with NBg+-uncertainty

void plotModel()
{
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - 3 sets for PR / 5 for NP
  TGraphErrors **graph_lth = new TGraphErrors*[3];
  TGraphErrors **graph_lthNP = new TGraphErrors*[3];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  graph_lthNP[0] = (TGraphErrors*)fIndB->Get("graph_lambda_NP");
  fIndB->Close();
  // 1 - get lambda values for linear mass bkg
  TFile *fIndmb = new TFile("../Mass_bkg_alt/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fIndmb->Get(Form("graph_lambda_J"));
  graph_lthNP[1] = (TGraphErrors*)fIndmb->Get("graph_lambda_NP");
  fIndmb->Close();
  // 2 - get lambda values for evt counting
  TFile *fIndct = new TFile("../fBg_evt_ct/PR_fit/files/finalFitRes.root");
  graph_lth[2] = (TGraphErrors*)fIndct->Get(Form("graph_lambda_J"));
  graph_lthNP[2] = (TGraphErrors*)fIndct->Get("graph_lambda_NP");
  fIndct->Close();
  // 3 - get lambda values for lt fits
  TFile *fIndlt1 = new TFile("../Lt_Nbg_minus/files/finalFitRes.root");
  graph_lth[3] = (TGraphErrors*)fIndlt1->Get(Form("graph_lambda_J"));
  fIndlt1->Close();
  TFile *fIndlt2 = new TFile("../Lt_Nbg_plus/files/finalFitRes.root");
  graph_lth[4] = (TGraphErrors*)fIndlt2->Get(Form("graph_lambda_J"));
  fIndlt2->Close();

  // get the differences
  double diff[4][nBinspT], diffNP[2][nBinspT], za[nBinspT];
  double err[4][nBinspT], errNP[2][nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    diff[0][i] = (graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[1][i] = (graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[2][i] = (graph_lth[3]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[3][i] = (graph_lth[4]->GetY()[i] - graph_lth[0]->GetY()[i]);
 
    double unc1 = graph_lth[0]->GetEY()[i];
    double unc2 = graph_lth[1]->GetEY()[i];
    err[0][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lth[2]->GetEY()[i];
    err[1][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lth[3]->GetEY()[i];
    err[2][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lth[4]->GetEY()[i];
    err[3][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    diffNP[0][i] = (graph_lthNP[1]->GetY()[i] - graph_lthNP[0]->GetY()[i]);
    diffNP[1][i] = (graph_lthNP[2]->GetY()[i] - graph_lthNP[0]->GetY()[i]);
 
    unc1 = graph_lthNP[0]->GetEY()[i];
    unc2 = graph_lthNP[1]->GetEY()[i];
    errNP[0][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    unc2 = graph_lthNP[2]->GetEY()[i];
    errNP[1][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    za[i] = 0;
  }
  TGraphErrors *g_lthMB = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[0], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthCt = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[1], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthLt1 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[2], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthLt2 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[3], graph_lth[0]->GetEX(), za);

  TGraphErrors *g_lthNPMB = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), diffNP[0], graph_lthNP[0]->GetEX(), za);
  TGraphErrors *g_lthNPCt = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), diffNP[1], graph_lthNP[0]->GetEX(), za);

  TGraphErrors *g_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_lth[0]->GetEY());
  TGraphErrors *g_uncNP = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), za, graph_lthNP[0]->GetEX(), graph_lthNP[0]->GetEY());
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  
  double d_lim = 60;

  // FIRST - draw the abs diff + Simult unc band for lin M bkg
  double da_lim = 0.3;
  
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#Delta#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("prompt #Delta#lambda_{#theta} (mass bkg subtraction)");
  
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

  // SECOND - same for evt ct check
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#Delta#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("prompt #Delta#lambda_{#theta} (event counting)");
  
  g_lthCt->SetLineColor(kRed);
  g_lthCt->SetMarkerColor(kRed);
  g_lthCt->SetMarkerStyle(20);
  g_lthCt->SetMarkerSize(.75);
  g_lthCt->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->Draw("ce3");

  c->SaveAs("plots/lth_absDiff_Ct.pdf");
  c->Clear();

  // THIRD - same for lt fits with new NBg
  TH1F *fl3 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl3->SetXTitle("p_{T} (GeV)");
  fl3->SetYTitle("#Delta#lambda_{#theta}");
  fl3->GetYaxis()->SetTitleOffset(1.3);
  fl3->GetYaxis()->SetLabelOffset(0.01);
  fl3->SetTitle("prompt #Delta#lambda_{#theta} (varying N_{bkg})");
  
  g_lthLt1->SetLineColor(kRed);
  g_lthLt1->SetMarkerColor(kRed);
  g_lthLt1->SetMarkerStyle(20);
  g_lthLt1->SetMarkerSize(.75);
  g_lthLt1->Draw("p same");

  g_lthLt2->SetLineColor(kRed);
  g_lthLt2->SetMarkerColor(kRed);
  g_lthLt2->SetMarkerStyle(24);
  g_lthLt2->SetMarkerSize(.75);
  g_lthLt2->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->Draw("ce3");

  c->SaveAs("plots/lth_absDiff_Lt.pdf");
  c->Clear();

  // now we do the same for lth NP
  // FIRST - lin M bkg (NP)
  TH1F *flNP1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  flNP1->SetXTitle("p_{T} (GeV)");
  flNP1->SetYTitle("#Delta#lambda_{#theta}");
  flNP1->GetYaxis()->SetTitleOffset(1.3);
  flNP1->GetYaxis()->SetLabelOffset(0.01);
  flNP1->SetTitle("non-prompt #Delta#lambda_{#theta} (mass bkg subtraction)");
  
  g_lthNPMB->SetLineColor(kRed);
  g_lthNPMB->SetMarkerColor(kRed);
  g_lthNPMB->SetMarkerStyle(20);
  g_lthNPMB->SetMarkerSize(.75);
  g_lthNPMB->Draw("p same");

  g_uncNP->SetLineColor(kBlack);
  g_uncNP->SetFillColorAlpha(kBlack, 0.1);
  g_uncNP->Draw("ce3");

  c->SaveAs("plots/lthNP_absDiff_M.pdf");
  c->Clear();

  // SECOND - evt ct check (NP)
  TH1F *flNP2 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  flNP2->SetXTitle("p_{T} (GeV)");
  flNP2->SetYTitle("#Delta#lambda_{#theta}");
  flNP2->GetYaxis()->SetTitleOffset(1.3);
  flNP2->GetYaxis()->SetLabelOffset(0.01);
  flNP2->SetTitle("non-prompt #Delta#lambda_{#theta} (event counting)");
  
  g_lthNPCt->SetLineColor(kRed);
  g_lthNPCt->SetMarkerColor(kRed);
  g_lthNPCt->SetMarkerStyle(20);
  g_lthNPCt->SetMarkerSize(.75);
  g_lthNPCt->Draw("p same");

  g_uncNP->SetLineColor(kBlack);
  g_uncNP->SetFillColorAlpha(kBlack, 0.1);
  g_uncNP->Draw("ce3");

  c->SaveAs("plots/lthNP_absDiff_Ct.pdf");
  c->Clear();
  c->Destructor();

  fIn->Close();


}
