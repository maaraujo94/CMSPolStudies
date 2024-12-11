// code to compare the deviations for several scenarios
// 1) mass fit w alpha+/-uncertainty
// 2) lifetime fit w NBg+/-uncertainty

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

void plotModel()
{
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - 5 (3) sets
  TGraphErrors **graph_lth = new TGraphErrors*[5];
  TGraphErrors **graph_lthNP = new TGraphErrors*[3];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  graph_lthNP[0] = (TGraphErrors*)fIndB->Get("graph_lambda_NP"); 
  fIndB->Close();
  // 1 - get results for mass fits
  TFile *fIndMPR1 = new TFile("../MPR_alpha_minus/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fIndMPR1->Get("graph_lambda_J");
  fIndMPR1->Close();
  TFile *fIndMPR2 = new TFile("../MPR_alpha_plus/files/finalFitRes.root");
  graph_lth[2] = (TGraphErrors*)fIndMPR2->Get("graph_lambda_J");
  fIndMPR2->Close();
  TFile *fIndMNP1 = new TFile("../MNP_alpha_minus/files/finalFitRes.root");
  graph_lthNP[1] = (TGraphErrors*)fIndMNP1->Get("graph_lambda_NP");
  fIndMNP1->Close();
  TFile *fIndMNP2 = new TFile("../MNP_alpha_plus/files/finalFitRes.root");
  graph_lthNP[2] = (TGraphErrors*)fIndMNP2->Get("graph_lambda_NP");
  fIndMNP2->Close();
  // 2 - get results for lifetime fits
  TFile *fIndLt1 = new TFile("../Lt_Nbg_minus/files/finalFitRes.root");
  graph_lth[3] = (TGraphErrors*)fIndLt1->Get("graph_lambda_J");
  fIndLt1->Close();
  TFile *fIndLt2 = new TFile("../Lt_Nbg_plus/files/finalFitRes.root");
  graph_lth[4] = (TGraphErrors*)fIndLt2->Get("graph_lambda_J");
  fIndLt2->Close();
  
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

    cout << "PR " << i << " : ";
    for(int j = 0; j < 4; j++) {
      cout << diff[j][i] << " ";
    }
    cout << endl;
    cout << "NP " << i << " : ";
    for(int j = 0; j < 2; j++) {
      cout << diffNP[j][i] << " ";
    }
    cout << endl << endl;

  }

     // unc plots should extend to the end of the pT bins, not the middle
  double xv[nBinspT+2], yv[nBinspT+2], xe[nBinspT+2], ye[2][nBinspT+2];
  for(int i = 0; i < nBinspT+2; i++) {
    if(i==0) {
      xv[i] = graph_lth[0]->GetX()[0]-graph_lth[0]->GetEX()[0];
      yv[i] = 0;
      xe[i] = 0;//graph_lth[0]->GetEX()[0];
      ye[0][i] = graph_lth[0]->GetEY()[0];
      ye[1][i] = graph_lthNP[0]->GetEY()[0];
    }
    else if(i==nBinspT+1) {
      xv[i] = graph_lth[0]->GetX()[nBinspT-1]+graph_lth[0]->GetEX()[nBinspT-1];
      yv[i] = 0;
      xe[i] = 0;//graph_lth[0]->GetEX()[nBinspT-1];
      ye[0][i] = graph_lth[0]->GetEY()[nBinspT-1];
      ye[1][i] = graph_lthNP[0]->GetEY()[nBinspT-1];
    }
    else {
      xv[i] = graph_lth[0]->GetX()[i-1];
      yv[i] = 0;
      xe[i] = graph_lth[0]->GetEX()[i-1];
      ye[0][i] = graph_lth[0]->GetEY()[i-1];
      ye[1][i] = graph_lthNP[0]->GetEY()[i-1];
    }
    cout << i << " " << xv[i] << " " << xe[i] << endl;
  }

  TGraphErrors *g_lthMmin = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[0], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthMplus = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[1], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthLtmin = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[2], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthLtplus = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[3], graph_lth[0]->GetEX(), za);
  
  TGraphErrors *g_lthNPMmin = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), diffNP[0], graph_lthNP[0]->GetEX(), za);
  TGraphErrors *g_lthNPMplus = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), diffNP[1], graph_lthNP[0]->GetEX(), za);
  
  // TGraphErrors *g_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_lth[0]->GetEY());
  //TGraphErrors *g_uncNP = new TGraphErrors(nBinspT, graph_lthNP[0]->GetX(), za, graph_lthNP[0]->GetEX(), graph_lthNP[0]->GetEY());
  TGraphErrors *g_unc = new TGraphErrors(nBinspT+2, xv, yv, xe, ye[0]);
  TGraphErrors *g_uncNP = new TGraphErrors(nBinspT+2, xv, yv, xe, ye[1]);

  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.02);

  double d_lim = 60;

  // PROMPT checks
  // FIRST - draw the mass fit check
  double da_lim = 0.19;
  
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT]+5, da_lim);
  fl1->SetXTitle("#it{p}_{T} (GeV)");
  fl1->SetYTitle("#Delta#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.4);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->GetXaxis()->SetTitleOffset(1.1);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->GetXaxis()->CenterTitle(true);
  
  g_lthMmin->SetLineColor(kRed);
  g_lthMmin->SetMarkerColor(kRed);
  g_lthMmin->SetMarkerStyle(20);
  g_lthMmin->SetMarkerSize(.75);
  g_lthMmin->Draw("p same");

  g_lthMplus->SetLineColor(kRed);
  g_lthMplus->SetMarkerColor(kRed);
  g_lthMplus->SetMarkerStyle(24);
  g_lthMplus->SetMarkerSize(.75);
  g_lthMplus->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->SetLineStyle(kDashed);
  g_unc->Draw("ce3");

  TLegend *legAlpha = new TLegend(0.7, 0.85, 1., 0.95);
  legAlpha->SetTextSize(0.03);
  legAlpha->SetBorderSize(0);
  legAlpha->SetFillColorAlpha(kWhite,0);
  legAlpha->AddEntry(g_lthMmin, "#alpha - #sigma_{#alpha}", "pl");
  legAlpha->AddEntry(g_lthMplus, "#alpha + #sigma_{#alpha}", "pl");
  legAlpha->Draw();

  TLatex lc;
  lc.SetTextSize(0.04);
  
  // draw CMS text
  double xp = getPos(pTBins[0]-5, pTBins[nBinspT]+5, 0.1, 0);
  double yp = getPos(-da_lim, da_lim, 0.9, 0);
  lc.DrawLatex(xp, yp, "#bf{prompt J/#psi}");

  c->SaveAs("plots/lth_absDiff_M.pdf");
  c->Clear();

  
  // SECOND - draw the lifetime fit check
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT]+5, da_lim);
  fl2->SetXTitle("#it{p}_{T} (GeV)");
  fl2->SetYTitle("#Delta#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.4);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->GetXaxis()->SetTitleOffset(1.1);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->GetXaxis()->CenterTitle(true);
  
  g_lthLtmin->SetLineColor(kRed);
  g_lthLtmin->SetMarkerColor(kRed);
  g_lthLtmin->SetMarkerStyle(20);
  g_lthLtmin->SetMarkerSize(.75);
  g_lthLtmin->Draw("p same");

  g_lthLtplus->SetLineColor(kRed);
  g_lthLtplus->SetMarkerColor(kRed);
  g_lthLtplus->SetMarkerStyle(24);
  g_lthLtplus->SetMarkerSize(.75);
  g_lthLtplus->Draw("p same");

  g_unc->Draw("ce3");

  TLegend *legNbg = new TLegend(0.7, 0.85, 1., 0.95);
  legNbg->SetTextSize(0.03);
  legNbg->SetBorderSize(0);
  legNbg->SetFillColorAlpha(kWhite,0);
  legNbg->AddEntry(g_lthLtmin, "N_{bkg} - #sigma_{N_{bkg}}", "pl");
  legNbg->AddEntry(g_lthLtplus, "N_{bkg} + #sigma_{N_{bkg}}", "pl");
  legNbg->Draw();

  lc.DrawLatex(xp, yp, "#bf{J/#psi}");

  c->SaveAs("plots/lth_absDiff_Lt.pdf");
  c->Clear();

  // now the same for NP
  // FIRST - draw the mass fit check
  TH1F *flNP1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT]+5, da_lim);
  flNP1->SetXTitle("#it{p}_{T} (GeV)");
  flNP1->SetYTitle("#Delta#lambda_{#theta}");
  flNP1->GetYaxis()->SetTitleOffset(1.4);
  flNP1->GetYaxis()->SetLabelOffset(0.01);
  flNP1->GetXaxis()->SetTitleOffset(1.1);
  flNP1->GetYaxis()->SetLabelOffset(0.01);
  flNP1->GetXaxis()->CenterTitle(true);
  
  g_lthNPMmin->SetLineColor(kRed);
  g_lthNPMmin->SetMarkerColor(kRed);
  g_lthNPMmin->SetMarkerStyle(20);
  g_lthNPMmin->SetMarkerSize(.75);
  g_lthNPMmin->Draw("p same");

  g_lthNPMplus->SetLineColor(kRed);
  g_lthNPMplus->SetMarkerColor(kRed);
  g_lthNPMplus->SetMarkerStyle(24);
  g_lthNPMplus->SetMarkerSize(.75);
  g_lthNPMplus->Draw("p same");

  g_uncNP->SetLineColor(kBlack);
  g_uncNP->SetFillColorAlpha(kBlack, 0.1);
  g_uncNP->SetLineStyle(kDashed);
  g_uncNP->Draw("ce3");

  legAlpha->Draw();

  lc.DrawLatex(xp, yp, "#bf{non-prompt J/#psi}");

  c->SaveAs("plots/lthNP_absDiff_M.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();


}
