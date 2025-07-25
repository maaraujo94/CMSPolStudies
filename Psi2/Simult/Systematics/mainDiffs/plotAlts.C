// code to compare the deviations for several scenarios
// 1) f(pT) - Run2 / Run2
// 2) 1/f(pT) - Run2 / Run2
// 3) pT_cut - Run2 / Run2
// 4) eta_cut - Run2 / Run2
// 5) dR cuts - Run 2 / Run2
// 6) phi reweighing - Run 2 / Run 2

// 0) 2017 - 2018 / Run2
// this case is studied separately because the differences are independent

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

void plotAlts()
{
  double pt_c = 50;
  
  // PART 1 - regular fine-binned results for all studies
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  rHist->SetDirectory(0);
  fIn->Close();
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();

  // get the fit results - 10 sets
  TGraphErrors **graph_lth = new TGraphErrors*[10];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  fIndB->Close();
  // 1 - get 2017 and 2018 results
  TFile *fInd17 = new TFile("../../../2017/PR_fit/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fInd17->Get(Form("graph_lambda_J"));
  fInd17->Close();
  TFile *fInd18 = new TFile("../../../2018/PR_fit/files/finalFitRes.root");
  graph_lth[2] = (TGraphErrors*)fInd18->Get(Form("graph_lambda_J"));
  fInd18->Close();
  // 2 - get lambda values for each eff model
  TFile *fIndf1 = new TFile("../../../Simult_eff1/PR_fit/files/finalFitRes.root");
  graph_lth[3] = (TGraphErrors*)fIndf1->Get(Form("graph_lambda_J"));
  fIndf1->Close();
  TFile *fIndf2 = new TFile("../../../Simult_eff2/PR_fit/files/finalFitRes.root");
  graph_lth[4] = (TGraphErrors*)fIndf2->Get(Form("graph_lambda_J"));
  fIndf2->Close();
  // 3 - get results with higher pT cut
  TFile *fIndpt = new TFile("../../../Simult_cutPt/PR_fit/files/finalFitRes.root");
  graph_lth[5] = (TGraphErrors*)fIndpt->Get(Form("graph_lambda_J"));
  fIndpt->Close();
  // 4 - get results with |eta| cut
  TFile *fIndEta = new TFile("../../../Simult_cutEta/PR_fit/files/finalFitRes.root");
  graph_lth[6] = (TGraphErrors*)fIndEta->Get(Form("graph_lambda_J"));
  fIndEta->Close();
  // 5 - get results with DeltaR_pT cut
  TFile *fIndR = new TFile("../deltaR/files/finalFitRes.root");
  graph_lth[7] = (TGraphErrors*)fIndR->Get(Form("graph_lambda_T"));
  graph_lth[8] = (TGraphErrors*)fIndR->Get(Form("graph_lambda_L"));
  graph_lth[9] = (TGraphErrors*)fIndR->Get(Form("graph_lambda_B"));
  fIndR->Close();
  
  // the phi reweighing is for PR or NP depending on the cuts
  TGraphErrors **graph_phi = new TGraphErrors*[6];

  // get lambda values for each phi reweighing model
  TFile *fIndph = new TFile("../../../Simult/PR_fit/files/finalFitRes.root");
  graph_phi[0] = (TGraphErrors*)fIndph->Get(Form("graph_lambda_J"));
  graph_phi[1] = (TGraphErrors*)fIndph->Get(Form("graph_lambda_NP"));
  fIndph->Close();
  TFile *fIndph1 = new TFile("../../../Simult_phi1/PR_fit/files/finalFitRes.root");
  graph_phi[2] = (TGraphErrors*)fIndph1->Get(Form("graph_lambda_J"));
  fIndph1->Close();
  TFile *fIndph2 = new TFile("../../../Simult_phi2/PR_fit/files/finalFitRes.root");
  graph_phi[3] = (TGraphErrors*)fIndph2->Get(Form("graph_lambda_J"));
  fIndph2->Close();
  TFile *fIndph3 = new TFile("../../../Simult_phi3/PR_fit/files/finalFitRes.root");
  graph_phi[4] = (TGraphErrors*)fIndph3->Get(Form("graph_lambda_NP"));
  fIndph3->Close();
  TFile *fIndph4 = new TFile("../../../Simult_phi4/PR_fit/files/finalFitRes.root");
  graph_phi[5] = (TGraphErrors*)fIndph4->Get(Form("graph_lambda_NP"));
  fIndph4->Close();

  // get the differences
  double diff[5][nBinspT], za[nBinspT], uncR[nBinspT];
  double diffP[4][nBinspT];
  double diffY[nBinspT], errY[nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    diff[0][i] = (graph_lth[3]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[1][i] = (graph_lth[4]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[2][i] = (graph_lth[5]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[3][i] = (graph_lth[6]->GetY()[i] - graph_lth[0]->GetY()[i]);
    // deltaR cuts apply to different pT regions
    if(graph_lth[0]->GetX()[i] < pt_c) {
      diff[4][i] = (graph_lth[7]->GetY()[i] - graph_lth[9]->GetY()[i]);
      double unc1 = graph_lth[7]->GetEY()[i];
      double unc2 = graph_lth[9]->GetEY()[i];
      uncR[i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));
    }
    else {
      diff[4][i] = (graph_lth[8]->GetY()[i] - graph_lth[9]->GetY()[i]);
      double unc1 = graph_lth[8]->GetEY()[i];
      double unc2 = graph_lth[9]->GetEY()[i];
      uncR[i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));
    }
    // 2017 - 2018 case can be calculated w their uncertainties
    double val7 = graph_lth[1]->GetY()[i];
    double val8 = graph_lth[2]->GetY()[i];
    diffY[i] = (val7 - val8);
    double err7 = graph_lth[1]->GetEY()[i];
    double err8 = graph_lth[2]->GetEY()[i];
    errY[i] = sqrt(pow(err7, 2) + pow(err8, 2));

    //phi part done separately
    diffP[0][i] = (graph_phi[2]->GetY()[i] - graph_phi[0]->GetY()[i]);
    diffP[1][i] = (graph_phi[3]->GetY()[i] - graph_phi[0]->GetY()[i]);
    diffP[2][i] = (graph_phi[4]->GetY()[i] - graph_phi[1]->GetY()[i]);
    diffP[3][i] = (graph_phi[5]->GetY()[i] - graph_phi[1]->GetY()[i]);

    cout << i << " " << diffP[0][i] << " " << diffP[3][i] << endl;
    
    za[i] = 0;
  }

  // unc plots should extend to the end of the pT bins, not the middle
  double xv[nBinspT+2], yv[nBinspT+2], xe[nBinspT+2], ye[2][nBinspT+2];
  for(int i = 0; i < nBinspT+2; i++) {
    if(i==0) {
      xv[i] = graph_lth[0]->GetX()[0]-graph_lth[0]->GetEX()[0];
      yv[i] = 0;
      xe[i] = 0;//graph_lth[0]->GetEX()[0];
      ye[0][i] = graph_lth[0]->GetEY()[0];
      ye[1][i] = graph_phi[1]->GetEY()[0];
    }
    else if(i==nBinspT+1) {
      xv[i] = graph_lth[0]->GetX()[nBinspT-1]+graph_lth[0]->GetEX()[nBinspT-1];
      yv[i] = 0;
      xe[i] = 0;//graph_lth[0]->GetEX()[nBinspT-1];
      ye[0][i] = graph_lth[0]->GetEY()[nBinspT-1];
      ye[1][i] = graph_phi[1]->GetEY()[nBinspT-1];
    }
    else {
      xv[i] = graph_lth[0]->GetX()[i-1];
      yv[i] = 0;
      xe[i] = graph_lth[0]->GetEX()[i-1];
      ye[0][i] = graph_lth[0]->GetEY()[i-1];
      ye[1][i] = graph_phi[1]->GetEY()[i-1];
    }
    cout << i << " " << xv[i] << " " << xe[i] << endl;
  }

  // all the tgraphs with the deviations
  TGraphErrors *g_lthY = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diffY, graph_lth[0]->GetEX(), errY);
  TGraphErrors *g_lthF = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[0], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthFI = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[1], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthpT = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[2], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthEta = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[3], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthR = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[4], graph_lth[0]->GetEX(), uncR);
  TGraphErrors *g_lthPhi1 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diffP[0], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthPhi2 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diffP[1], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthPhiNP1 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diffP[2], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthPhiNP2 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diffP[3], graph_lth[0]->GetEX(), za);
  
  //TGraphErrors *g_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_lth[0]->GetEY());
  TGraphErrors *g_unc = new TGraphErrors(nBinspT+2, xv, yv, xe, ye[0]);
  // TGraphErrors *g_uncNP = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_phi[1]->GetEY());
  TGraphErrors *g_uncNP = new TGraphErrors(nBinspT+2, xv, yv, xe, ye[1]);
 
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.02);
  
  // FIRST - draw the abs diff + Simult unc band
  double da_lim = 0.29;
  
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT]+5, da_lim);
  fl1->SetXTitle("#it{p}_{T} (GeV)");
  fl1->SetYTitle("#Delta#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.4);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->GetXaxis()->SetTitleOffset(1.1);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->GetXaxis()->CenterTitle(true);
  
  g_lthF->SetLineColor(kBlue);
  g_lthF->SetMarkerColor(kBlue);
  g_lthF->SetMarkerStyle(24);
  g_lthF->SetMarkerSize(.75);
  g_lthF->Draw("p same");

  g_lthFI->SetLineColor(kBlue);
  g_lthFI->SetMarkerColor(kBlue);
  g_lthFI->SetMarkerStyle(22);
  g_lthFI->SetMarkerSize(.75);
  g_lthFI->Draw("p same");

  g_lthpT->SetLineColor(kGreen+1);
  g_lthpT->SetMarkerColor(kGreen+1);
  g_lthpT->SetMarkerStyle(20);
  g_lthpT->SetMarkerSize(.75);
  g_lthpT->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->SetLineStyle(kDashed);
  g_unc->Draw("ce3");

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT]+5, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);

  TLegend *legmu = new TLegend(0.7, 0.825, 1., 0.975);
  legmu->SetTextSize(0.03);
  legmu->SetBorderSize(0);
  legmu->SetFillColorAlpha(kWhite,0);
  legmu->AddEntry(g_lthF, "weight f(#it{p}_{T})", "pl");
  legmu->AddEntry(g_lthFI, "weight 1/f(#it{p}_{T})", "pl");
  legmu->AddEntry(g_lthpT, "#it{p}_{T} cut", "pl");
  legmu->Draw();

  TLatex lc;
  lc.SetTextSize(0.04);
  
  // draw CMS text
  double xp = getPos(pTBins[0]-5, pTBins[nBinspT]+5, 0.1, 0);
  double yp = getPos(-da_lim, da_lim, 0.9, 0);;
  lc.DrawLatex(xp, yp, "#bf{#psi(2S)}");

  c->SaveAs("plots/lth_absDiff_muEff.pdf");
  c->Clear();

  // eta variation
  TH1F *flE = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT]+5, da_lim);
  flE->SetXTitle("#it{p}_{T} (GeV)");
  flE->SetYTitle("#Delta#lambda_{#theta}");
  flE->GetYaxis()->SetTitleOffset(1.4);
  flE->GetYaxis()->SetLabelOffset(0.01);
  flE->GetXaxis()->SetTitleOffset(1.1);
  flE->GetYaxis()->SetLabelOffset(0.01);
  flE->GetXaxis()->CenterTitle(true);
  
  g_lthEta->SetLineColor(kRed);
  g_lthEta->SetMarkerColor(kRed);
  g_lthEta->SetMarkerStyle(20);
  g_lthEta->SetMarkerSize(.75);
  g_lthEta->Draw("p same");

  g_unc->Draw("ce3");

  TLegend *legWheel = new TLegend(0.7, 0.9, 1., 0.95);
  legWheel->SetTextSize(0.03);
  legWheel->SetBorderSize(0);
  legWheel->SetFillColorAlpha(kWhite,0);
  legWheel->AddEntry(g_lthEta, "#eta cut", "pl");
  legWheel->Draw();

  // draw CMS text
  lc.DrawLatex(xp, yp, "#bf{#psi(2S)}");

  c->SaveAs("plots/lth_absDiff_eta.pdf");
  c->Clear();

  // FIRST (2) - draw just deltaR (larger difference)
  TH1F *fl12 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT]+5, da_lim);
  fl12->SetXTitle("#it{p}_{T} (GeV)");
  fl12->SetYTitle("#Delta#lambda_{#theta}");
  fl12->GetYaxis()->SetTitleOffset(1.4);
  fl12->GetYaxis()->SetLabelOffset(0.01);
  fl12->GetXaxis()->SetTitleOffset(1.1);
  fl12->GetYaxis()->SetLabelOffset(0.01);
  fl12->GetXaxis()->CenterTitle(true);
  
  g_lthR->SetLineColor(kBlue);
  g_lthR->SetMarkerColor(kBlue);
  g_lthR->SetMarkerStyle(20);
  g_lthR->SetMarkerSize(.75);
  g_lthR->Draw("p same");

  zero->Draw();

  TLine *trans22A = new TLine(pt_c, -da_lim, pt_c, da_lim);
  trans22A->SetLineColor(kBlack);
  trans22A->SetLineStyle(kDashed);
  trans22A->Draw();

  TLegend *legrho = new TLegend(0.6, 0.9, 0.9, 0.95);
  legrho->SetTextSize(0.03);
  legrho->SetBorderSize(0);
  legrho->SetFillColorAlpha(kWhite,0);
  legrho->AddEntry(g_lthR, "#DeltaR_{#it{p}_{T}} cut", "pl");
  legrho->Draw();

  // draw CMS text
  yp = getPos(-da_lim, da_lim, 0.9, 0);
  lc.DrawLatex(xp, yp, "#bf{#psi(2S)}");

  c->SaveAs("plots/lth_absDiff_rho.pdf");
  c->Clear();
  
  // SECOND - draw 2017-2018 with the calculated uncertainty
  da_lim = 0.59;
  c->SetTopMargin(0.02);

  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT]+5, da_lim);
  fl2->SetXTitle("#it{p}_{T} (GeV)");
  fl2->SetYTitle("#Delta#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.4);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->GetXaxis()->SetTitleOffset(1.1);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->GetXaxis()->CenterTitle(true);
  
  g_lthY->SetLineColor(kBlack);
  g_lthY->SetMarkerColor(kBlack);
  g_lthY->SetMarkerStyle(20);
  g_lthY->SetMarkerSize(.75);
  g_lthY->Draw("p same");

  zero->Draw();
  
  TLegend* legY = new TLegend(0.7, 0.9, 1, 0.95);
  legY->SetTextSize(0.03);
  legY->SetBorderSize(0);
  legY->SetFillColorAlpha(kWhite,0);
  legY->AddEntry(g_lthY, "2017-2018", "pl");
  legY->Draw();

  // draw CMS text
  xp = getPos(pTBins[0]-5, pTBins[nBinspT]+5, 0.1, 0);
  yp = 0.45;
  lc.DrawLatex(xp, yp, "#bf{#psi(2S)}");

  c->SaveAs("plots/lth_absDiff_Y.pdf");
  c->Clear();

  // phi-based variation
  da_lim = 0.29;

  TH1F *fl3 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT]+5, da_lim);
  fl3->SetXTitle("#it{p}_{T} (GeV)");
  fl3->SetYTitle("#Delta#lambda_{#theta}");
  fl3->GetYaxis()->SetTitleOffset(1.4);
  fl3->GetYaxis()->SetLabelOffset(0.01);
  fl3->GetXaxis()->SetTitleOffset(1.1);
  fl3->GetYaxis()->SetLabelOffset(0.01);
  fl3->GetXaxis()->CenterTitle(true);
  
  g_lthPhi2->SetLineColor(kBlue);
  g_lthPhi2->SetMarkerColor(kBlue);
  g_lthPhi2->SetMarkerStyle(20);
  g_lthPhi2->SetMarkerSize(.75);
  g_lthPhi2->Draw("p same");

  g_lthPhi1->SetLineColor(kRed);
  g_lthPhi1->SetMarkerColor(kRed);
  g_lthPhi1->SetMarkerStyle(20);
  g_lthPhi1->SetMarkerSize(.75);
  g_lthPhi1->Draw("p same");

  g_unc->Draw("ce3");

  TLegend *legp = new TLegend(0.7, 0.85, 1., 0.95);
  legp->SetTextSize(0.03);
  legp->SetBorderSize(0);
  legp->SetFillColorAlpha(kWhite,0);
  legp->AddEntry(g_lthPhi2, "#beta = -0.01", "pl");
  legp->AddEntry(g_lthPhi1, "#beta = -0.02", "pl");
  legp->Draw();

  // draw CMS text
  xp = getPos(pTBins[0]-5, pTBins[nBinspT]+5, 0.1, 0);
  yp = getPos(-da_lim, da_lim, 0.9, 0);;
  lc.DrawLatex(xp, yp, "#bf{prompt #psi(2S)}");

  c->SaveAs("plots/lth_absDiff_phi.pdf");
  c->Clear();

  TH1F *fl4 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT]+5, da_lim);
  fl4->SetXTitle("#it{p}_{T} (GeV)");
  fl4->SetYTitle("#Delta#lambda_{#theta}");
  fl4->GetYaxis()->SetTitleOffset(1.4);
  fl4->GetYaxis()->SetLabelOffset(0.01);
  fl4->GetXaxis()->SetTitleOffset(1.1);
  fl4->GetYaxis()->SetLabelOffset(0.01);
  fl4->GetXaxis()->CenterTitle(true);
  
  g_lthPhiNP1->SetLineColor(kBlue);
  g_lthPhiNP1->SetMarkerColor(kBlue);
  g_lthPhiNP1->SetMarkerStyle(20);
  g_lthPhiNP1->SetMarkerSize(.75);
  g_lthPhiNP1->Draw("p same");

  g_lthPhiNP2->SetLineColor(kRed);
  g_lthPhiNP2->SetMarkerColor(kRed);
  g_lthPhiNP2->SetMarkerStyle(20);
  g_lthPhiNP2->SetMarkerSize(.75);
  g_lthPhiNP2->Draw("p same");

  g_uncNP->SetLineColor(kBlack);
  g_uncNP->SetFillColorAlpha(kBlack, 0.1);
  g_uncNP->SetLineStyle(kDashed);
  g_uncNP->Draw("ce3");

  TLegend *legpn = new TLegend(0.7, 0.85, 1., 0.95);
  legpn->SetTextSize(0.03);
  legpn->SetBorderSize(0);
  legpn->SetFillColorAlpha(kWhite,0);
  legpn->AddEntry(g_lthPhiNP1, "#beta = 0.005", "pl");
  legpn->AddEntry(g_lthPhiNP2, "#beta = 0.025", "pl");
  legpn->Draw();

  // draw CMS text
  lc.DrawLatex(xp, yp, "#bf{non-prompt #psi(2S)}");

  c->SaveAs("plots/lthNP_absDiff_phi.pdf");
  c->Clear();

  c->Destructor();

}
