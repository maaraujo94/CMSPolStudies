// code to compare the deviations for several scenarios
// 1) f(pT) - Run2 / Run2
// 2) 1/f(pT) - Run2 / Run2
// 3) pT_cut - Run2 / Run2
// 4) eta_cut - Run2 / Run2
// 5) dR cuts - Run 2 / Run2

// 0) 2017 - 2018 / Run2
// this case is studied separately because the differences are independent

void plotRes()
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
  TGraphErrors **graph_lth = new TGraphErrors*[9];
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
  TFile *fIndR1 = new TFile("../../../Simult_dR1/PR_fit/files/finalFitRes.root");
  graph_lth[7] = (TGraphErrors*)fIndR1->Get(Form("graph_lambda_J"));
  fIndR1->Close();
  TFile *fIndR2 = new TFile("../../../Simult_dR2/PR_fit/files/finalFitRes.root");
  graph_lth[8] = (TGraphErrors*)fIndR2->Get(Form("graph_lambda_J"));
  fIndR2->Close();

  // get the differences
  double diff[5][nBinspT], za[nBinspT], uncR[nBinspT];
  double diffY[3][nBinspT], errY[nBinspT], sigY[nBinspT];
  double v_int = 0;
  for(int i = 0; i < nBinspT; i++) {
    diff[0][i] = (graph_lth[3]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[1][i] = (graph_lth[4]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[2][i] = (graph_lth[5]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[3][i] = (graph_lth[6]->GetY()[i] - graph_lth[0]->GetY()[i]);
    // deltaR cuts apply to different pT regions
    if(graph_lth[0]->GetX()[i] < 66) {
      diff[4][i] = (graph_lth[7]->GetY()[i] - graph_lth[0]->GetY()[i]);
      uncR[i] = graph_lth[7]->GetEY()[i];
    }
    else {
      diff[4][i] = (graph_lth[8]->GetY()[i] - graph_lth[0]->GetY()[i]);
      uncR[i] = graph_lth[8]->GetEY()[i];
    }
    // 2017 - 2018 case can be calculated w their uncertainties
    double val7 = graph_lth[1]->GetY()[i];
    double val8 = graph_lth[2]->GetY()[i];
    diffY[0][i] = (val7 - val8);
    double err7 = graph_lth[1]->GetEY()[i];
    double err8 = graph_lth[2]->GetEY()[i];
    errY[i] = sqrt(pow(err7, 2) + pow(err8, 2));
    // also calculate 2017-2018 sigmas / pulls
    sigY[i] = diffY[0][i]/errY[i];
    // 2017-run2 and 2018-run2
    diffY[1][i] = val7-graph_lth[0]->GetY()[i];
    diffY[2][i] = val8-graph_lth[0]->GetY()[i];
    
    za[i] = 0;
    v_int += abs(diffY[0][i])/nBinspT;
  }
  v_int /= 2.;
  double d_v[nBinspT];
  for(int i = 0; i < nBinspT; i++) d_v[i] = v_int;

   
  // PART 2 - coarse-binned results for rho study
  // get the histo limits
  TFile *fInc = new TFile("../deltaR/files/chistStore.root");
  TH2D* cHist;
  fInc->GetObject("cHistB", cHist);
  cHist->SetDirectory(0);
  fInc->Close();
  
  int nBinspT_c = cHist->GetNbinsY();
  const double *pTBins_c = cHist->GetYaxis()->GetXbins()->GetArray();

  // get the fit results - 3 sets
  TGraphErrors **graph_lth_c = new TGraphErrors*[3];
  // 0 - get Run2 results
  TFile *fIndB_c = new TFile("../deltaR/files/finalFitRes.root");
  graph_lth_c[0] = (TGraphErrors*)fIndB_c->Get("graph_lambda_B");
  graph_lth_c[1] = (TGraphErrors*)fIndB_c->Get(Form("graph_lambda_T"));
  graph_lth_c[2] = (TGraphErrors*)fIndB_c->Get(Form("graph_lambda_L"));
  fIndB_c->Close();

  // get the differences
  double diff_c[nBinspT], za_c[nBinspT];
  for(int i = 0; i < nBinspT_c; i++) {
    // deltaR cuts apply to different pT regions
    if(graph_lth_c[0]->GetX()[i] < 66) {
      diff_c[i] = (graph_lth_c[1]->GetY()[i] - graph_lth_c[0]->GetY()[i]);
    }
    else {
      diff_c[i] = (graph_lth_c[2]->GetY()[i] - graph_lth_c[0]->GetY()[i]);
    }

    za_c[i] = 0;
  }

  TGraphErrors *g_lthY = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diffY[0], graph_lth[0]->GetEX(), errY);
  TGraphErrors *g_lthY7 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diffY[1], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthY8 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diffY[2], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthF = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[0], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthFI = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[1], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthpT = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[2], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthEta = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[3], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthR = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[4], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthR_c = new TGraphErrors(nBinspT_c, graph_lth_c[0]->GetX(), diff_c, graph_lth_c[0]->GetEX(), za_c);
  
  TGraphErrors *g_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_lth[0]->GetEY());
  TGraph *g_sigY = new TGraph(nBinspT, graph_lth[0]->GetX(), sigY);
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  double d_lim = 60;

  // FIRST - draw the abs diff + Simult unc band
  double da_lim = 0.3;
  
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#Delta#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("#Delta#lambda_{#theta} (muon eff |#eta|<0.2)");
  
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
  g_unc->SetFillColorAlpha(kBlack, 0.3);
  g_unc->Draw("ce3");

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);

  TLegend *leg = new TLegend(0.65, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(g_lthF, "weight f(p_{T})", "pl");
  leg->AddEntry(g_lthFI, "weight 1/f(p_{T})", "pl");
  leg->AddEntry(g_lthpT, "p_{T} > 6.7 GeV", "pl");
  leg->Draw();

  c->SaveAs("lth_absDiff_band.pdf");
  c->Clear();

  TH1F *flE = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  flE->SetXTitle("p_{T} (GeV)");
  flE->SetYTitle("#Delta#lambda_{#theta}");
  flE->GetYaxis()->SetTitleOffset(1.3);
  flE->GetYaxis()->SetLabelOffset(0.01);
  flE->SetTitle("#Delta#lambda_{#theta} (muon eff 0.2<|#eta|<0.3)");
  
  g_lthEta->SetLineColor(kRed);
  g_lthEta->SetMarkerColor(kRed);
  g_lthEta->SetMarkerStyle(20);
  g_lthEta->SetMarkerSize(.75);
  g_lthEta->Draw("p same");

  g_unc->Draw("ce3");

  c->SaveAs("lth_absDiff_eta.pdf");
  c->Clear();

  // FIRST (2) - draw just deltaR (larger difference)
  da_lim = 0.3;
  
  TH1F *fl12 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl12->SetXTitle("p_{T} (GeV)");
  fl12->SetYTitle("#Delta#lambda_{#theta}");
  fl12->GetYaxis()->SetTitleOffset(1.3);
  fl12->GetYaxis()->SetLabelOffset(0.01);
  fl12->SetTitle("#Delta#lambda_{#theta} (#rho factor)");
  
  g_lthR->SetLineColor(kBlue);
  g_lthR->SetMarkerColor(kBlue);
  g_lthR->SetMarkerStyle(20);
  g_lthR->SetMarkerSize(.75);
  g_lthR->Draw("p same");

  g_lthR_c->SetLineColor(kRed);
  g_lthR_c->SetMarkerColor(kRed);
  g_lthR_c->SetMarkerStyle(20);
  g_lthR_c->SetMarkerSize(.75);
  g_lthR_c->Draw("p same");

  g_unc->Draw("ce3");

  TLine *trans22A = new TLine(66, -da_lim, 66, da_lim);
  trans22A->SetLineColor(kBlack);
  trans22A->SetLineStyle(kDashed);
  trans22A->Draw();

  TLegend *legR = new TLegend(0.6, 0.75, 0.9, 0.9);
  legR->SetTextSize(0.03);
  legR->AddEntry(g_lthR, "#DeltaR cut (fine)", "pl");
  legR->AddEntry(g_lthR_c, "#DeltaR cut (coarse)", "pl");
  legR->Draw();

  c->SaveAs("lth_absDiff_rho.pdf");
  c->Clear();

  
  // SECOND - draw 2017-2018 with the calculated uncertainty
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#Delta#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("#Delta#lambda_{#theta} (2017-2018)");
  
  g_lthY->SetLineColor(kBlack);
  g_lthY->SetMarkerColor(kBlack);
  g_lthY->SetMarkerStyle(20);
  g_lthY->SetMarkerSize(.75);
  g_lthY->Draw("p same");

  zero->Draw();
  
  c->SaveAs("lth_absDiff_Y.pdf");

  TGraphErrors *g_year = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), d_v);
  g_year->SetLineColor(kBlack);
  g_year->SetFillColorAlpha(kBlack, 0.3);
  g_year->Draw("ce3");

  c->SaveAs("lth_Y_band.pdf");
  c->Clear();
  
  // get 2017, 2018 avg diff from run2
  double vy_sum[2] = {0,0};
  double vy_int[2] = {0,0};
  for(int i = 0; i < nBinspT; i++) {
    vy_sum[0] += diffY[1][i]/nBinspT;
    vy_sum[1] += diffY[2][i]/nBinspT;
    vy_int[0] += diffY[1][i]*(pTBins[i+1]-pTBins[i])/(pTBins[nBinspT]-pTBins[0]);
    vy_int[1] += diffY[2][i]*(pTBins[i+1]-pTBins[i])/(pTBins[nBinspT]-pTBins[0]);
  }
  cout << "2017: avg = " << vy_sum[0] << "; int = " << vy_int[0] << endl;
  cout << "2018: avg = " << vy_sum[1] << "; int = " << vy_int[1] << endl;

  cout << "symm integrated difference: " << 0.5*(abs(vy_int[0])+abs(vy_int[1])) << endl;
  
  //draw 2017-run2, 2018-run2
  TH1F *fl78 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl78->SetXTitle("p_{T} (GeV)");
  fl78->SetYTitle("#Delta#lambda_{#theta}");
  fl78->GetYaxis()->SetTitleOffset(1.3);
  fl78->GetYaxis()->SetLabelOffset(0.01);
  fl78->SetTitle("#Delta#lambda_{#theta} (2017,2018)");
  
  g_lthY7->SetLineColor(kBlue);
  g_lthY7->SetMarkerColor(kBlue);
  g_lthY7->SetMarkerStyle(20);
  g_lthY7->SetMarkerSize(.75);
  g_lthY7->Draw("p same");

  g_lthY8->SetLineColor(kBlack);
  g_lthY8->SetMarkerColor(kBlack);
  g_lthY8->SetMarkerStyle(20);
  g_lthY8->SetMarkerSize(.75);
  g_lthY8->Draw("p same");

  zero->Draw();

  /*  TF1 *favg17 = new TF1("favg17", "[0]", pTBins[0]-5, pTBins[nBinspT]);
  favg17->SetParameter(0, vy_sum[0]);
  favg17->SetLineColor(kBlue);
  favg17->Draw("lsame");
  TF1 *favg18 = new TF1("favg18", "[0]", pTBins[0]-5, pTBins[nBinspT]);
  favg18->SetParameter(0, vy_sum[1]);
  favg18->SetLineColor(kRed);  
  favg18->Draw("lsame");*/
 
  TF1 *fint17 = new TF1("fint17", "[0]", pTBins[0]-5, pTBins[nBinspT]);
  fint17->SetParameter(0, vy_int[0]);
  fint17->SetLineStyle(kDashed);
  fint17->SetLineColor(kBlue);  
  // fint17->Draw("lsame");  
  TF1 *fint18 = new TF1("fint18", "[0]", pTBins[0]-5, pTBins[nBinspT]);
  fint18->SetParameter(0, vy_int[1]); 
  fint18->SetLineStyle(kDashed);
  fint18->SetLineColor(kRed);   
  // fint18->Draw("lsame");

  TLegend *leg78 = new TLegend(0.75, 0.75, 0.9, 0.9);
  leg78->SetTextSize(0.03);
  leg78->AddEntry(g_lthY7, "2017", "pl");
  leg78->AddEntry(g_lthY8, "2018", "pl");
  leg78->Draw();

  c->SaveAs("lth_absDiff_Y78.pdf");
  c->Clear();
    
  // THIRD - Draw 2017-2018 as pulls / nr sigma (divide values by unc)
  TH1F *fl3 = c->DrawFrame(pTBins[0]-5, -3, pTBins[nBinspT], 3);
  fl3->SetXTitle("p_{T} (GeV)");
  fl3->SetYTitle("pulls");
  fl3->GetYaxis()->SetTitleOffset(1.3);
  fl3->GetYaxis()->SetLabelOffset(0.01);
  fl3->SetTitle("#Delta#lambda_{#theta} pulls (2017-2018)");
  
  g_sigY->SetLineColor(kBlack);
  g_sigY->SetMarkerColor(kBlack);
  g_sigY->SetMarkerStyle(20);
  g_sigY->SetMarkerSize(1);
  g_sigY->Draw("p same");

  TLine *trans1S = new TLine(46, -3, 46, 3);
  trans1S->SetLineColor(kBlack);
  trans1S->SetLineStyle(kDashed);
  //trans1S->Draw();
  TLine *trans2S = new TLine(66, -3, 66, 3);
  trans2S->SetLineColor(kBlack);
  trans2S->SetLineStyle(kDashed);
  //trans2S->Draw();
  zero->Draw();
  TLine *oneM = new TLine(pTBins[0]-5, -1, pTBins[nBinspT], -1);
  oneM->SetLineColor(kBlack);
  oneM->SetLineStyle(kDashDotted);
  oneM->Draw();
  TLine *oneP = new TLine(pTBins[0]-5, 1, pTBins[nBinspT], 1);
  oneP->SetLineColor(kBlack);
  oneP->SetLineStyle(kDashDotted);
  oneP->Draw();

  TLine *twoM = new TLine(pTBins[0]-5, -2, pTBins[nBinspT], -2);
  twoM->SetLineColor(kBlack);
  twoM->SetLineStyle(kDashDotted);
  twoM->Draw();
  TLine *twoP = new TLine(pTBins[0]-5, 2, pTBins[nBinspT], 2);
  twoP->SetLineColor(kBlack);
  twoP->SetLineStyle(kDashDotted);
  twoP->Draw();

  c->SaveAs("lth_pulls_Y.pdf");
  c->Clear();

  c->Destructor();

}
