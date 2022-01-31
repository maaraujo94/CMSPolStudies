// code to compare the deviations for several scenarios
// 1) 2017 - 2018 / Run2
// 2) f(pT) - 1/f(pT) / Run2
// 3) pT_cut - Run2 / Run2
// 4) eta_cut - Run2 / Run2

void plotRes()
{
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - 6 sets
  TGraphErrors **graph_lth = new TGraphErrors*[7];
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

  // get the differences and the unc band
  double rdiff[4][nBinspT], diff[4][nBinspT], za[nBinspT], rerr[nBinspT];
  for(int i = 0; i < nBinspT; i++) { 
    diff[0][i] = (graph_lth[1]->GetY()[i] - graph_lth[2]->GetY()[i]);
    diff[1][i] = (graph_lth[3]->GetY()[i] - graph_lth[4]->GetY()[i]);
    diff[2][i] = (graph_lth[5]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[3][i] = (graph_lth[6]->GetY()[i] - graph_lth[0]->GetY()[i]);
    cout << i; 
    for(int j = 0; j < 4; j++){
      rdiff[j][i] = diff[j][i]/graph_lth[0]->GetY()[i];
      rdiff[j][i] *= 100.;
      cout << " " << diff[j][i];
    }
    cout << endl;
    za[i] = 0;
    rerr[i] = graph_lth[0]->GetEY()[i]/graph_lth[0]->GetY()[i]*100;
  }
  TGraphErrors *g_lthY = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[0], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthEff = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[1], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthpT = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[2], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthEta = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[3], graph_lth[0]->GetEX(), za);

  TGraphErrors *gr_lthY = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), rdiff[0], graph_lth[0]->GetEX(), za);
  TGraphErrors *gr_lthEff = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), rdiff[1], graph_lth[0]->GetEX(), za);
  TGraphErrors *gr_lthpT = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), rdiff[2], graph_lth[0]->GetEX(), za);
  TGraphErrors *gr_lthEta = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), rdiff[3], graph_lth[0]->GetEX(), za);

  TGraphErrors *g_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_lth[0]->GetEY());
  TGraphErrors *gr_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), rerr);
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  double d_lim = 60;

  // FIRST - draw the relative difference in %
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -d_lim, pTBins[nBinspT], d_lim);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("rel #delta#lambda_{#theta} (%)");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("Relative difference #delta#lambda_{#theta} (prompt J/#psi)");
  
  gr_lthY->SetLineColor(kBlack);
  gr_lthY->SetMarkerColor(kBlack);
  gr_lthY->SetMarkerStyle(20);
  gr_lthY->SetMarkerSize(.5);
  gr_lthY->Draw("pl same");

  gr_lthEff->SetLineColor(kBlue);
  gr_lthEff->SetMarkerColor(kBlue);
  gr_lthEff->SetMarkerStyle(20);
  gr_lthEff->SetMarkerSize(.5);
  gr_lthEff->Draw("pl same");

  gr_lthpT->SetLineColor(kGreen+1);
  gr_lthpT->SetMarkerColor(kGreen+1);
  gr_lthpT->SetMarkerStyle(20);
  gr_lthpT->SetMarkerSize(.5);
  gr_lthpT->Draw("pl same");

  gr_lthEta->SetLineColor(kViolet-1);
  gr_lthEta->SetMarkerColor(kViolet-1);
  gr_lthEta->SetMarkerStyle(20);
  gr_lthEta->SetMarkerSize(.5);
  gr_lthEta->Draw("pl same");

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  TLine *trans1D = new TLine(46, -d_lim, 46, d_lim);
  trans1D->SetLineColor(kBlack);
  trans1D->SetLineStyle(kDashed);
  trans1D->Draw();
  TLine *trans2D = new TLine(66, -d_lim, 66, d_lim);
  trans2D->SetLineColor(kBlack);
  trans2D->SetLineStyle(kDashed);
  trans2D->Draw();

  TLegend *leg = new TLegend(0.675, 0.75, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(gr_lthY, "2017 - 2018", "pl");
  leg->AddEntry(gr_lthEff, "f(p_{T}) - 1/f(p_{T})", "pl");
  leg->AddEntry(gr_lthpT, "p_{T} cut - base", "pl");
  leg->AddEntry(gr_lthEta, "#eta cut - base", "pl");
  leg->Draw();
  
  c->SaveAs("lth_relDiff.pdf");

  gr_unc->SetLineColor(kRed);
  gr_unc->SetFillColorAlpha(kRed, 0.4);
  gr_unc->Draw("ce3");

  c->SaveAs("lth_relDiff_band.pdf");
  c->Clear();

  // SECOND - draw the relative difference in %
  double da_lim = 0.2;
  
  TH1F *fl3 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl3->SetXTitle("p_{T} (GeV)");
  fl3->SetYTitle("#delta#lambda_{#theta}");
  fl3->GetYaxis()->SetTitleOffset(1.3);
  fl3->GetYaxis()->SetLabelOffset(0.01);
  fl3->SetTitle("Difference #delta#lambda_{#theta} (prompt J/#psi)");
  
  g_lthY->SetLineColor(kBlack);
  g_lthY->SetMarkerColor(kBlack);
  g_lthY->SetMarkerStyle(20);
  g_lthY->SetMarkerSize(.5);
  g_lthY->Draw("pl same");

  g_lthEff->SetLineColor(kBlue);
  g_lthEff->SetMarkerColor(kBlue);
  g_lthEff->SetMarkerStyle(20);
  g_lthEff->SetMarkerSize(.5);
  g_lthEff->Draw("pl same");

  g_lthpT->SetLineColor(kGreen+1);
  g_lthpT->SetMarkerColor(kGreen+1);
  g_lthpT->SetMarkerStyle(20);
  g_lthpT->SetMarkerSize(.5);
  g_lthpT->Draw("pl same");

  g_lthEta->SetLineColor(kViolet-1);
  g_lthEta->SetMarkerColor(kViolet-1);
  g_lthEta->SetMarkerStyle(20);
  g_lthEta->SetMarkerSize(.5);
  g_lthEta->Draw("pl same");

  g_unc->SetLineColor(kRed);
  g_unc->SetFillColorAlpha(kRed, 0.4);
  g_unc->Draw("ce3");

  zero->Draw();
  TLine *trans1A = new TLine(46, -da_lim, 46, da_lim);
  trans1A->SetLineColor(kBlack);
  trans1A->SetLineStyle(kDashed);
  trans1A->Draw();
  TLine *trans2A = new TLine(66, -da_lim, 66, da_lim);
  trans2A->SetLineColor(kBlack);
  trans2A->SetLineStyle(kDashed);
  trans2A->Draw();

  /*  TLegend *leg = new TLegend(0.675, 0.75, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(g_lthY, "2017-2018", "pl");
  leg->AddEntry(g_lthEff, "f(p_{T})-1/f(p_{T})", "pl");
  leg->AddEntry(g_lthpT, "p_{T} cut - base", "pl");
  leg->AddEntry(g_lthEta, "#eta cut - base", "pl");*/
  leg->Draw();
  
  c->SaveAs("lth_absDiff_band.pdf");
  c->Clear();

  c->Destructor();
  
  fIn->Close();


}
