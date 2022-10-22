// code to plot the fit results

void plotLth_NP()
{
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();

  // get the histo limits - coarse version
  TFile *fIn_c = new TFile("files/chistStore.root");
  TH2D* rHist_c;
  fIn_c->GetObject("cHistNPB", rHist_c);
  
  int nBinspT_c = rHist_c->GetNbinsY();
  const double *pTBins_c = rHist_c->GetYaxis()->GetXbins()->GetArray();

  // get the fit results
  // get lambda values for each eff model
  TGraphErrors **graph_lth = new TGraphErrors*[3];
  TGraphErrors **graph_lth_c = new TGraphErrors*[3];
  TFile *fInd = new TFile("files/finalFitRes_NP.root");
  string lbl[] = {"B", "T", "L"};
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lth[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
    graph_lth_c[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_c_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);

  // draw just final lambda_th(pT) - comp btw std, alt
  double val[2][nBinspT], unc[2][nBinspT];
  for(int j = 0; j < 2; j++) {
    for(int i = 0; i < nBinspT; i++) { 
      val[j][i] = graph_lth[j+1]->GetY()[i] - graph_lth[0]->GetY()[i];
      double unc1 = graph_lth[j+1]->GetEY()[i];
      double unc2 = graph_lth[0]->GetEY()[i];
      unc[j][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));
    }
  }
  TGraphErrors *g_lthD1 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), val[0], graph_lth[0]->GetEX(), unc[0]);
  TGraphErrors *g_lthD2 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), val[1], graph_lth[0]->GetEX(), unc[1]);

  // diffs for coarse bins
  double val_c[2][nBinspT_c], unc_c[2][nBinspT_c];
  for(int j = 0; j < 2; j++) {
    for(int i = 0; i < nBinspT_c; i++) { 
      val_c[j][i] = graph_lth_c[j+1]->GetY()[i] - graph_lth_c[0]->GetY()[i];
      double unc1 = graph_lth_c[j+1]->GetEY()[i];
      double unc2 = graph_lth_c[0]->GetEY()[i];
      unc_c[j][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));
    }
  }
  TGraphErrors *g_lthD1_c = new TGraphErrors(nBinspT_c, graph_lth_c[0]->GetX(), val_c[0], graph_lth_c[0]->GetEX(), unc_c[0]);
  TGraphErrors *g_lthD2_c = new TGraphErrors(nBinspT_c, graph_lth_c[0]->GetX(), val_c[1], graph_lth_c[0]->GetEX(), unc_c[1]);

  double d_lim = 0.4;

  // only plot the fiducial cut plots in specific ranges
  // [0]: dR > 0.17 -> 25 to 70 GeV
  // [1]: dR > 0.15 -> 70 to 120 GeV
  double* binsX = graph_lth[0]->GetX();
  int cut_val;
  for(int i = 1; i < nBinspT; i++) {
    if(binsX[i] > 70 && binsX[i-1] < 70) cut_val = i;
  }

  graph_lth[0]->RemovePoint(8);
  for(int i = cut_val; i < nBinspT; i++){
    g_lthD1->RemovePoint(cut_val);
    graph_lth[1]->RemovePoint(cut_val);
  }
  g_lthD1->RemovePoint(8);
  graph_lth[1]->RemovePoint(8);
  for(int i = 0; i < cut_val; i++) {
    g_lthD2->RemovePoint(0);
    graph_lth[2]->RemovePoint(0);
  }

  //same thing for coarse bins
  graph_lth_c[0]->RemovePoint(1);
  for(int i = 3; i < nBinspT_c; i++){
    g_lthD1_c->RemovePoint(3);
    graph_lth_c[1]->RemovePoint(3);
  }
  g_lthD1_c->RemovePoint(1);
  graph_lth_c[1]->RemovePoint(1);
  for(int i = 0; i < 3; i++) {
    g_lthD2_c->RemovePoint(0);
    graph_lth_c[2]->RemovePoint(0);
  }

  // draw the differences - fine bins
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -d_lim, pTBins[nBinspT], d_lim);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#delta#lambda_{#theta} (dev - base)");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("Run 2 #delta#lambda_{#theta} (NP J/#psi)");
  
  g_lthD1->SetLineColor(kBlue);
  g_lthD1->SetMarkerColor(kBlue);
  g_lthD1->SetMarkerStyle(20);
  g_lthD1->SetMarkerSize(.5);
  g_lthD1->Draw("p same");

  g_lthD2->SetLineColor(kRed);
  g_lthD2->SetMarkerColor(kRed);
  g_lthD2->SetMarkerStyle(20);
  g_lthD2->SetMarkerSize(.5);
  g_lthD2->Draw("p same");

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  TLine *trans1D = new TLine(45, -d_lim, 45, d_lim);
  trans1D->SetLineColor(kBlack);
  trans1D->SetLineStyle(kDashed);
  //trans1D->Draw();
  TLine *trans2D = new TLine(70, -d_lim, 70, d_lim);
  trans2D->SetLineColor(kBlack);
  trans2D->SetLineStyle(kDashed);
  trans2D->Draw();
  
  TLegend *leg = new TLegend(0.65, 0.15, 0.9, 0.3);
  leg->SetTextSize(0.03);
  leg->AddEntry(g_lthD1, "#DeltaR > 0.17", "pl");
  leg->AddEntry(g_lthD2, "#DeltaR > 0.15", "pl");
  leg->Draw();
  
  c->SaveAs("plots/par_dlth_NP.pdf");
  c->Clear();

  // draw the differences - coarse bins
  TH1F *fl2_c = c->DrawFrame(pTBins[0]-5, -d_lim, pTBins[nBinspT], d_lim);
  fl2_c->SetXTitle("p_{T} (GeV)");
  fl2_c->SetYTitle("#delta#lambda_{#theta} (dev - base)");
  fl2_c->GetYaxis()->SetTitleOffset(1.3);
  fl2_c->GetYaxis()->SetLabelOffset(0.01);
  fl2_c->SetTitle("Run 2 #delta#lambda_{#theta} (NP J/#psi)");
  
  g_lthD1_c->SetLineColor(kBlue);
  g_lthD1_c->SetMarkerColor(kBlue);
  g_lthD1_c->SetMarkerStyle(20);
  g_lthD1_c->SetMarkerSize(.5);
  g_lthD1_c->Draw("p same");

  g_lthD2_c->SetLineColor(kRed);
  g_lthD2_c->SetMarkerColor(kRed);
  g_lthD2_c->SetMarkerStyle(20);
  g_lthD2_c->SetMarkerSize(.5);
  g_lthD2_c->Draw("p same");

  zero->Draw();
  trans2D->Draw();
  
  leg->Draw();
  
  c->SaveAs("plots/par_c_dlth_NP.pdf");
  c->Clear();

  // now compare absolute values
  double a_lim = 1;
  
  TH1F *fl3 = c->DrawFrame(pTBins[0]-5, -a_lim, pTBins[nBinspT], a_lim);
  fl3->SetXTitle("p_{T} (GeV)");
  fl3->SetYTitle("#lambda_{#theta}");
  fl3->GetYaxis()->SetTitleOffset(1.3);
  fl3->GetYaxis()->SetLabelOffset(0.01);
  fl3->SetTitle("Run 2 #lambda_{#theta} (NP J/#psi)");
    
  graph_lth[1]->SetLineColor(kBlue);
  graph_lth[1]->SetMarkerColor(kBlue);
  graph_lth[1]->SetMarkerStyle(20);
  graph_lth[1]->SetMarkerSize(.5);
  graph_lth[1]->Draw("p same");

  graph_lth[2]->SetLineColor(kRed);
  graph_lth[2]->SetMarkerColor(kRed);
  graph_lth[2]->SetMarkerStyle(20);
  graph_lth[2]->SetMarkerSize(.5);
  graph_lth[2]->Draw("p same");

  graph_lth[0]->SetLineColor(kBlack);
  graph_lth[0]->SetMarkerColor(kBlack);
  graph_lth[0]->SetMarkerStyle(20);
  graph_lth[0]->SetMarkerSize(.5);
  graph_lth[0]->Draw("p same");
  
  zero->Draw();
  TLine *trans1A = new TLine(46, -a_lim, 46, a_lim);
  trans1A->SetLineColor(kBlack);
  trans1A->SetLineStyle(kDashed);
  //trans1A->Draw();
  TLine *trans2A = new TLine(70, -a_lim, 70, a_lim);
  trans2A->SetLineColor(kBlack);
  trans2A->SetLineStyle(kDashed);
  trans2A->Draw();

  leg->AddEntry(graph_lth[0], "no cut", "pl");
  leg->Draw();
  
  c->SaveAs("plots/par_lth_NP.pdf");
  c->Clear();

  // same with coarse bins
  TH1F *fl3_c = c->DrawFrame(pTBins[0]-5, -a_lim, pTBins[nBinspT], a_lim);
  fl3_c->SetXTitle("p_{T} (GeV)");
  fl3_c->SetYTitle("#lambda_{#theta}");
  fl3_c->GetYaxis()->SetTitleOffset(1.3);
  fl3_c->GetYaxis()->SetLabelOffset(0.01);
  fl3_c->SetTitle("Run 2 #lambda_{#theta} (NP J/#psi)");
    
  graph_lth_c[1]->SetLineColor(kBlue);
  graph_lth_c[1]->SetMarkerColor(kBlue);
  graph_lth_c[1]->SetMarkerStyle(20);
  graph_lth_c[1]->SetMarkerSize(.5);
  graph_lth_c[1]->Draw("p same");

  graph_lth_c[2]->SetLineColor(kRed);
  graph_lth_c[2]->SetMarkerColor(kRed);
  graph_lth_c[2]->SetMarkerStyle(20);
  graph_lth_c[2]->SetMarkerSize(.5);
  graph_lth_c[2]->Draw("p same");

  graph_lth_c[0]->SetLineColor(kBlack);
  graph_lth_c[0]->SetMarkerColor(kBlack);
  graph_lth_c[0]->SetMarkerStyle(20);
  graph_lth_c[0]->SetMarkerSize(.5);
  graph_lth_c[0]->Draw("p same");
  
  zero->Draw();
  trans2A->Draw();

  leg->Draw();
  
  c->SaveAs("plots/par_c_lth_NP.pdf");
  c->Clear();

  c->Destructor();
  
  fIn->Close();
  fIn_c->Close();


}
