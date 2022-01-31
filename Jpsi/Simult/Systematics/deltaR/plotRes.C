// code to plot the fit results

void plotRes()
{
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get lambda values for each eff model
  TGraphErrors **graph_lth = new TGraphErrors*[2];
  TFile *fInd1 = new TFile("../../../Simult_dR1/PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fInd1->Get(Form("graph_lambda_J"));
  fInd1->Close();
  TFile *fInd2 = new TFile("../../../Simult_dR2/PR_fit/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fInd2->Get(Form("graph_lambda_J"));
  fInd2->Close();
 
  // get the final lth from the base SB/MC fit
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  TGraphErrors *graph_lthBase = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  fIndB->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);

  // draw just final lambda_th(pT) - comp btw std, alt
  double val[2][nBinspT];
  for(int j = 0; j < 2; j++) {
    for(int i = 0; i < nBinspT; i++) { 
      val[j][i] = graph_lth[j]->GetY()[i] - graph_lthBase->GetY()[i];
    }
  }
  TGraphErrors *g_lthD1 = new TGraphErrors(nBinspT, graph_lthBase->GetX(), val[0], graph_lthBase->GetEX(), graph_lth[0]->GetEY());
  TGraphErrors *g_lthD2 = new TGraphErrors(nBinspT, graph_lthBase->GetX(), val[1], graph_lthBase->GetEX(), graph_lth[0]->GetEY());

  double d_lim = 0.6;
  
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -d_lim, pTBins[nBinspT], d_lim);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#delta#lambda_{#theta} (dev - base)");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("Run 2 #delta#lambda_{#theta} (prompt J/#psi)");
  
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
  TLine *trans1D = new TLine(46, -d_lim, 46, d_lim);
  trans1D->SetLineColor(kBlack);
  trans1D->SetLineStyle(kDashed);
  trans1D->Draw();
  TLine *trans2D = new TLine(66, -d_lim, 66, d_lim);
  trans2D->SetLineColor(kBlack);
  trans2D->SetLineStyle(kDashed);
  trans2D->Draw();

  TLegend *leg = new TLegend(0.65, 0.15, 0.9, 0.3);
  leg->SetTextSize(0.03);
  leg->AddEntry(g_lthD1, "#deltaR > 0.17", "pl");
  leg->AddEntry(g_lthD2, "#deltaR > 0.15", "pl");
  leg->Draw();
  
  c->SaveAs("par_lth_F.pdf");
  c->Clear();

  // now compare absolute values
  double a_lim = 1;
  
  TH1F *fl3 = c->DrawFrame(pTBins[0]-5, -a_lim, pTBins[nBinspT], a_lim);
  fl3->SetXTitle("p_{T} (GeV)");
  fl3->SetYTitle("#lambda_{#theta}");
  fl3->GetYaxis()->SetTitleOffset(1.3);
  fl3->GetYaxis()->SetLabelOffset(0.01);
  fl3->SetTitle("Run 2 #lambda_{#theta} (prompt J/#psi)");
  
  graph_lth[0]->SetLineColor(kBlue);
  graph_lth[0]->SetMarkerColor(kBlue);
  graph_lth[0]->SetMarkerStyle(20);
  graph_lth[0]->SetMarkerSize(.5);
  graph_lth[0]->Draw("p same");
  
  graph_lth[1]->SetLineColor(kRed);
  graph_lth[1]->SetMarkerColor(kRed);
  graph_lth[1]->SetMarkerStyle(20);
  graph_lth[1]->SetMarkerSize(.5);
  graph_lth[1]->Draw("p same");

  graph_lthBase->SetLineColor(kBlack);
  graph_lthBase->SetMarkerColor(kBlack);
  graph_lthBase->SetMarkerStyle(20);
  graph_lthBase->SetMarkerSize(.5);
  graph_lthBase->Draw("p same");
  
  zero->Draw();
  TLine *trans1A = new TLine(46, -a_lim, 46, a_lim);
  trans1A->SetLineColor(kBlack);
  trans1A->SetLineStyle(kDashed);
  trans1A->Draw();
  TLine *trans2A = new TLine(66, -a_lim, 66, a_lim);
  trans2A->SetLineColor(kBlack);
  trans2A->SetLineStyle(kDashed);
  trans2A->Draw();

  leg->AddEntry(graph_lthBase, "no cut", "pl");
  leg->Draw();
  
  c->SaveAs("par_lth_abs.pdf");
  c->Clear();

  c->Destructor();
  
  fIn->Close();


}
