// code to plot the fit results

void plotRes()
{
  // get the histo limits
  TFile *fIn = new TFile("files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get A, lambda, chiProb values for each bin
  string lbl[] = {"Data", "NP", "PR", "J"};
  TFile *fInd = new TFile("files/finalFitRes.root");
  TGraphErrors **graph_A = new TGraphErrors*[4];
  TGraphErrors **graph_B = new TGraphErrors*[4];
  TGraph **graph_chi = new TGraph*[4];
  for(int i_t = 0; i_t < 4; i_t++) {
    graph_A[i_t] = (TGraphErrors*)fInd->Get(Form("graph_A_%s", lbl[i_t].c_str()));
    graph_B[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
    graph_chi[i_t] = (TGraph*)fInd->Get(Form("graph_chiP_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();

  TFile *fIndth = new TFile("../../Simult/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth = new TGraphErrors*[4];
  for(int i_t = 0; i_t < 4; i_t++) {
    graph_lth[i_t] = (TGraphErrors*)fIndth->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }
  fIndth->Close();

  // the results for lbd_phi are actually for B = 2lph / (3+lth)
  // resolve to get lbd_phi
  double lth_B[4][nBinspT], elth_B[4][nBinspT];
  double lth[4][nBinspT], elth[4][nBinspT];
  for(int i = 0; i < nBinspT; i++)
    for(int j = 0; j < 4; j++) {
      lth_B[j][i] = graph_B[j]->GetY()[i];
      elth_B[j][i] = graph_B[j]->GetEY()[i];

      lth[j][i] = graph_lth[j]->GetY()[i];
      elth[j][i] = graph_lth[j]->GetEY()[i];
    }

  double lph[4][nBinspT];
  double elph[4][nBinspT];
  for(int i = 0; i < nBinspT; i++)
    for(int j = 0; j < 4; j++) {
      lph[j][i] = lth_B[j][i]/2.*(3.+lth[j][i]);
      elph[j][i] = sqrt(pow((3.+lth[j][i])/2.,2)*elth_B[j][i]*elth_B[j][i] + pow(lth_B[j][i]/2.,2)*elth[j][i]*elth[j][i]);
    }
  TGraphErrors **graph_lph = new TGraphErrors*[4];
  for(int j = 0; j < 4; j++)
    graph_lph[j] = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), lph[j], graph_lth[0]->GetEX(), elph[j]);
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0]-5, -0.25, pTBins[nBinspT], 0.25);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#phi}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  //fl->SetTitle("Run 2 #lambda_{#phi} (PR)");
  fl->SetTitle("Run 2 #beta (PR)");

  int col[] = {kViolet, kRed, kBlack, kBlue};
  for(int i = 0; i < 4; i++) {
    //    graph_B[i]->SetLineStyle(kDashed);
    graph_B[i]->SetLineColor(col[i]);
    graph_B[i]->SetMarkerColor(col[i]);
    graph_B[i]->Draw("p same");

    /* graph_lph[i]->SetLineColor(col[i]);
    graph_lph[i]->SetMarkerColor(col[i]);
    graph_lph[i]->Draw("p same");*/
  }

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg = new TLegend(0.77, 0.7, 0.97, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(graph_B[0], "total", "pl");
  leg->AddEntry(graph_B[1], "NP contrib", "pl");
  leg->AddEntry(graph_B[2], "prompt", "pl");
  leg->AddEntry(graph_B[3], "prompt #psi(2S)", "pl");
  leg->Draw();
  
  c->SaveAs("plots/ratioFinal/par_lth.pdf");
  c->Clear();

  // draw just final lambda_th(pT)
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -0.25, pTBins[nBinspT], 0.25);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#lambda_{#phi}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  //fl2->SetTitle("Run 2 #lambda_{#phi} (prompt #psi(2S))");
  fl2->SetTitle("Run 2 #beta (prompt #psi(2S))");

  /*  graph_lph[3]->SetLineColor(kBlack);
  graph_lph[3]->SetMarkerColor(kBlack);
  graph_lph[3]->Draw("p same");*/

  //graph_B[3]->SetLineStyle(kDashed);
  graph_B[3]->SetLineColor(kBlue);
  graph_B[3]->SetMarkerColor(kBlue);
  graph_B[3]->Draw("p same");

  graph_B[1]->SetLineColor(kRed);
  graph_B[1]->SetMarkerColor(kRed);
  graph_B[1]->Draw("p same");

  zero->Draw();

  TF1 *fcon = new TF1("fc", "[0]", pTBins[0], pTBins[nBinspT]);
  fcon->SetParameter(0, -0.01);
  //graph_B[1]->Fit(fcon);
  //graph_B[3]->Fit(fcon);

  c->SaveAs("plots/ratioFinal/par_lth_F.pdf");
  c->Clear();

  // draw A(pT)
  c->SetLogy();
  TH1F *fa = c->DrawFrame(pTBins[0], 1e-2, pTBins[nBinspT], 6e-1);
  fa->SetXTitle("p_{T} (GeV)");
  fa->SetYTitle("A");
  fa->GetYaxis()->SetTitleOffset(1.3);
  fa->GetYaxis()->SetLabelOffset(0.01);
  fa->SetTitle("Run 2 A");

  // combine both lambda_th distributions
  for(int i = 0; i < 4; i++) {
    graph_A[i]->SetLineColor(col[i]);
    graph_A[i]->SetMarkerColor(col[i]);
    graph_A[i]->Draw("p same");
  }

  TLine *trans1_A = new TLine(46, 1e-2, 46, 6e-1);
  trans1_A->SetLineColor(kBlack);
  trans1_A->SetLineStyle(kDashed);
  trans1_A->Draw();
  TLine *trans2_A = new TLine(66, 1e-2, 66, 6e-1);
  trans2_A->SetLineColor(kBlack);
  trans2_A->SetLineStyle(kDashed);
  trans2_A->Draw();

  c->SaveAs("plots/ratioFinal/par_A.pdf");
  c->Clear();

  // draw chiProb(pT)
  c->SetLogy(0);
  TH1F *fc = c->DrawFrame(pTBins[0], 0, pTBins[nBinspT], 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("P(#chi^{2}, ndf)");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);
  fc->SetTitle("Run 2 P(#chi^{2}, ndf)");

  // combine both lambda_th distributions
  for(int i = 0; i < 4; i++) {
    graph_chi[i]->SetLineColor(col[i]);
    graph_chi[i]->SetMarkerColor(col[i]);
    graph_chi[i]->SetMarkerStyle(20);
    graph_chi[i]->SetMarkerSize(.75);
    graph_chi[i]->Draw("p same");
  }

  TLine *trans1_C = new TLine(46, 0, 46, 1);
  trans1_C->SetLineColor(kBlack);
  trans1_C->SetLineStyle(kDashed);
  trans1_C->Draw();
  TLine *trans2_C = new TLine(66, 0, 66, 1);
  trans2_C->SetLineColor(kBlack);
  trans2_C->SetLineStyle(kDashed);
  trans2_C->Draw();

  c->SaveAs("plots/ratioFinal/par_chiP.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();


}
