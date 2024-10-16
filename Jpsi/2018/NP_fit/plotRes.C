// code to plot the fit results

void plotRes()
{
  // get the histo limits
  TFile *fIn = new TFile("files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_NP", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get A, lambda, chiProb values for each bin
  string lbl[] = {"NP", "NPc"};
  TFile *fInd = new TFile("files/finalFitRes.root");
  TGraphErrors **graph_A = new TGraphErrors*[2];
  TGraphErrors **graph_lth = new TGraphErrors*[2];
  TGraph **graph_chi = new TGraph*[2];
  for(int i_t = 0; i_t < 2; i_t++) {
    graph_A[i_t] = (TGraphErrors*)fInd->Get(Form("graph_A_%s", lbl[i_t].c_str()));
    graph_lth[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
    graph_chi[i_t] = (TGraph*)fInd->Get(Form("graph_chiP_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.02);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("");
    
  int col[] = {kRed+3, kRed};
  for(int i = 0; i < 2; i++) {
    graph_lth[i]->SetLineColor(col[i]);
    graph_lth[i]->SetMarkerColor(col[i]);
    graph_lth[i]->Draw("p same");
  }

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg = new TLegend(0.65, 0.85, 0.95, 0.95);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite, 0);
  leg->AddEntry(graph_lth[0], "NPS", "pl");
  leg->AddEntry(graph_lth[1], "non-prompt J/#psi", "pl");
  leg->Draw();
  
  c->SaveAs("plots/ratioFinal/par_lth.pdf");
  c->Clear();

  // draw just final lambda_th(pT)
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("");

  graph_lth[1]->SetLineColor(kBlack);
  graph_lth[1]->SetMarkerColor(kBlack);
  graph_lth[1]->Draw("p same");


  zero->Draw();
  
  c->SaveAs("plots/ratioFinal/par_lth_F.pdf");
  c->Clear();

  // draw A(pT)
  c->SetTopMargin(0.1);
    
  c->SetLogy();
  TH1F *fa = c->DrawFrame(pTBins[0], 1e-2, pTBins[nBinspT], 6e-1);
  fa->SetXTitle("p_{T} (GeV)");
  fa->SetYTitle("A");
  fa->GetYaxis()->SetTitleOffset(1.3);
  fa->GetYaxis()->SetLabelOffset(0.01);
  fa->SetTitle("2018 A");

  // combine both lambda_th distributions
  for(int i = 0; i < 2; i++) {
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
  fc->SetTitle("2018 P(#chi^{2}, ndf)");

  // combine both lambda_th distributions
  for(int i = 0; i < 2; i++) {
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
