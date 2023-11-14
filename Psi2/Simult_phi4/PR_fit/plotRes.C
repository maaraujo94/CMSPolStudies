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
  TGraphErrors **graph_lth = new TGraphErrors*[4];
  TGraph **graph_chi = new TGraph*[4];
  for(int i_t = 0; i_t < 4; i_t++) {
    graph_A[i_t] = (TGraphErrors*)fInd->Get(Form("graph_A_%s", lbl[i_t].c_str()));
    graph_lth[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
    graph_chi[i_t] = (TGraph*)fInd->Get(Form("graph_chiP_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("");

  int col[] = {kViolet, kRed, kBlack, kBlue};
  for(int i = 0; i < 4; i++) {
    graph_lth[i]->SetLineColor(col[i]);
    graph_lth[i]->SetMarkerColor(col[i]);
    graph_lth[i]->Draw("p same");
  }

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg = new TLegend(0.65, 0.12, 0.95, 0.32);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  leg->AddEntry(graph_lth[0], "total", "pl");
  leg->AddEntry(graph_lth[1], "non-prompt #psi(2S)", "pl");
  leg->AddEntry(graph_lth[2], "prompt", "pl");
  leg->AddEntry(graph_lth[3], "prompt #psi(2S)", "pl");
  leg->Draw();
  
  c->SaveAs("plots/ratioFinal/par_lth.pdf");
  c->Clear();

  // draw just final lambda_th(pT)
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#lambda_{#theta}^{PR}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("");

  graph_lth[3]->SetLineColor(kBlack);
  graph_lth[3]->SetMarkerColor(kBlack);
  graph_lth[3]->Draw("p same");

  zero->Draw();
  //trans1->Draw();
  
  c->SaveAs("plots/ratioFinal/par_lth_F.pdf");
  c->Clear();

  // draw A(pT)
  c->SetLogy();
  TH1F *fa = c->DrawFrame(pTBins[0], 1e-2, pTBins[nBinspT], 6e-1);
  fa->SetXTitle("p_{T} (GeV)");
  fa->SetYTitle("A");
  fa->GetYaxis()->SetTitleOffset(1.3);
  fa->GetYaxis()->SetLabelOffset(0.01);
  fa->SetTitle("");

  // combine both lambda_th distributions
  for(int i = 0; i < 4; i++) {
    graph_A[i]->SetLineColor(col[i]);
    graph_A[i]->SetMarkerColor(col[i]);
    graph_A[i]->Draw("p same");
  }

  c->SaveAs("plots/ratioFinal/par_A.pdf");
  c->Clear();

  // draw chiProb(pT)
  c->SetLogy(0);
  TH1F *fc = c->DrawFrame(pTBins[0], 0, pTBins[nBinspT], 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("P(#chi^{2}, ndf)");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);
  fc->SetTitle("");

  // combine both lambda_th distributions
  for(int i = 0; i < 4; i++) {
    graph_chi[i]->SetLineColor(col[i]);
    graph_chi[i]->SetMarkerColor(col[i]);
    graph_chi[i]->SetMarkerStyle(20);
    graph_chi[i]->SetMarkerSize(.75);
    graph_chi[i]->Draw("p same");
  }

  c->SaveAs("plots/ratioFinal/par_chiP.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();


}
