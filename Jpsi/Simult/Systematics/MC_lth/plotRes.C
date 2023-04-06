// code to plot the fit results

void plotRes()
{
  // get the histo limits
  TFile *fIn = new TFile("files/histoStore.root");
  TH2D* rHist;
  fIn->GetObject("MCH", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get A, lambda, chiProb values for each bin
  string lbl[] = {"MCp", "MCm"};
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
  c->SetTopMargin(0.015);
  c->SetLeftMargin(0.1);
  
  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.4);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->GetXaxis()->SetLabelOffset(0.015);
  fl->GetXaxis()->SetTitleOffset(1.3);
  fl->SetTitle("");

  int col[] = {kBlue, kRed};
  int st[] = {20, 25};
  for(int i = 0; i < 2; i++) {
    graph_lth[i]->SetLineColor(col[i]);
    graph_lth[i]->SetMarkerColor(col[i]);
    graph_lth[i]->SetMarkerStyle(st[i]);
    graph_lth[i]->SetMarkerSize(.75);
    graph_lth[i]->Draw("p same");
  }

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  leg->AddEntry(graph_lth[0], "MC (#lambda = +0.4)", "pl");
  leg->AddEntry(graph_lth[1], "MC (#lambda = -0.1)", "pl");
  leg->Draw();
  
  c->SaveAs("plots/ratioFinal/par_lth.pdf");
  c->Clear();

  // also draw the deviation from expected (+0.4, -0.1)
  double dv[2][nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    dv[0][i] = graph_lth[0]->GetY()[i]-0.4;
    dv[1][i] = graph_lth[1]->GetY()[i]+0.1;
  }
  TGraphErrors **graph_dlth = new TGraphErrors*[2];
  for(int i = 0; i < 2; i++) {
    graph_dlth[i] = new TGraphErrors(nBinspT, graph_lth[i]->GetX(), dv[i], graph_lth[i]->GetEX(), graph_lth[i]->GetEY());
  }
  
  TH1F *fdl = c->DrawFrame(pTBins[0]-5, -0.0299, pTBins[nBinspT], 0.0299);
  fdl->SetXTitle("p_{T} (GeV)");
  fdl->SetYTitle("#Delta#lambda_{#theta}");
  fdl->GetYaxis()->SetTitleOffset(1.4);
  fdl->GetYaxis()->SetLabelOffset(0.01);
  fdl->GetXaxis()->SetLabelOffset(0.015);
  fdl->GetXaxis()->SetTitleOffset(1.3);
  fdl->SetTitle("");

  for(int i = 0; i < 2; i++) {
    graph_dlth[i]->SetLineColor(col[i]);
    graph_dlth[i]->SetMarkerColor(col[i]);
    graph_dlth[i]->SetMarkerStyle(st[i]);
    graph_dlth[i]->SetMarkerSize(.75);
    graph_dlth[i]->Draw("p same");
  }

  zero->Draw();

  leg->Draw();
  
  c->SaveAs("plots/ratioFinal/par_dlth.pdf");
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
  for(int i = 0; i < 2; i++) {
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
  fc->SetTitle("Run 2 P(#chi^{2}, ndf)");

  // combine both lambda_th distributions
  for(int i = 0; i < 2; i++) {
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
