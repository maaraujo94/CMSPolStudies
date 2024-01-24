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
  string lbl[] = {"Data", "NP", "J"};
  TFile *fInd = new TFile("files/finalFitRes.root");
  TGraphErrors **graph_A = new TGraphErrors*[3];
  TGraphErrors **graph_B = new TGraphErrors*[3];
  TGraph **graph_chi = new TGraph*[3];
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_A[i_t] = (TGraphErrors*)fInd->Get(Form("graph_A_%s", lbl[i_t].c_str()));
    graph_B[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
    graph_chi[i_t] = (TGraph*)fInd->Get(Form("graph_chiP_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();

  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.11);
  c->SetTopMargin(0.015);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0]-5, -0.25, pTBins[nBinspT], 0.25);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#beta");
  fl->GetYaxis()->SetTitleOffset(1.6);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("");

  int col[] = {kViolet, kRed, kBlue};
  for(int i = 0; i < 3; i++) {
    graph_B[i]->SetLineColor(col[i]);
    graph_B[i]->SetMarkerColor(col[i]);
    graph_B[i]->Draw("p same");
  }

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg = new TLegend(0.65, 0.75, 0.95, 0.95);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  leg->AddEntry(graph_B[0], "total", "pl");
  leg->AddEntry(graph_B[1], "non-prompt #psi(2S)", "pl");
  leg->AddEntry(graph_B[2], "prompt #psi(2S)", "pl");
  leg->Draw();
  
  c->SaveAs("plots/ratioFinal/par_lth.pdf");
  c->Clear();
  
  // draw just final lambda_th(pT)
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -0.25, pTBins[nBinspT], 0.25);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#beta");
  fl2->GetYaxis()->SetTitleOffset(1.6);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("");

  graph_B[2]->SetLineColor(kBlue);
  graph_B[2]->SetMarkerColor(kBlue);
  graph_B[2]->Draw("p same");

  graph_B[1]->SetLineColor(kRed);
  graph_B[1]->SetMarkerColor(kRed);
  graph_B[1]->Draw("p same");

  zero->Draw();

  TLegend *leg2 = new TLegend(0.65, 0.85, 0.95, 0.95);
  leg2->SetTextSize(0.03);
  leg2->SetBorderSize(0);
  leg2->SetFillColorAlpha(kWhite,0);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(graph_B[2], "prompt #psi(2S)", "pl");
  leg2->AddEntry(graph_B[1], "non-prompt #psi(2S)", "pl");
  leg2->Draw();

  c->SaveAs("plots/ratioFinal/par_lth_F.pdf");

  // draw the bands around dists
  int n = nBinspT+1;
  double xv[n], xe[n], yv1[n], ye1[n], yv2[n], ye2[n];
  for(int i = 0; i < n; i++) {
    if (i == 0) {
      xv[i] = pTBins[0];
      xe[i] = 0;
    }
    else if (i == n-1) {
      xv[i] = pTBins[n-1];
      xe[i] = 0;
    }
    else {
      xv[i] = 0.5*(pTBins[i-1]+pTBins[i]);
      xe[i] = 0.5*(pTBins[i]-pTBins[i-1]);
    }
    yv1[i] = -0.015;
    ye1[i] = 0.005;
    yv2[i] = 0.015;
    ye2[i] = 0.01;
  }

  TGraphErrors *g_1 = new TGraphErrors(n, xv, yv1, xe, ye1);
  g_1->SetFillColorAlpha(kBlue, 0.3);
  g_1->Draw("e3");
  TGraphErrors *g_2 = new TGraphErrors(n, xv, yv2, xe, ye2);
  g_2->SetFillColorAlpha(kRed, 0.3);
  g_2->Draw("e3");

  c->SaveAs("plots/ratioFinal/par_lth_band.pdf");
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
  for(int i = 0; i < 3; i++) {
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
  for(int i = 0; i < 3; i++) {
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
