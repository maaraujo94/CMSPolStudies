// code to plot the fit results

void plotLth()
{
  // get the histo limits
  TFile *fIn = new TFile("../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get lambda values for each bin
  string lbl[] = {"Data", "NP", "PR", "J"};
  TFile *fInd = new TFile("../PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth = new TGraphErrors*[4];
  for(int i_t = 0; i_t < 4; i_t++) {
    graph_lth[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);

  // draw lambda_th(pT) - just peak
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("2017 #lambda_{#theta}");

  graph_lth[0]->SetLineColor(kBlack);
  graph_lth[0]->SetMarkerColor(kBlack);
  graph_lth[0]->Draw("p same");

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  
  c->SaveAs("plots/ratioFinal/lth1.pdf");
  c->Clear();

  // draw lambda_th(pT) - peak + np + pr
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("2017 #lambda_{#theta}");

  graph_lth[0]->SetLineColor(kBlack);
  graph_lth[0]->SetMarkerColor(kBlack);
  graph_lth[0]->Draw("p same");

  graph_lth[2]->SetLineColor(kViolet);
  graph_lth[2]->SetMarkerColor(kViolet);
  graph_lth[2]->Draw("p same");

  zero->Draw();
  
  c->SaveAs("plots/ratioFinal/lth2.pdf");

  graph_lth[1]->SetLineColor(kBlue);
  graph_lth[1]->SetMarkerColor(kBlue);
  graph_lth[1]->Draw("p same");
  
  c->SaveAs("plots/ratioFinal/lth2_np.pdf");
  c->Clear();

  // draw lambda_th(pT) - peak + pr + jpsi
  TH1F *fl3 = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl3->SetXTitle("p_{T} (GeV)");
  fl3->SetYTitle("#lambda_{#theta}");
  fl3->GetYaxis()->SetTitleOffset(1.3);
  fl3->GetYaxis()->SetLabelOffset(0.01);
  fl3->SetTitle("2017 #lambda_{#theta}");

  graph_lth[0]->SetLineColor(kBlack);
  graph_lth[0]->SetMarkerColor(kBlack);
  graph_lth[0]->Draw("p same");

  graph_lth[2]->SetLineColor(kViolet);
  graph_lth[2]->SetMarkerColor(kViolet);
  graph_lth[2]->Draw("p same");
  
  graph_lth[3]->SetLineColor(kRed);
  graph_lth[3]->SetMarkerColor(kRed);
  graph_lth[3]->Draw("p same");

  zero->Draw();
  
  c->SaveAs("plots/ratioFinal/lth3.pdf");
  c->Clear();

  c->Destructor();
  
  fIn->Close();


}
