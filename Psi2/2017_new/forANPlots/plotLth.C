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
  int cols[] = {kViolet-1, kRed, kBlack, kBlue, kGreen};

  // draw lambda_th(pT) - just peak
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("2017 #lambda_{#theta}");

  graph_lth[0]->SetLineColor(cols[0]);
  graph_lth[0]->SetMarkerColor(cols[0]);
  graph_lth[0]->Draw("p same");

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  
  c->SaveAs("plots/ratioFinal/lth1_2017.pdf");

  // add pr lambda_th
  graph_lth[2]->SetLineColor(cols[2]);
  graph_lth[2]->SetMarkerColor(cols[2]);
  graph_lth[2]->Draw("p same");

  c->SaveAs("plots/ratioFinal/lth2_2017.pdf");
  
  // add prompt jpsi lambda_th
  graph_lth[3]->SetLineColor(cols[3]);
  graph_lth[3]->SetMarkerColor(cols[3]);
  graph_lth[3]->Draw("p same");
  
  c->SaveAs("plots/ratioFinal/lth3_2017.pdf");
  c->Destructor();
  
  fIn->Close();


}
