// code to plot the fit results
void plotRes()
{
  // get the histo limits
  TFile *fIn = new TFile("files/chistStore.root");
  TH2D* rHist;
  fIn->GetObject("cHistB", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get A, lambda, chiProb values for each bin
  string lbl[] = {"B", "L", "T"};
  TFile *fInd = new TFile("files/finalFitRes.root");
  TGraphErrors **graph_lth = new TGraphErrors*[3];
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lth[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();

  // also get full result for comparison
  // get the fit results
  // get lambda values for each eff model
  TGraphErrors **graph_lthF = new TGraphErrors*[3];
  // get the final lth from the base SB/MC fit
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lthF[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  fIndB->Close();
  TFile *fInd1 = new TFile("../../../Simult_dR1/PR_fit/files/finalFitRes.root");
  graph_lthF[2] = (TGraphErrors*)fInd1->Get(Form("graph_lambda_J"));
  fInd1->Close();
  TFile *fInd2 = new TFile("../../../Simult_dR2/PR_fit/files/finalFitRes.root");
  graph_lthF[1] = (TGraphErrors*)fInd2->Get(Form("graph_lambda_J"));
  fInd2->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle(Form("Run 2 #lambda_{#theta} (coarse)"));

  int col[] = {kBlack, kBlue, kRed};
  for(int i = 0; i < 3; i++) {
    graph_lth[i]->SetLineColor(col[i]);
    graph_lth[i]->SetMarkerColor(col[i]);
    graph_lth[i]->SetMarkerStyle(20);
    graph_lth[i]->SetMarkerSize(.5);
    graph_lth[i]->Draw("p same");
  }

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  TLine *trans1 = new TLine(46, -1, 46, 1);
  trans1->SetLineColor(kBlack);
  trans1->SetLineStyle(kDashed);
  trans1->Draw();
  TLine *trans2 = new TLine(66, -1, 66, 1);
  trans2->SetLineColor(kBlack);
  trans2->SetLineStyle(kDashed);
  trans2->Draw();

  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(graph_lth[0], "baseline", "pl");
  leg->AddEntry(graph_lth[1], "#DeltaR>0.15", "pl");
  leg->AddEntry(graph_lth[2], "#DeltaR>0.17", "pl");
  leg->Draw();
  
  c->SaveAs("plots/ratioFinal/par_lth.pdf");
  //c->Clear();

  for(int i = 0; i < 3; i++) {
    graph_lth[i]->SetLineStyle(kDashDotted);

    graph_lthF[i]->SetLineColor(col[i]);
    graph_lthF[i]->SetMarkerColor(col[i]);
    graph_lthF[i]->SetMarkerStyle(20);
    graph_lthF[i]->SetMarkerSize(.75);
    graph_lthF[i]->Draw("p same");
  }

  c->SaveAs("plots/ratioFinal/par_lth_F.pdf");
  c->Clear();
  
  c->Destructor();
  
  fIn->Close();


}
