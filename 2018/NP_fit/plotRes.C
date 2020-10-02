// code to plot the fit results

void plotRes()
{
  // get the histo limits
  TFile *fIn = new TFile("files/ratioHist.root");
  TH2D* rHist;
  TH2D* rHist_hpt;
  fIn->GetObject("ratioHist_ab", rHist);
  fIn->GetObject("ratioHist_ab_hpt", rHist_hpt);

  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  int nBinspT_hpt = rHist_hpt->GetNbinsY();
  const double *pTBins_hpt = rHist_hpt->GetYaxis()->GetXbins()->GetArray();

  
  // get the fit results
  // get A, lambda, chiProb values for each bin
  TFile *fInd = new TFile("files/fit_res_1d.root");
  TGraphErrors* graph_A = (TGraphErrors*)fInd->Get("graph_A");
  TGraphErrors* graph_lth = (TGraphErrors*)fInd->Get("graph_lambda");
  TGraphErrors* graph_chi = (TGraphErrors*)fInd->Get("graph_chiP");
  TGraphErrors* graph_A_hpt = (TGraphErrors*)fInd->Get("graph_A_hpt");
  TGraphErrors* graph_lth_hpt = (TGraphErrors*)fInd->Get("graph_lambda_hpt");
  TGraphErrors* graph_chi_hpt = (TGraphErrors*)fInd->Get("graph_chiP_hpt");

  double* chisquare = ((TGraphErrors*)fInd->Get("graph_chisquare"))->GetY();
  double* ndf = ((TGraphErrors*)fInd->Get("graph_NDF"))->GetY();
  double* chisquare_hpt = ((TGraphErrors*)fInd->Get("graph_chisquare_hpt"))->GetY();
  double* ndf_hpt = ((TGraphErrors*)fInd->Get("graph_NDF_hpt"))->GetY();
  
  TH1D* pHist[nBinspT+nBinspT_hpt];
  for(int i = 0; i < nBinspT+nBinspT_hpt; i++)
    pHist[i] = (TH1D*)fInd->Get(Form("fine_bin%d_1d_min", i+1));

  // draw the fit results
  TCanvas *c = new TCanvas("name", "title", 700, 700);

  // draw the distributions
  for(int i = 0; i < nBinspT+nBinspT_hpt; i++) {
    double chis = i < nBinspT ? chisquare[i] : chisquare_hpt[i-nBinspT];
    int nf = i < nBinspT ? (int)ndf[i] : (int)ndf_hpt[i-nBinspT];
    double aVal = i < nBinspT ? graph_A->GetY()[i] : graph_A_hpt->GetY()[i-nBinspT];
    
    pHist[i]->SetStats(0);
    pHist[i]->SetMaximum(aVal*1.4);
    pHist[i]->Draw();
    
    TLatex lc;
    lc.SetTextSize(0.03);
    lc.DrawLatex(0.1, aVal*0.3, Form("#chi^{2}/ndf = %.0f/%d", chis, nf));
    lc.DrawLatex(0.1, aVal*0.2, Form("P(#chi^{2},ndf) = %.1f%%", 100*TMath::Prob(chis, nf)));

    c->SaveAs(Form("plots/ratio_final/bin_%d.pdf", i+1));
  }
 
  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0], -1, pTBins_hpt[nBinspT_hpt], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("#lambda_{#theta} (NP)");

  // combine both lambda_th distributions
  graph_lth->SetLineColor(kBlack);
  graph_lth->SetMarkerColor(kBlack);
  graph_lth->Draw("p");
  graph_lth_hpt->SetLineColor(kBlack);
  graph_lth_hpt->SetMarkerColor(kBlack);
  graph_lth_hpt->Draw("p");

  TLine *zero = new TLine(pTBins[0], 0, pTBins_hpt[nBinspT_hpt], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  TLine *trans = new TLine(pTBins_hpt[0], -1, pTBins_hpt[0], 1);
  trans->SetLineColor(kBlack);
  trans->SetLineStyle(kDashed);
  trans->Draw();

  c->SaveAs("plots/ratio_final/lth.pdf");
  c->Clear();

  // draw A(pT)
  c->SetLogy();
  TH1F *fa = c->DrawFrame(pTBins[0], 0.3, pTBins_hpt[nBinspT_hpt], 9);
  fa->SetXTitle("p_{T} (GeV)");
  fa->SetYTitle("A");
  fa->GetYaxis()->SetTitleOffset(1.3);
  fa->GetYaxis()->SetLabelOffset(0.01);
  fa->SetTitle("A (NP)");

  // combine both lambda_th distributions
  graph_A->SetLineColor(kBlack);
  graph_A->SetMarkerColor(kBlack);
  graph_A->Draw("p");
  graph_A_hpt->SetLineColor(kBlack);
  graph_A_hpt->SetMarkerColor(kBlack);
  graph_A_hpt->Draw("p");

  TLine *trans_A = new TLine(pTBins_hpt[0], 0.3, pTBins_hpt[0], 9);
  trans_A->SetLineColor(kBlack);
  trans_A->SetLineStyle(kDashed);
  trans_A->Draw();

  c->SaveAs("plots/ratio_final/A.pdf");
  c->Clear();

  // draw chiProb(pT)
  c->SetLogy(0);
  TH1F *fc = c->DrawFrame(pTBins[0], 0, pTBins_hpt[nBinspT_hpt], 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("P(#chi^{2}, ndf)");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);
  fc->SetTitle("P(#chi^{2}, ndf) (NP)");

  // combine both lambda_th distributions
  graph_chi->SetLineColor(kBlack);
  graph_chi->SetMarkerColor(kBlack);
  graph_chi->SetMarkerStyle(20);
  graph_chi->SetMarkerSize(0.75);
  graph_chi->Draw("p");
  graph_chi_hpt->SetLineColor(kBlack);
  graph_chi_hpt->SetMarkerColor(kBlack);
  graph_chi_hpt->SetMarkerStyle(20);
  graph_chi_hpt->SetMarkerSize(0.75);
  graph_chi_hpt->Draw("p");

  TLine *trans_C = new TLine(pTBins_hpt[0], 0, pTBins_hpt[0], 1);
  trans_C->SetLineColor(kBlack);
  trans_C->SetLineStyle(kDashed);
  trans_C->Draw();

  c->SaveAs("plots/ratio_final/chiP.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();
  fInd->Close();


}

