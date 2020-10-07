// code to plot the fit results

void plotRes()
{
  // get the histo limits
  TFile *fIn = new TFile("files/ratioHist.root");
  TH2D* rHist;
  TH2D* rHist_hpt;
  TH2D* rHist_vhpt;
  fIn->GetObject("ratioHist_ab", rHist);
  fIn->GetObject("ratioHist_ab_hpt", rHist_hpt);
  fIn->GetObject("ratioHist_ab_vhpt", rHist_vhpt);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  int nBinspT_hpt = rHist_hpt->GetNbinsY();
  const double *pTBins_hpt = rHist_hpt->GetYaxis()->GetXbins()->GetArray();
  int nBinspT_vhpt = rHist_vhpt->GetNbinsY();
  const double *pTBins_vhpt = rHist_vhpt->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get A, lambda, chiProb values for each bin
  TFile *fInd = new TFile("files/fit_res_1d.root");
  TGraphErrors* graph_A = (TGraphErrors*)fInd->Get("graph_A");
  TGraphErrors* graph_lth = (TGraphErrors*)fInd->Get("graph_lambda");
  TGraphErrors* graph_chi = (TGraphErrors*)fInd->Get("graph_chiP");
  TGraphErrors* graph_A_hpt = (TGraphErrors*)fInd->Get("graph_A_hpt");
  TGraphErrors* graph_lth_hpt = (TGraphErrors*)fInd->Get("graph_lambda_hpt");
  TGraphErrors* graph_chi_hpt = (TGraphErrors*)fInd->Get("graph_chiP_hpt");
  TGraphErrors* graph_A_vhpt = (TGraphErrors*)fInd->Get("graph_A_vhpt");
  TGraphErrors* graph_lth_vhpt = (TGraphErrors*)fInd->Get("graph_lambda_vhpt");
  TGraphErrors* graph_chi_vhpt = (TGraphErrors*)fInd->Get("graph_chiP_vhpt");

  double* chisquare = ((TGraphErrors*)fInd->Get("graph_chisquare"))->GetY();
  double* ndf = ((TGraphErrors*)fInd->Get("graph_NDF"))->GetY();
  double* chisquare_hpt = ((TGraphErrors*)fInd->Get("graph_chisquare_hpt"))->GetY();
  double* ndf_hpt = ((TGraphErrors*)fInd->Get("graph_NDF_hpt"))->GetY();
  double* chisquare_vhpt = ((TGraphErrors*)fInd->Get("graph_chisquare_vhpt"))->GetY();
  double* ndf_vhpt = ((TGraphErrors*)fInd->Get("graph_NDF_vhpt"))->GetY();

  int ntot = nBinspT+nBinspT_hpt+nBinspT_vhpt;
  TH1D* pHist[ntot];
  for(int i = 0; i < ntot; i++)
    pHist[i] = (TH1D*)fInd->Get(Form("fine_bin%d_1d_min", i+1));

  // draw the fit results
  TCanvas *c = new TCanvas("name", "title", 700, 700);

  // draw the distributions
  for(int i = 0; i < ntot; i++) {
    double chis = i < nBinspT ? chisquare[i] : i < nBinspT+nBinspT_hpt ? chisquare_hpt[i-nBinspT] : chisquare_vhpt[i-nBinspT-nBinspT_hpt];
    int nf = i < nBinspT ? (int)ndf[i] : i < nBinspT+nBinspT_hpt ? (int)ndf_hpt[i-nBinspT] : (int)ndf_vhpt[i-nBinspT-nBinspT_hpt];
    double aVal = i < nBinspT ? graph_A->GetY()[i] : i < nBinspT+nBinspT_hpt ? graph_A_hpt->GetY()[i-nBinspT] : graph_A_vhpt->GetY()[i-nBinspT-nBinspT_hpt];

    cout << i << " " << aVal << endl;
    
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
  TH1F *fl = c->DrawFrame(pTBins[0], -1, pTBins_vhpt[nBinspT_vhpt], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("#lambda_{#theta} (PR)");

  // combine both lambda_th distributions
  graph_lth->SetLineColor(kBlack);
  graph_lth->SetMarkerColor(kBlack);
  graph_lth->Draw("p");
  graph_lth_hpt->SetLineColor(kBlack);
  graph_lth_hpt->SetMarkerColor(kBlack);
  graph_lth_hpt->Draw("p");
  graph_lth_vhpt->SetLineColor(kBlack);
  graph_lth_vhpt->SetMarkerColor(kBlack);
  graph_lth_vhpt->Draw("p");

  TLine *zero = new TLine(pTBins[0], 0, pTBins_vhpt[nBinspT_vhpt], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  TLine *trans1 = new TLine(pTBins_hpt[0], -1, pTBins_hpt[0], 1);
  trans1->SetLineColor(kBlack);
  trans1->SetLineStyle(kDashed);
  trans1->Draw();
  TLine *trans2 = new TLine(pTBins_vhpt[0], -1, pTBins_vhpt[0], 1);
  trans2->SetLineColor(kBlack);
  trans2->SetLineStyle(kDashed);
  trans2->Draw();

  c->SaveAs("plots/ratio_final/lth.pdf");
  c->Clear();

  // draw A(pT)
  c->SetLogy();
  TH1F *fa = c->DrawFrame(pTBins[0], 0.04, pTBins_vhpt[nBinspT_vhpt], 6);
  fa->SetXTitle("p_{T} (GeV)");
  fa->SetYTitle("A");
  fa->GetYaxis()->SetTitleOffset(1.3);
  fa->GetYaxis()->SetLabelOffset(0.01);
  fa->SetTitle("A (PR)");

  // combine both lambda_th distributions
  graph_A->SetLineColor(kBlack);
  graph_A->SetMarkerColor(kBlack);
  graph_A->Draw("p");
  graph_A_hpt->SetLineColor(kBlack);
  graph_A_hpt->SetMarkerColor(kBlack);
  graph_A_hpt->Draw("p");
  graph_A_vhpt->SetLineColor(kBlack);
  graph_A_vhpt->SetMarkerColor(kBlack);
  graph_A_vhpt->Draw("p");

  TLine *trans1_A = new TLine(pTBins_hpt[0], 0.04, pTBins_hpt[0], 6);
  trans1_A->SetLineColor(kBlack);
  trans1_A->SetLineStyle(kDashed);
  trans1_A->Draw();
  TLine *trans2_A = new TLine(pTBins_vhpt[0], 0.04, pTBins_vhpt[0], 6);
  trans2_A->SetLineColor(kBlack);
  trans2_A->SetLineStyle(kDashed);
  trans2_A->Draw();

  c->SaveAs("plots/ratio_final/A.pdf");
  c->Clear();

  // draw chiProb(pT)
  c->SetLogy(0);
  TH1F *fc = c->DrawFrame(pTBins[0], 0, pTBins_vhpt[nBinspT_vhpt], 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("P(#chi^{2}, ndf)");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);
  fc->SetTitle("P(#chi^{2}, ndf) (PR)");

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
  graph_chi_vhpt->SetLineColor(kBlack);
  graph_chi_vhpt->SetMarkerColor(kBlack);
  graph_chi_vhpt->SetMarkerStyle(20);
  graph_chi_vhpt->SetMarkerSize(0.75);
  graph_chi_vhpt->Draw("p");

  TLine *trans1_C = new TLine(pTBins_hpt[0], 0, pTBins_hpt[0], 1);
  trans1_C->SetLineColor(kBlack);
  trans1_C->SetLineStyle(kDashed);
  trans1_C->Draw();
  TLine *trans2_C = new TLine(pTBins_vhpt[0], 0, pTBins_vhpt[0], 1);
  trans2_C->SetLineColor(kBlack);
  trans2_C->SetLineStyle(kDashed);
  trans2_C->Draw();

  c->SaveAs("plots/ratio_final/chiP.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();
  fInd->Close();


}

