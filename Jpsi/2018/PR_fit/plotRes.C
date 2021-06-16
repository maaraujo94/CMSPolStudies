// code to plot the fit results

void plotRes()
{
  // get the histo limits
  TFile *fIn = new TFile("files/ratioHist.root");
  TH2D* rHist;
  fIn->GetObject("ratioHist_ab", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get A, lambda, chiProb values for each bin
  TFile *fInd = new TFile("files/fit_res_1d.root");
  TGraphErrors* graph_A = (TGraphErrors*)fInd->Get("graph_A");
  TGraphErrors* graph_lth = (TGraphErrors*)fInd->Get("graph_lambda");
  TGraphErrors* graph_chi = (TGraphErrors*)fInd->Get("graph_chiP");
  
  double* chisquare = ((TGraphErrors*)fInd->Get("graph_chisquare"))->GetY();
  double* ndf = ((TGraphErrors*)fInd->Get("graph_NDF"))->GetY();

  TH1D* pHist[nBinspT];
  for(int i = 0; i < nBinspT; i++)
    pHist[i] = (TH1D*)fInd->Get(Form("fine_bin%d_1d_min", i+1));

  // draw the fit results
  TCanvas *c = new TCanvas("name", "title", 700, 700);

  // draw the distributions
  for(int i = 0; i < nBinspT; i++) {
    double aVal = graph_A->GetY()[i];
    
    pHist[i]->SetStats(0);
    pHist[i]->SetMinimum(0);
    pHist[i]->SetMaximum(aVal*1.6);
    pHist[i]->SetTitle(Form("2018 PR/MC %s", pHist[i]->GetTitle()));
    pHist[i]->Draw();
    
    TLatex lc;
    lc.SetTextSize(0.03);
    lc.DrawLatex(0.1, aVal*0.3, Form("#chi^{2}/ndf = %.0f/%.0f", chisquare[i], ndf[i]));
    lc.DrawLatex(0.1, aVal*0.2, Form("P(#chi^{2},ndf) = %.1f%%", 100*TMath::Prob(chisquare[i], (int)ndf[i])));

    c->SaveAs(Form("plots/ratio_final/bin_%d.pdf", i+1));
  }
 
  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("2018 #lambda_{#theta} (PR)");

  graph_lth->SetLineColor(kBlack);
  graph_lth->SetMarkerColor(kBlack);
  graph_lth->Draw("p");

  TF1 *cons = new TF1("constant", "[0]", pTBins[0], pTBins[nBinspT]);
  cons->SetLineColor(kBlue);
  graph_lth->Fit(cons);
  
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

  TLatex lc;
  lc.SetTextSize(0.03);
  lc.DrawLatex(70, 0.85, Form("#lambda_{#theta}^{PR} = %.3f #pm %.3f", cons->GetParameter(0), cons->GetParError(0)));
  lc.DrawLatex(70, 0.7, Form("#chi^{2}/ndf = %.0f/%d", cons->GetChisquare(), cons->GetNDF()));
  lc.DrawLatex(70, 0.55, Form("P(#chi^{2},ndf) = %.1f%%", 100*TMath::Prob(cons->GetChisquare(), cons->GetNDF())));
  
  c->SaveAs("plots/ratio_final/lth.pdf");
  c->Clear();
  
  // draw A(pT)
  c->SetLogy();
  TH1F *fa = c->DrawFrame(pTBins[0], 0.04, pTBins[nBinspT], 6);
  fa->SetXTitle("p_{T} (GeV)");
  fa->SetYTitle("A");
  fa->GetYaxis()->SetTitleOffset(1.3);
  fa->GetYaxis()->SetLabelOffset(0.01);
  fa->SetTitle("A (PR)");

  // combine both lambda_th distributions
  graph_A->SetLineColor(kBlack);
  graph_A->SetMarkerColor(kBlack);
  graph_A->Draw("p");

  TLine *trans1_A = new TLine(46, 0.04, 46, 6);
  trans1_A->SetLineColor(kBlack);
  trans1_A->SetLineStyle(kDashed);
  trans1_A->Draw();
  TLine *trans2_A = new TLine(66, 0.04, 66, 6);
  trans2_A->SetLineColor(kBlack);
  trans2_A->SetLineStyle(kDashed);
  trans2_A->Draw();

  c->SaveAs("plots/ratio_final/A.pdf");
  c->Clear();

  // draw chiProb(pT)
  c->SetLogy();
  TH1F *fc = c->DrawFrame(pTBins[0], 0.02, pTBins[nBinspT], 1);
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

  TLine *trans1_C = new TLine(46, 0.02, 46, 1);
  trans1_C->SetLineColor(kBlack);
  trans1_C->SetLineStyle(kDashed);
  trans1_C->Draw();
  TLine *trans2_C = new TLine(66, 0.02, 66, 1);
  trans2_C->SetLineColor(kBlack);
  trans2_C->SetLineStyle(kDashed);
  trans2_C->Draw();

  c->SaveAs("plots/ratio_final/chiP.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();
  fInd->Close();


}

