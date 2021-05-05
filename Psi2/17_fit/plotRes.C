// code to plot the fit results

void plotRes()
{
  // get the histo limits
  TFile *fIn = new TFile("files/ratioHist.root");
  TH2D* rHist;
  fIn->GetObject(Form("ratioHist_ab_S"), rHist);
  TH2D *pS = (TH2D*)fIn->Get("PSig_ab");
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get A, lambda, chiProb values for each bin
  TFile *fInd = new TFile("files/fit_res_1d.root");
  TGraphErrors* graph_A     = (TGraphErrors*)fInd->Get("graph_A");
  TGraphErrors* graph_lth   = (TGraphErrors*)fInd->Get("graph_lambda");
  TGraphErrors* graph_chiP  = (TGraphErrors*)fInd->Get("graph_chiP");
  TGraphErrors* graph_chi   = (TGraphErrors*)fInd->Get("graph_chisquare");
  TGraphErrors* graph_chi_c = (TGraphErrors*)fInd->Get("graph_c_chisquare");

  double* chisquare = ((TGraphErrors*)fInd->Get("graph_chisquare"))->GetY();
  double* ndf = ((TGraphErrors*)fInd->Get("graph_NDF"))->GetY();

  int ntot = nBinspT;
  TH1D* pHist[ntot];
  TH1D* pcHist[ntot];
  for(int i = 0; i < ntot; i++) {
    pHist[i] = (TH1D*)fInd->Get(Form("fine_bin%d_1d_min", i+1));
    pcHist[i] = (TH1D*)fInd->Get(Form("fine_bin%d_1d_c", i+1));
  }
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);

  // draw the distributions
  for(int i = 0; i < ntot; i++) {
    double chis = chisquare[i];
    int nf =  (int)ndf[i];
    double aVal = graph_A->GetY()[i];
    
    pHist[i]->SetStats(0);
    pHist[i]->SetMinimum(0);
    pHist[i]->SetMaximum(aVal*1.6);
    pHist[i]->SetMarkerColor(kBlack);
    pHist[i]->SetLineColor(kBlack);
    pHist[i]->SetTitle(Form("2017 PR/MC %s", pHist[i]->GetTitle()));
    pHist[i]->Draw();

    TF1 *func_c = pcHist[i]->GetFunction("fit f 1d");
    func_c->SetLineColor(kBlue);
    func_c->Draw("lsame");
    double chis_c = func_c->GetChisquare();
    int nf_c = func_c->GetNDF();
    
    TLatex lc;
    lc.SetTextSize(0.03);
    lc.DrawLatex(0.1, aVal*0.3, Form("#chi^{2}/ndf = %.0f/%d", chis, nf));
    lc.DrawLatex(0.1, aVal*0.2, Form("P(#chi^{2},ndf) = %.1f%%", 100*TMath::Prob(chis, nf)));
    lc.SetTextColor(kBlue);
    lc.DrawLatex(0.5, aVal*0.3, Form("#chi^{2}/ndf = %.0f/%d", chis_c, nf_c));
    lc.DrawLatex(0.5, aVal*0.2, Form("P(#chi^{2},ndf) = %.1f%%", 100*TMath::Prob(chis_c, nf_c)));
    
    c->SaveAs(Form("plots/ratio_final/bin_%d.pdf", i+1));
  }
 
  int N = graph_lth->GetN();
  double X[N], Y[N], EX[N], EY[N];
  for(int i = 0; i < N; i++) {
    X[i] = graph_lth->GetX()[i];
    Y[i] = graph_lth->GetY()[i];
    EX[i] = graph_lth->GetEX()[i];
    EY[i] = graph_lth->GetEY()[i];
  }
  
  double minX = X[0]-EX[0];
  double maxX = X[N-1]+EX[N-1];
  
  TH1F *fl = c->DrawFrame(20, -1, maxX, 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("2017 #lambda_{#theta} (PR)");

  TGraphErrors* lth_full = new TGraphErrors(N, X, Y, EX, EY);
  lth_full->SetLineColor(kBlack);
  lth_full->SetMarkerColor(kBlack);
  lth_full->Draw("p");

  TLine *zero = new TLine(20, 0, maxX, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TF1 *cons = new TF1("constant", "[0]", minX, maxX);
  cons->SetLineColor(kBlue);
  lth_full->Fit(cons);
  
  TLatex lc;
  lc.SetTextSize(0.03);
  lc.DrawLatex(70, 0.85, Form("#lambda_{#theta}^{PR} = %.3f #pm %.3f", cons->GetParameter(0), cons->GetParError(0)));
  lc.DrawLatex(70, 0.7, Form("#chi^{2}/ndf = %.0f/%d", cons->GetChisquare(), cons->GetNDF()));
  lc.DrawLatex(70, 0.55, Form("P(#chi^{2},ndf) = %.1f%%", 100*TMath::Prob(cons->GetChisquare(), cons->GetNDF())));
  
  c->SaveAs("plots/ratio_final/lth.pdf");

  // draw A(pT)
  c->SetLogy(0);
  TH1F *fa = c->DrawFrame(20, 0, pTBins[nBinspT], 0.15);
  fa->SetXTitle("p_{T} (GeV)");
  fa->SetYTitle("A");
  fa->GetYaxis()->SetTitleOffset(1.3);
  fa->GetYaxis()->SetLabelOffset(0.01);
  fa->SetTitle("2017 A (PR)");

  // combine both lambda_th distributions
  graph_A->SetLineColor(kBlack);
  graph_A->SetMarkerColor(kBlack);
  graph_A->Draw("p");

  c->SaveAs("plots/ratio_final/A.pdf");
  c->Clear();

  // draw chiProb(pT)
  c->SetLogy();
  TH1F *fc = c->DrawFrame(20, 0.02, pTBins[nBinspT], 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("P(#chi^{2}, ndf)");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);
  fc->SetTitle("2017 P(#chi^{2}, ndf) (PR)");

  // combine both lambda_th distributions
  graph_chiP->SetLineColor(kBlack);
  graph_chiP->SetMarkerColor(kBlack);
  graph_chiP->SetMarkerStyle(20);
  graph_chiP->SetMarkerSize(0.75);
  graph_chiP->Draw("p");
 
  c->SaveAs("plots/ratio_final/chiP.pdf");
  c->Clear();

  // draw chisquare(pT)
  c->SetLogy(0);
  TH1F *fcc = c->DrawFrame(20, 0, pTBins[nBinspT], 40);
  fcc->SetXTitle("p_{T} (GeV)");
  fcc->SetYTitle("#chi^{2}");
  fcc->GetYaxis()->SetTitleOffset(1.3);
  fcc->GetYaxis()->SetLabelOffset(0.01);
  fcc->SetTitle("2017 #chi^{2} (PR)");

  // combine both lambda_th distributions
  graph_chi->SetLineColor(kRed);
  graph_chi->SetMarkerColor(kRed);
  graph_chi->SetMarkerStyle(20);
  graph_chi->SetMarkerSize(0.75);
  graph_chi->Draw("pl");
  graph_chi_c->SetLineColor(kBlue);
  graph_chi_c->SetMarkerColor(kBlue);
  graph_chi_c->SetMarkerStyle(20);
  graph_chi_c->SetMarkerSize(0.75);
  graph_chi_c->Draw("pl");

  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(graph_chi, "free #lambda_{#theta}", "l");
  leg->AddEntry(graph_chi_c, "#lambda_{#theta} = 0", "l");
  leg->Draw();
 
  c->SaveAs("plots/ratio_final/chisquare.pdf");
  c->Clear();

  // plot the 2D pure signal / MC dist
  pS->SetStats(0);
  pS->SetTitle("2017 Data/MC (Pure Signal)");
  pS->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  pS->GetYaxis()->SetTitle("p_{T} (GeV)");
  pS->Draw("COLZ");
  c->SaveAs("plots/PSig_2d_abs.pdf");
  c->Clear();

  c->Destructor();
  
  fIn->Close();
  fInd->Close();
}

