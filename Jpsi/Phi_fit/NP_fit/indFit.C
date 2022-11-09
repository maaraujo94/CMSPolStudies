// code to do the individual fit (1d phi maps)

// main
void indFit()
{
  // read the histos from subtraction
  TFile *infile = new TFile("files/bkgSubRes.root");
  TH2D **h_fit = new TH2D*[3];
  string lbl[] = {"NP", "NPc", "SB"};
  for(int i = 0; i < 3; i++) {
    infile->GetObject(Form("h_%s", lbl[i].c_str()), h_fit[i]);
    h_fit[i]->SetDirectory(0);
  }
  infile->Close();

  // get the binning
  int nBinsX = h_fit[0]->GetNbinsX(), nBinsY = h_fit[0]->GetNbinsY();
  const double *yBins = h_fit[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *pHist[3][nBinsY];
  for(int i_t = 0; i_t < 3; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      pHist[i_t][i-1] = h_fit[i_t]->ProjectionX(Form("bin%d_%d", i, i_t+1), i, i);
      pHist[i_t][i-1]->SetTitle(Form("%s bin %d: [%.1f, %.1f] GeV", lbl[i_t].c_str(), i, yBins[i-1], yBins[i]));
    }
  }
  
  // the fit function to be used
  TF1 **fit1d = new TF1*[2];
  for(int i = 0; i < 2; i++) {
    fit1d[i] = new TF1(Form("fit_%d", i), "[0]*(1+[1]*cos(2.*x*3.1415/180.))", -180, 180);
    fit1d[i]->SetParNames("A", "l_phi");
  }
   
  // the cycle to fit each bin and store fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  TFile *outfile = new TFile("files/finalFitRes.root", "recreate");

  double parA[2][nBinsY], eparA[2][nBinsY];
  double parL[2][nBinsY], eparL[2][nBinsY];
  double chi2[2][nBinsY], ndf[2][nBinsY], chiP[2][nBinsY];
  double pt[nBinsY], ept[nBinsY];
  
  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_fit[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_fit[0]->GetYaxis()->GetBinUpEdge(i+1);
    pt[i] = (pMax+pMin)/2.;
    ept[i] = (pMax-pMin)/2.;

    // fit the 2 functions
    for(int i_t = 0; i_t < 2; i_t++) {
      fit1d[i_t]->SetParameters(pHist[i_t][i]->GetBinContent(1)*1.1, 0.1);

      pHist[i_t][i]->Fit(fit1d[i_t], "R0");

      parA[i_t][i] = fit1d[i_t]->GetParameter(0);
      eparA[i_t][i] = fit1d[i_t]->GetParError(0);
      parL[i_t][i] = fit1d[i_t]->GetParameter(1);
      eparL[i_t][i] = fit1d[i_t]->GetParError(1);
      chi2[i_t][i] = fit1d[i_t]->GetChisquare();
      ndf[i_t][i] = fit1d[i_t]->GetNDF();
      chiP[i_t][i] = TMath::Prob(chi2[i_t][i], ndf[i_t][i]);
    }

    // plotting everything
    pHist[0][i]->SetTitle(Form("data/MC #phi (%.1f < p_{T} < %.1f GeV)", pMin, pMax));
    pHist[0][i]->SetStats(0);
    pHist[0][i]->SetLineColor(kRed);
    pHist[0][i]->SetMarkerColor(kRed);
    pHist[0][i]->SetMinimum(0);
    pHist[0][i]->SetMaximum(pHist[0][i]->GetBinContent(1)*1.5);
    pHist[0][i]->GetXaxis()->SetTitle("#phi_{HX}");
    pHist[0][i]->Draw("error");
    fit1d[0]->SetLineColor(kRed);
    fit1d[0]->SetLineStyle(kDashed);
    fit1d[0]->Draw("same");

    pHist[1][i]->SetLineColor(kRed+3);
    pHist[1][i]->SetMarkerColor(kRed+3);
    pHist[1][i]->Draw("same");
    fit1d[1]->SetLineColor(kRed+3);
    fit1d[1]->SetLineStyle(kDashed);
    fit1d[1]->Draw("same");

    pHist[2][i]->SetLineColor(kGreen);
    pHist[2][i]->SetMarkerColor(kGreen);
    pHist[2][i]->Draw("same");
    
    TLatex lc;
    lc.SetTextSize(0.03);
    lc.DrawLatex(-150, pHist[0][i]->GetMaximum()*0.9, Form("#lambda_{#phi}^{NP} = %.3f #pm %.3f", parL[0][i], eparL[0][i]));
    lc.DrawLatex(-150, pHist[0][i]->GetMaximum()*0.8, Form("#lambda_{#phi}^{pure NP} = %.3f #pm %.3f", parL[1][i], eparL[1][i]));
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.97, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(pHist[0][i], "NP", "pl");
    leg->AddEntry(pHist[2][i], "bkg", "pl");
    leg->AddEntry(pHist[1][i], "non-prompt J/#psi", "pl");
    leg->Draw();

    for(int i_t = 0; i_t < 3; i_t++)
      pHist[i_t][i]->Write();
 
    c->SaveAs(Form("plots/ratioFinal/bin_%d.pdf", i));
    c->Clear();
    cout << endl << endl;
  }

  for(int i_t = 0; i_t < 2; i_t++) {
    // make and save the TGraph with the fit results and max costh used
    TGraphErrors *graphA = new TGraphErrors(nBinsY, pt, parA[i_t], ept, eparA[i_t]);
    TGraphErrors *graphL = new TGraphErrors(nBinsY, pt, parL[i_t], ept, eparL[i_t]);
    TGraph *graphC = new TGraph(nBinsY, pt, chi2[i_t]);
    TGraph *graphN = new TGraph(nBinsY, pt, ndf[i_t]);
    TGraph *graphP = new TGraph(nBinsY, pt, chiP[i_t]);

    graphA->SetName(Form("graph_A_%s", lbl[i_t].c_str()));
    graphL->SetName(Form("graph_lambda_%s", lbl[i_t].c_str()));
    graphC->SetName(Form("graph_chisquare_%s", lbl[i_t].c_str()));
    graphN->SetName(Form("graph_NDF_%s", lbl[i_t].c_str()));
    graphP->SetName(Form("graph_chiP_%s", lbl[i_t].c_str()));
 
    graphA->Write();
    graphL->Write();
    graphC->Write();
    graphN->Write();
    graphP->Write();
  }
  outfile->Close();

  c->Destructor();
}
