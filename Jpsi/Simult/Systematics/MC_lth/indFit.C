#import "../../cosMax/imp_jumpF.C"

// code to do the individual fit (1d costheta maps)

// main
void indFit()
{
  // read the histos from subtraction
  TFile *infile = new TFile("files/histoStore.root");
  TH2D **h_fit = new TH2D*[21];
  for(int i = 0; i < 21; i++) {
    infile->GetObject(Form("rMCH_%d", i), h_fit[i]);
    h_fit[i]->SetDirectory(0);
  }
  infile->Close();

  // get the binning
  int nBinsX = h_fit[0]->GetNbinsX(), nBinsY = h_fit[0]->GetNbinsY();
  const double *yBins = h_fit[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *pHist[21][nBinsY];
  for(int i_t = 0; i_t < 21; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      pHist[i_t][i-1] = h_fit[i_t]->ProjectionX(Form("bin%d_%d", i, i_t+1), i, i);
      pHist[i_t][i-1]->SetTitle(Form("%d bin %d: [%.1f, %.1f] GeV", i_t, i, yBins[i-1], yBins[i]));
    }
  }
  
  // the fit function to be used
  TF1 **fit1d = new TF1*[21];
  for(int i = 0; i < 21; i++) {
    fit1d[i] = new TF1(Form("fit_%d", i), "[0]*(1+[1]*x*x)", 0, 1);
    fit1d[i]->SetParNames("A", "l_th");
  }
  
  // get the fit range from our cosmax(pT)
  ifstream in;
  string dataS;
  in.open("../../cosMax/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", yBins[0]-10, yBins[nBinsY]+10);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);
 
  // the cycle to fit each bin and store fit results
  TCanvas *c = new TCanvas("", "", 700, 700);    
  TFile *outfile = new TFile("files/finalFitRes.root", "recreate");

  double parA[21][nBinsY], eparA[21][nBinsY];
  double parL[21][nBinsY], eparL[21][nBinsY];
  double chi2[21][nBinsY], ndf[21][nBinsY], chiP[21][nBinsY];
  double pt[nBinsY], ept[nBinsY];
  
  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_fit[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_fit[0]->GetYaxis()->GetBinUpEdge(i+1);
    pt[i] = (pMax+pMin)/2.;
    ept[i] = (pMax-pMin)/2.;

    // get max costheta
    double cMaxVal = jumpF(cosMax->Integral(pMin, pMax)/(pMax-pMin));

    // fit the 4 functions
    for(int i_t = 0; i_t < 21; i_t++) {
      fit1d[i_t]->SetRange(0, cMaxVal);
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
    pHist[0][i]->SetTitle(Form("data/MC |cos#theta| (%.0f < p_{T} < %.0f GeV)", pMin, pMax));
    pHist[0][i]->SetStats(0);
    pHist[0][i]->SetLineColor(kBlue);
    pHist[0][i]->SetMarkerColor(kBlue);
    pHist[0][i]->SetMinimum(0);
    pHist[0][i]->SetMaximum(pHist[0][i]->GetBinContent(1)*1.5);
    pHist[0][i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    pHist[0][i]->Draw("error");
    fit1d[0]->SetLineColor(kBlue);
    fit1d[0]->SetLineStyle(kDashed);
    fit1d[0]->Draw("same");

    for(int j = 1; j < 21; j++) {
      pHist[j][i]->SetLineColor(kBlue);
      pHist[j][i]->SetMarkerColor(kBlue);
      pHist[j][i]->Draw("same");
      fit1d[j]->SetLineColor(kBlue);
      fit1d[j]->SetLineStyle(kDashed);
      fit1d[j]->Draw("same");
    }

    TLatex lc;
    lc.SetTextSize(0.03);
    //lc.DrawLatex(0.1, pHist[0][i]->GetMaximum()*0.9, Form("#lambda_{#theta}^{+0.4} = %.3f #pm %.3f", parL[0][i], eparL[0][i]));
    //lc.DrawLatex(0.1, pHist[0][i]->GetMaximum()*0.8, Form("#lambda_{#theta}^{-0.1} = %.3f #pm %.3f", parL[1][i], eparL[1][i]));
    
    TLine *c_lim = new TLine(cMaxVal, 0, cMaxVal, pHist[0][i]->GetMaximum());
    c_lim->SetLineStyle(kDashed);
    c_lim->SetLineColor(kBlack);
    c_lim->Draw();

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(pHist[0][i], "MC (#lambda = +0.4)", "pl");
    leg->AddEntry(pHist[1][i], "MC (#lambda = -0.1)", "pl");
    //leg->Draw();

    for(int i_t = 0; i_t < 21; i_t++)
      pHist[i_t][i]->Write();
 
    c->SaveAs(Form("plots/ratioFinal/bin_%d.pdf", i));
    c->Clear();
    cout << endl << endl;
  }

  for(int i_t = 0; i_t < 21; i_t++) {
    // make and save the TGraph with the fit results and max costh used
    TGraphErrors *graphA = new TGraphErrors(nBinsY, pt, parA[i_t], ept, eparA[i_t]);
    TGraphErrors *graphL = new TGraphErrors(nBinsY, pt, parL[i_t], ept, eparL[i_t]);
    TGraph *graphC = new TGraph(nBinsY, pt, chi2[i_t]);
    TGraph *graphN = new TGraph(nBinsY, pt, ndf[i_t]);
    TGraph *graphP = new TGraph(nBinsY, pt, chiP[i_t]);

    graphA->SetName(Form("graph_A_%d", i_t));
    graphL->SetName(Form("graph_lambda_%d", i_t));
    graphC->SetName(Form("graph_chisquare_%d", i_t));
    graphN->SetName(Form("graph_NDF_%d", i_t));
    graphP->SetName(Form("graph_chiP_%d", i_t));
 
    graphA->Write();
    graphL->Write();
    graphC->Write();
    graphN->Write();
    graphP->Write();
  }
  outfile->Close();

  c->Destructor();
}
