#import "../../cosMax/imp_jumpF.C"

// code to do the individual fit (1d costheta maps)

// main
void indFit()
{
  // read the histos from subtraction
  TFile *infile = new TFile("files/bkgSubRes.root");
  TH2D **h_fit = new TH2D*[5];
  string lbl[] = {"Data", "NP", "PR", "J", "SB"};
  for(int i = 0; i < 5; i++) {
    infile->GetObject(Form("h_%s", lbl[i].c_str()), h_fit[i]);
    h_fit[i]->SetDirectory(0);
  }
  infile->Close();

  // get the binning
  int nBinsX = h_fit[0]->GetNbinsX(), nBinsY = h_fit[0]->GetNbinsY();
  const double *yBins = h_fit[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *pHist[5][nBinsY];
  for(int i_t = 0; i_t < 5; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      pHist[i_t][i-1] = h_fit[i_t]->ProjectionX(Form("bin%d_%d", i, i_t+1), i, i);
      pHist[i_t][i-1]->SetTitle(Form("%s bin %d: [%.1f, %.1f] GeV", lbl[i_t].c_str(), i, yBins[i-1], yBins[i]));
    }
  }
  
  // the fit function to be used - only on total and prompt J/psi
  TF1 **fit1d = new TF1*[3];
  for(int i = 0; i < 3; i++) {
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

  double parA[3][nBinsY], eparA[3][nBinsY];
  double parL[3][nBinsY], eparL[3][nBinsY];
  double chi2[3][nBinsY], ndf[3][nBinsY], chiP[3][nBinsY];
  double pt[nBinsY], ept[nBinsY], cMax[nBinsY];
  
  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_fit[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_fit[0]->GetYaxis()->GetBinUpEdge(i+1);
    pt[i] = (pMax+pMin)/2.;
    ept[i] = (pMax-pMin)/2.;

    // get max costheta
    double cMaxVal = jumpF(cosMax->Integral(pMin, pMax)/(pMax-pMin));
    cMax[i] = cMaxVal;
    
    // fit the 3 functions
    for(int i_t = 0; i_t < 2; i_t++) {
      int i_fit = 3*i_t; // fit 0 and 3

      fit1d[i_t]->SetRange(0, cMaxVal);
      fit1d[i_t]->SetParameters(pHist[i_fit][i]->GetBinContent(1)*1.1, 0.1);

      pHist[i_fit][i]->Fit(fit1d[i_t], "R0");

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
    pHist[0][i]->SetLineColor(kViolet);
    pHist[0][i]->SetMarkerColor(kViolet);
    pHist[0][i]->SetMinimum(0);
    pHist[0][i]->SetMaximum(pHist[0][i]->GetBinContent(cMaxVal*nBinsX)*1.8);
    pHist[0][i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    pHist[0][i]->Draw("error");
    fit1d[0]->SetLineColor(kViolet);
    fit1d[0]->SetLineStyle(kDashed);
    fit1d[0]->Draw("same");

    pHist[1][i]->SetLineColor(kRed);
    pHist[1][i]->SetMarkerColor(kRed);
    pHist[1][i]->Draw("same");

    pHist[2][i]->SetLineColor(kBlack);
    pHist[2][i]->SetMarkerColor(kBlack);
    pHist[2][i]->Draw("same");

    pHist[4][i]->SetLineColor(kGreen);
    pHist[4][i]->SetMarkerColor(kGreen);
    pHist[4][i]->Draw("same");

    pHist[3][i]->SetLineColor(kBlue);
    pHist[3][i]->SetMarkerColor(kBlue);
    pHist[3][i]->Draw("same");
    fit1d[1]->SetLineColor(kBlue);
    fit1d[1]->SetLineStyle(kDashed);
    fit1d[1]->Draw("same");

    TLatex lc;
    lc.SetTextSize(0.03);
    lc.DrawLatex(0.1, pHist[0][i]->GetMaximum()*0.9, Form("#lambda_{#theta}^{total} = %.3f #pm %.3f", parL[0][i], eparL[0][i]));
    lc.DrawLatex(0.1, pHist[0][i]->GetMaximum()*0.8, Form("#lambda_{#theta}^{prompt #psi(2S)} = %.3f #pm %.3f", parL[1][i], eparL[1][i]));
    
    TLine *c_lim = new TLine(cMaxVal, 0, cMaxVal, pHist[0][i]->GetMaximum());
    c_lim->SetLineStyle(kDashed);
    c_lim->SetLineColor(kBlack);
    c_lim->Draw();

    TLegend *leg = new TLegend(0.74, 0.7, 1.04, 0.9);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(kWhite,0);
    leg->AddEntry(pHist[0][i], "total", "pl");
    leg->AddEntry(pHist[1][i], "NP contrib", "pl");
    leg->AddEntry(pHist[2][i], "prompt", "pl");
    leg->AddEntry(pHist[4][i], "SB contrib", "pl");
    leg->AddEntry(pHist[3][i], "prompt #psi(2S)", "pl");
    leg->Draw();

    for(int i_t = 0; i_t < 5; i_t++)
      pHist[i_t][i]->Write();
 
    c->SaveAs(Form("plots/ratioFinal/bin_%d.pdf", i));
    c->Clear();
    cout << endl << endl;
  }
  outfile->Close();

  // NP needs to come from the corresponding input
  TFile *infile_NP = new TFile("files/bkgSubRes_NP.root");
  infile_NP->GetObject(Form("h_NPc"), h_fit[1]);
  h_fit[1]->SetDirectory(0);
  infile_NP->Close();

  // get the 1d plots
  for(int i_pt = 1; i_pt <= nBinsY; i_pt++) {
    pHist[1][i_pt-1] = h_fit[1]->ProjectionX(Form("bin%d_%d", i_pt, 2), i_pt, i_pt);
    pHist[1][i_pt-1]->SetTitle(Form("NP bin %d: [%.0f, %.0f] GeV", i_pt, yBins[i_pt-1], yBins[i_pt]));
  }

  // new cycle, now for NP fit
  for(int i = 0; i < nBinsY; i++) {
    double pMin = h_fit[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_fit[0]->GetYaxis()->GetBinUpEdge(i+1);

    // need to fit NP first
    fit1d[2]->SetRange(0, cMax[i]);
    fit1d[2]->SetParameters(pHist[1][i]->GetBinContent(1)*1.1, 0.1);
    pHist[1][i]->Fit(fit1d[2], "R0");
    parA[2][i] = fit1d[2]->GetParameter(0);
    eparA[2][i] = fit1d[2]->GetParError(0);
    parL[2][i] = fit1d[2]->GetParameter(1);
    eparL[2][i] = fit1d[2]->GetParError(1);
    chi2[2][i] = fit1d[2]->GetChisquare();
    ndf[2][i] = fit1d[2]->GetNDF();
    chiP[2][i] = TMath::Prob(chi2[2][i], ndf[2][i]);

  }
  
  string lbl_s[] = {"Data", "J", "NP"};
  TFile *outfile2 = new TFile("files/finalFitRes.root", "update");
  for(int i_t = 0; i_t < 3; i_t++) {
    // make and save the TGraph with the fit results and max costh used
    TGraphErrors *graphA = new TGraphErrors(nBinsY, pt, parA[i_t], ept, eparA[i_t]);
    TGraphErrors *graphL = new TGraphErrors(nBinsY, pt, parL[i_t], ept, eparL[i_t]);
    TGraph *graphC = new TGraph(nBinsY, pt, chi2[i_t]);
    TGraph *graphN = new TGraph(nBinsY, pt, ndf[i_t]);
    TGraph *graphP = new TGraph(nBinsY, pt, chiP[i_t]);

    graphA->SetName(Form("graph_A_%s", lbl_s[i_t].c_str()));
    graphL->SetName(Form("graph_lambda_%s", lbl_s[i_t].c_str()));
    graphC->SetName(Form("graph_chisquare_%s", lbl_s[i_t].c_str()));
    graphN->SetName(Form("graph_NDF_%s", lbl_s[i_t].c_str()));
    graphP->SetName(Form("graph_chiP_%s", lbl_s[i_t].c_str()));
 
    graphA->Write();
    graphL->Write();
    graphC->Write();
    graphN->Write();
    graphP->Write();
  }
  outfile2->Close();

  c->Destructor();
}
