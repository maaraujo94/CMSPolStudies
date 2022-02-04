#import "../cosMax/imp_jumpF.C"
#import "../rcut.C"

// code to do the individual fit (1d costheta maps)

// aux func for costheta_min
double cminf(double pt, double a, double b, double c, double d)
{
  if (pt < d) return 0;
  else
    return a*(1.-exp(b+c*pt));
}

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
      pHist[i_t][i-1]->SetTitle(Form("%s bin %d: [%.0f, %.0f] GeV", lbl[i_t].c_str(), i, yBins[i-1], yBins[i]));
    }
  }
  
  // the fit function to be used
  TF1 **fit1d = new TF1*[2];
  for(int i = 0; i < 2; i++) {
    fit1d[i] = new TF1(Form("fit_%d", i), "[0]*(1+[1]*x*x)", 0, 1);
    fit1d[i]->SetParNames("A", "l_th");
  }
  
  // get the fit range from our cosmax(pT), cosmin(pT)
  ifstream in;
  string dataS;
  in.open("../cosMax/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", yBins[0]-10, yBins[nBinsY]+10);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);

  in.open("../cosMax/cosMinFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double minPar[4];
  in >> minPar[0] >> aux >> minPar[1] >> aux >> minPar[2] >> aux >> minPar[3];
  in.close();
  
  TF1 *cosMin = new TF1("cosMin", "cminf(x, [0], [1], [2], [3])", yBins[0]-10, yBins[nBinsY]+10);
  cosMin->SetParameters(minPar[0], minPar[1], minPar[2], minPar[3]);

 
  // the cycle to fit each bin and store fit results
  TCanvas *c = new TCanvas("", "", 700, 700);    
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

    // get max costheta
    double cMaxVal = jumpF(cosMax->Eval(pMin));
    double cMinVal = jumpF(cosMin->Eval(pMax));

    // fit the 2 functions
    for(int i_t = 0; i_t < 2; i_t++) {
      fit1d[i_t]->SetRange(cMinVal, cMaxVal);
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
    pHist[0][i]->SetTitle(Form("#DeltaR>%.2f NP |cos#theta| (%.0f < p_{T} < %.0f GeV)", r_cut, pMin, pMax));
    pHist[0][i]->SetStats(0);
    pHist[0][i]->SetLineColor(kRed);
    pHist[0][i]->SetMarkerColor(kRed);
    pHist[0][i]->SetMinimum(0);
    //pHist[0][i]->SetMaximum(pHist[0][i]->GetBinContent(1)*1.5);
    pHist[0][i]->SetMaximum(parA[0][i]*1.5);
    if(i==nBinsY-1) pHist[0][i]->SetMaximum(pHist[0][i]->GetMaximum()*1.5);
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
    lc.DrawLatex(0.1, pHist[0][i]->GetMaximum()*0.9, Form("#lambda_{#theta}^{NP} = %.3f #pm %.3f", parL[0][i], eparL[0][i]));
    lc.DrawLatex(0.1, pHist[0][i]->GetMaximum()*0.8, Form("#lambda_{#theta}^{pure NP} = %.3f #pm %.3f", parL[1][i], eparL[1][i]));
    
    TLine *c_lim_Max = new TLine(cMaxVal, 0, cMaxVal, pHist[0][i]->GetMaximum());
    c_lim_Max->SetLineStyle(kDashed);
    c_lim_Max->SetLineColor(kBlack);
    c_lim_Max->Draw();
    TLine *c_lim_Min = new TLine(cMinVal, 0, cMinVal, pHist[0][i]->GetMaximum());
    c_lim_Min->SetLineStyle(kDashed);
    c_lim_Min->SetLineColor(kBlack);
    c_lim_Min->Draw();

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(pHist[0][i], "NP", "pl");
    leg->AddEntry(pHist[2][i], "SB contrib", "pl");
    leg->AddEntry(pHist[1][i], "pure NP", "pl");
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
