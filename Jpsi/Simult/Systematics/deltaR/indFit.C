#import "../../../Simult_dR1/cosMax/imp_jumpF.C"

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
  string loc = "/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi";

  // read the three coarse histos
  TFile *infile = new TFile("files/chistStore.root");
  TH2D **h_fit = new TH2D*[3];
  string lbl[] = {"B", "L", "T"};
  for(int i = 0; i < 3; i++) {
    infile->GetObject(Form("cHist%s", lbl[i].c_str()), h_fit[i]);
    h_fit[i]->SetDirectory(0);
  }
  infile->Close();

  // get the binning
  int nBinsX = h_fit[0]->GetNbinsX(), nBinsY = h_fit[0]->GetNbinsY();
  const double *yBins = h_fit[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *pHist[3][nBinsY];
  string t_name[3] = {"baseline", "#DeltaR>0.15", "#DeltaR>0.17"};
  for(int i_t = 0; i_t < 3; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      pHist[i_t][i-1] = h_fit[i_t]->ProjectionX(Form("bin%d_%d", i, i_t+1), i, i);
      pHist[i_t][i-1]->SetTitle(Form("%s bin %d: [%.0f, %.0f] GeV", t_name[i_t].c_str(), i, yBins[i-1], yBins[i]));
    }
  }
  
  // the fit function to be used
  TF1 **fit1d = new TF1*[3];
  for(int i = 0; i < 3; i++) {
    fit1d[i] = new TF1(Form("fit_%d", i), "[0]*(1+[1]*x*x)", 0, 1);
    fit1d[i]->SetParNames("A", "l_th");
  }
  
  // get the fit range from our cosmax(pT), cosmin(pT) - 3 different ranges
  ifstream inB, inT, inL;
  string dataS;
  double maxPar[3][3], aux;
  double minPar[3][4];
  // baseline
  inB.open(Form("%s/Simult/cosMax/cosMaxFitRes.txt", loc.c_str()));
  getline(inB, dataS);
  getline(inB, dataS);
  inB >> maxPar[0][0] >> aux >> maxPar[0][1] >> aux >> maxPar[0][2];
  inB.close();
  // loose cut
  inL.open(Form("%s/Simult_dR2/cosMax/cosMaxFitRes.txt", loc.c_str()));
  getline(inL, dataS);
  getline(inL, dataS);
  inL >> maxPar[1][0] >> aux >> maxPar[1][1] >> aux >> maxPar[1][2];
  inL.close();
  inL.open(Form("%s/Simult_dR2/cosMax/cosMinFitRes.txt", loc.c_str()));
  getline(inL, dataS);
  getline(inL, dataS);
  inL >> minPar[1][0] >> aux >> minPar[1][1] >> aux >> minPar[1][2] >> aux >> minPar[1][3];
  inL.close();
  // tight cut
  inT.open(Form("%s/Simult_dR1/cosMax/cosMaxFitRes.txt", loc.c_str()));
  getline(inT, dataS);
  getline(inT, dataS);
  inT >> maxPar[2][0] >> aux >> maxPar[2][1] >> aux >> maxPar[2][2];
  inT.close();
  inT.open(Form("%s/Simult_dR1/cosMax/cosMinFitRes.txt", loc.c_str()));
  getline(inT, dataS);
  getline(inT, dataS);
  inT >> minPar[2][0] >> aux >> minPar[2][1] >> aux >> minPar[2][2] >> aux >> minPar[2][3];
  inT.close();

  // define max and min functions
  // params to be set as the fit runs
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", yBins[0]-10, yBins[nBinsY]+10);
  TF1 *cosMin = new TF1("cosMin", "cminf(x, [0], [1], [2], [3])", yBins[0]-10, yBins[nBinsY]+10);
 
  // the cycle to fit each bin and store fit results
  TCanvas *c = new TCanvas("", "", 700, 700);    
  TFile *outfile = new TFile("files/finalFitRes.root", "recreate");

  double parA[3][nBinsY], eparA[3][nBinsY];
  double parL[3][nBinsY], eparL[3][nBinsY];
  double chi2[3][nBinsY], ndf[3][nBinsY], chiP[3][nBinsY];
  double pt[nBinsY], ept[nBinsY];
  double cMinVal[3], cMaxVal[3];
  
  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_fit[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_fit[0]->GetYaxis()->GetBinUpEdge(i+1);
    pt[i] = (pMax+pMin)/2.;
    ept[i] = (pMax-pMin)/2.;

    // fit the 3 cases
    for(int i_t = 0; i_t < 3; i_t++) {
      // get max costheta
      cosMax->SetParameters(maxPar[i_t][0], maxPar[i_t][1], maxPar[i_t][2]);
      cMaxVal[i_t] = jumpF(cosMax->Eval(pMin));

      if(i_t>0) {
	cosMin->SetParameters(minPar[i_t][0], minPar[i_t][1], minPar[i_t][2], minPar[i_t][3]);
	cMinVal[i_t] = jumpF(cosMin->Eval(pMax));
      }
      else {
	if(pt[i] < 66) {
	  cosMin->SetParameters(minPar[2][0], minPar[2][1], minPar[2][2], minPar[2][3]);
	  cMinVal[i_t] = jumpF(cosMin->Eval(pMax));
	}
	else 
	  
	  cMinVal[i_t] = 0;
      }
	
      fit1d[i_t]->SetRange(cMinVal[i_t], cMaxVal[i_t]);
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
    pHist[0][i]->SetTitle(Form("|cos#theta| (%.0f < p_{T} < %.0f GeV)", pMin, pMax));
    pHist[0][i]->SetStats(0);
    pHist[0][i]->SetLineColor(kBlack);
    pHist[0][i]->SetMarkerColor(kBlack);
    pHist[0][i]->SetMinimum(0);
    pHist[0][i]->SetMaximum(pHist[0][i]->GetBinContent(1)*1.5);
    pHist[0][i]->Draw("error");
    fit1d[0]->SetLineColor(kBlack);
    fit1d[0]->SetLineStyle(kDashed);
    fit1d[0]->Draw("same");

    pHist[1][i]->SetLineColor(kBlue);
    pHist[1][i]->SetMarkerColor(kBlue);
    pHist[1][i]->Draw("same");
    fit1d[1]->SetLineColor(kBlue);
    fit1d[1]->SetLineStyle(kDashed);
    fit1d[1]->Draw("same");

    pHist[2][i]->SetLineColor(kRed);
    pHist[2][i]->SetMarkerColor(kRed);
    pHist[2][i]->Draw("same");
    fit1d[2]->SetLineColor(kRed);
    fit1d[2]->SetLineStyle(kDashed);
    fit1d[2]->Draw("same");

    TLatex lc;
    lc.SetTextSize(0.03);
    lc.DrawLatex(0.1, pHist[0][i]->GetMaximum()*0.9, Form("#lambda_{#theta}^{base} = %.3f #pm %.3f", parL[0][i], eparL[0][i]));
    lc.DrawLatex(0.1, pHist[0][i]->GetMaximum()*0.8, Form("#lambda_{#theta}^{loose} = %.3f #pm %.3f", parL[1][i], eparL[1][i]));
    lc.DrawLatex(0.1, pHist[0][i]->GetMaximum()*0.7, Form("#lambda_{#theta}^{tight} = %.3f #pm %.3f", parL[2][i], eparL[2][i]));

    int col[3] = {kBlack, kBlue, kRed};
    for(int i_t = 0; i_t < 3; i_t++) {
      TLine *c_lim_Max = new TLine(cMaxVal[i_t], 0, cMaxVal[i_t], pHist[0][i]->GetMaximum());
      c_lim_Max->SetLineStyle(kDashed);
      c_lim_Max->SetLineColor(col[i_t]);
      c_lim_Max->Draw();
      TLine *c_lim_Min = new TLine(cMinVal[i_t], 0, cMinVal[i_t], pHist[0][i]->GetMaximum());
      c_lim_Min->SetLineStyle(kDashed);
      c_lim_Min->SetLineColor(col[i_t]);
      c_lim_Min->Draw();
      cout << endl << i_t << " " << cMinVal[i_t] << " " << cMaxVal[i_t] << endl;
    }
    cout << endl << endl << cMinVal[2] << endl << endl;
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(pHist[0][i], "baseline", "pl");
    leg->AddEntry(pHist[1][i], "#DeltaR>0.15", "pl");
    leg->AddEntry(pHist[2][i], "#DeltaR>0.17", "pl");
    leg->Draw();

    for(int i_t = 0; i_t < 3; i_t++)
      pHist[i_t][i]->Write();
 
    c->SaveAs(Form("plots/ratioFinal/bin_%d.pdf", i));
    c->Clear();
    cout << endl << endl;
  }

  for(int i_t = 0; i_t < 3; i_t++) {
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
