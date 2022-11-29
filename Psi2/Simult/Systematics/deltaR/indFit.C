#import "../../../Simult_dR1/cosMax/imp_jumpF.C"

// code to do the individual fit (1d costheta maps)
// on the fine bins (using same fit |costheta| region)

// aux func for costheta_min
double cminf(double pt, double a, double b, double c)
{
 double val = 0;
 
 if (pt < c) return 0;
 else
   val = a + b * pt;
 if(val > 0) return val;
 else return 0;
}

// main
void indFit()
{
  double pt_c = 50;
  string loc = "/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2";

  // PART 1 - FINE BIN FITS
  // read the three fine histos
  TH2D **h_fit_f = new TH2D*[3];
  
  TFile *infile_f0 = new TFile(Form("%s/Simult/PR_fit/files/bkgSubRes.root", loc.c_str()));
  h_fit_f[0] = (TH2D*)infile_f0->Get("h_J");
  h_fit_f[0]->SetName("h_JB");
  h_fit_f[0]->SetDirectory(0);
  infile_f0->Close();

  TFile *infile_f1 = new TFile(Form("%s/Simult_dR2/PR_fit/files/bkgSubRes.root", loc.c_str()));
  h_fit_f[1] = (TH2D*)infile_f1->Get("h_J");
  h_fit_f[1]->SetName("h_JL");
  h_fit_f[1]->SetDirectory(0);
  infile_f1->Close();

  TFile *infile_f2 = new TFile(Form("%s/Simult_dR1/PR_fit/files/bkgSubRes.root", loc.c_str()));
  h_fit_f[2] = (TH2D*)infile_f2->Get("h_J");
  h_fit_f[2]->SetName("h_JT");
  h_fit_f[2]->SetDirectory(0);
  infile_f2->Close();

  // get the binning
  int nBinsX = h_fit_f[0]->GetNbinsX();
  int nBinsY = h_fit_f[0]->GetNbinsY();
  const double *yBins_f = h_fit_f[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *pHist_f[3][nBinsY];
  string t_name[3] = {"baseline", "#DeltaR>0.15", "#DeltaR>0.17"};
  for(int i_t = 0; i_t < 3; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      pHist_f[i_t][i-1] = h_fit_f[i_t]->ProjectionX(Form("bin%d_%d", i, i_t+1), i, i);
      pHist_f[i_t][i-1]->SetTitle(Form("%s bin %d: [%.1f, %.1f] GeV", t_name[i_t].c_str(), i, yBins_f[i-1], yBins_f[i]));
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
  double minPar[3][3];
  // loose cut
  inL.open(Form("%s/Simult_dR2/cosMax/cosMaxFitRes.txt", loc.c_str()));
  getline(inL, dataS);
  getline(inL, dataS);
  inL >> maxPar[1][0] >> aux >> maxPar[1][1] >> aux >> maxPar[1][2];
  inL.close();
  inL.open(Form("%s/Simult_dR2/cosMax/cosMinFitRes.txt", loc.c_str()));
  getline(inL, dataS);
  getline(inL, dataS);
  inL >> minPar[1][0] >> aux >> minPar[1][1] >> aux >> minPar[1][2];
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
  inT >> minPar[2][0] >> aux >> minPar[2][1] >> aux >> minPar[2][2];
  inT.close();

  // define max and min functions
  // params to be set as the fit runs
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", yBins_f[0]-10, yBins_f[nBinsY]+10);
  TF1 *cosMin = new TF1("cosMin", "cminf(x, [0], [1], [2])", yBins_f[0]-10, yBins_f[nBinsY]+10);

  // the cycle to fit each bin and store fit results
  TCanvas *c = new TCanvas("", "", 700, 700);    
  TFile *outfile_f = new TFile("files/finalFitRes.root", "recreate");

  double parA_f[3][nBinsY], eparA_f[3][nBinsY];
  double parL_f[3][nBinsY], eparL_f[3][nBinsY];
  double chi2_f[3][nBinsY], ndf_f[3][nBinsY], chiP_f[3][nBinsY];
  double pt_f[nBinsY], ept_f[nBinsY];
  double cMinVal[3], cMaxVal[3];

  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_fit_f[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_fit_f[0]->GetYaxis()->GetBinUpEdge(i+1);
    pt_f[i] = (pMax+pMin)/2.;
    ept_f[i] = (pMax-pMin)/2.;

    // fit the 3 cases
    for(int i_t = 0; i_t < 3; i_t++) {
      // get max costheta
      if(i_t>0) {
	cosMax->SetParameters(maxPar[i_t][0], maxPar[i_t][1], maxPar[i_t][2]);
	cMaxVal[i_t] = jumpF(cosMax->Eval(pMin));
	cosMin->SetParameters(minPar[i_t][0], minPar[i_t][1], minPar[i_t][2]);
	cMinVal[i_t] = jumpF(cosMin->Eval(pMax));
      }
      else {
	if(pt_f[i] < pt_c) {
	  cosMax->SetParameters(maxPar[2][0], maxPar[2][1], maxPar[2][2]);
	  cMaxVal[i_t] = jumpF(cosMax->Eval(pMin));
	  cosMin->SetParameters(minPar[2][0], minPar[2][1], minPar[2][2]);
	  cMinVal[i_t] = jumpF(cosMin->Eval(pMax));
	}
	else {
	  cosMax->SetParameters(maxPar[1][0], maxPar[1][1], maxPar[1][2]);
	  cMaxVal[i_t] = jumpF(cosMax->Eval(pMin));
	  cosMin->SetParameters(minPar[1][0], minPar[1][1], minPar[1][2]);
	  cMinVal[i_t] = jumpF(cosMin->Eval(pMax));
	}
      }
	
      fit1d[i_t]->SetRange(cMinVal[i_t], cMaxVal[i_t]);
      fit1d[i_t]->SetParameters(pHist_f[i_t][i]->GetBinContent(1)*1.1, 0.1);

      pHist_f[i_t][i]->Fit(fit1d[i_t], "R0");

      parA_f[i_t][i] = fit1d[i_t]->GetParameter(0);
      eparA_f[i_t][i] = fit1d[i_t]->GetParError(0);
      parL_f[i_t][i] = fit1d[i_t]->GetParameter(1);
      eparL_f[i_t][i] = fit1d[i_t]->GetParError(1);
      chi2_f[i_t][i] = fit1d[i_t]->GetChisquare();
      ndf_f[i_t][i] = fit1d[i_t]->GetNDF();
      chiP_f[i_t][i] = TMath::Prob(chi2_f[i_t][i], ndf_f[i_t][i]);
    }

    // plotting everything
    pHist_f[0][i]->SetTitle(Form("|cos#theta| (%.1f < p_{T} < %.1f GeV)", pMin, pMax));
    pHist_f[0][i]->SetStats(0);
    pHist_f[0][i]->SetLineColor(kBlack);
    pHist_f[0][i]->SetMarkerColor(kBlack);
    pHist_f[0][i]->SetMinimum(0);
    pHist_f[0][i]->SetMaximum(pHist_f[0][i]->GetBinContent(1)*1.5);
    pHist_f[0][i]->Draw("error");
    fit1d[0]->SetLineColor(kBlack);
    fit1d[0]->SetLineStyle(kDashed);
    fit1d[0]->Draw("same");

    pHist_f[1][i]->SetLineColor(kBlue);
    pHist_f[1][i]->SetMarkerColor(kBlue);
    pHist_f[1][i]->Draw("same");
    fit1d[1]->SetLineColor(kBlue);
    fit1d[1]->SetLineStyle(kDashed);
    fit1d[1]->Draw("same");

    pHist_f[2][i]->SetLineColor(kRed);
    pHist_f[2][i]->SetMarkerColor(kRed);
    pHist_f[2][i]->Draw("same");
    fit1d[2]->SetLineColor(kRed);
    fit1d[2]->SetLineStyle(kDashed);
    fit1d[2]->Draw("same");

    TLatex lc;
    lc.SetTextSize(0.03);
    lc.DrawLatex(0.1, pHist_f[0][i]->GetMaximum()*0.9, Form("#lambda_{#theta}^{base} = %.3f #pm %.3f", parL_f[0][i], eparL_f[0][i]));
    lc.DrawLatex(0.1, pHist_f[0][i]->GetMaximum()*0.8, Form("#lambda_{#theta}^{loose} = %.3f #pm %.3f", parL_f[1][i], eparL_f[1][i]));
    lc.DrawLatex(0.1, pHist_f[0][i]->GetMaximum()*0.7, Form("#lambda_{#theta}^{tight} = %.3f #pm %.3f", parL_f[2][i], eparL_f[2][i]));

    int col[3] = {kBlack, kBlue, kRed};
    for(int i_t = 0; i_t < 3; i_t++) {
      TLine *c_lim_Max = new TLine(cMaxVal[i_t], 0, cMaxVal[i_t], pHist_f[0][i]->GetMaximum());
      c_lim_Max->SetLineStyle(kDashed);
      c_lim_Max->SetLineColor(col[i_t]);
      c_lim_Max->Draw();
      TLine *c_lim_Min = new TLine(cMinVal[i_t], 0, cMinVal[i_t], pHist_f[0][i]->GetMaximum());
      c_lim_Min->SetLineStyle(kDashed);
      c_lim_Min->SetLineColor(col[i_t]);
      c_lim_Min->Draw();
      cout << endl << i_t << " " << cMinVal[i_t] << " " << cMaxVal[i_t] << endl;
    }
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(pHist_f[0][i], "baseline", "pl");
    leg->AddEntry(pHist_f[1][i], "#DeltaR>0.15", "pl");
    leg->AddEntry(pHist_f[2][i], "#DeltaR>0.17", "pl");
    leg->Draw();

    for(int i_t = 0; i_t < 3; i_t++)
      pHist_f[i_t][i]->Write();
 
    c->SaveAs(Form("plots/ratioFinal/bin_%d.pdf", i));
    c->Clear();
    cout << endl << endl;
  }

  string lbl[] = {"B", "L", "T"};
  for(int i_t = 0; i_t < 3; i_t++) {
    // make and save the TGraph with the fit results and max costh used
    TGraphErrors *graphA = new TGraphErrors(nBinsY, pt_f, parA_f[i_t], ept_f, eparA_f[i_t]);
    TGraphErrors *graphL = new TGraphErrors(nBinsY, pt_f, parL_f[i_t], ept_f, eparL_f[i_t]);
    TGraph *graphC = new TGraph(nBinsY, pt_f, chi2_f[i_t]);
    TGraph *graphN = new TGraph(nBinsY, pt_f, ndf_f[i_t]);
    TGraph *graphP = new TGraph(nBinsY, pt_f, chiP_f[i_t]);

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
  outfile_f->Close();

  // PART 2: the coarse-binned results
  
  // read the three coarse histos
  /*  TFile *infile = new TFile("files/chistStore.root");
  TH2D **h_fit = new TH2D*[3];
  for(int i = 0; i < 3; i++) {
    infile->GetObject(Form("cHist%s", lbl[i].c_str()), h_fit[i]);
    h_fit[i]->SetDirectory(0);
  }
  infile->Close();

  // get the binning
  nBinsX = h_fit[0]->GetNbinsX();
  nBinsY = h_fit[0]->GetNbinsY();
  const double *yBins = h_fit[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *pHist[3][nBinsY];
  for(int i_t = 0; i_t < 3; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      pHist[i_t][i-1] = h_fit[i_t]->ProjectionX(Form("bin%d_%d", i, i_t+1), i, i);
      pHist[i_t][i-1]->SetTitle(Form("%s bin %d: [%.1f, %.1f] GeV", t_name[i_t].c_str(), i, yBins[i-1], yBins[i]));
    }
  }
   
  // the cycle to fit each bin and store fit results
  TFile *outfile = new TFile("files/finalFitRes.root", "update");

  double parA[3][nBinsY], eparA[3][nBinsY];
  double parL[3][nBinsY], eparL[3][nBinsY];
  double chi2[3][nBinsY], ndf[3][nBinsY], chiP[3][nBinsY];
  double pt[nBinsY], ept[nBinsY];
  
  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_fit[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_fit[0]->GetYaxis()->GetBinUpEdge(i+1);
    pt[i] = (pMax+pMin)/2.;
    ept[i] = (pMax-pMin)/2.;

    // fit the 3 cases
    for(int i_t = 0; i_t < 3; i_t++) {
      // get max costheta
      if(i_t>0) {
	cosMax->SetParameters(maxPar[i_t][0], maxPar[i_t][1], maxPar[i_t][2]);
	cMaxVal[i_t] = jumpF(cosMax->Eval(pMin));
	cosMin->SetParameters(minPar[i_t][0], minPar[i_t][1], minPar[i_t][2], minPar[i_t][3]);
	cMinVal[i_t] = jumpF(cosMin->Eval(pMax));
      }
      else {
	if(pt[i] < pt_c) {
	  cosMax->SetParameters(maxPar[2][0], maxPar[2][1], maxPar[2][2]);
	  cMaxVal[i_t] = jumpF(cosMax->Eval(pMin));
	  cosMin->SetParameters(minPar[2][0], minPar[2][1], minPar[2][2], minPar[2][3]);
	  cMinVal[i_t] = jumpF(cosMin->Eval(pMax));
	}
	else {
	  cosMax->SetParameters(maxPar[1][0], maxPar[1][1], maxPar[1][2]);
	  cMaxVal[i_t] = jumpF(cosMax->Eval(pMin));
	  cosMin->SetParameters(minPar[1][0], minPar[1][1], minPar[1][2], minPar[1][3]);
	  cMinVal[i_t] = jumpF(cosMin->Eval(pMax));
	}
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
    pHist[0][i]->SetTitle(Form("|cos#theta| (%.1f < p_{T} < %.1f GeV)", pMin, pMax));
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
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(pHist[0][i], "baseline", "pl");
    leg->AddEntry(pHist[1][i], "#DeltaR>0.15", "pl");
    leg->AddEntry(pHist[2][i], "#DeltaR>0.17", "pl");
    leg->Draw();

    for(int i_t = 0; i_t < 3; i_t++)
      pHist[i_t][i]->Write();
 
    c->SaveAs(Form("plots/ratioFinal/bin_c_%d.pdf", i));
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

    graphA->SetName(Form("graph_A_c_%s", lbl[i_t].c_str()));
    graphL->SetName(Form("graph_lambda_c_%s", lbl[i_t].c_str()));
    graphC->SetName(Form("graph_chisquare_c_%s", lbl[i_t].c_str()));
    graphN->SetName(Form("graph_NDF_c_%s", lbl[i_t].c_str()));
    graphP->SetName(Form("graph_chiP_c_%s", lbl[i_t].c_str()));
 
    graphA->Write();
    graphL->Write();
    graphC->Write();
    graphN->Write();
    graphP->Write();
  }
  outfile->Close();*/

  c->Destructor();
}
