#import "../cosMax/imp_jumpF.C"

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
      pHist[i_t][i-1]->SetTitle(Form("%s bin %d: [%.0f, %.0f] GeV", lbl[i_t].c_str(), i, yBins[i-1], yBins[i]));
    }
  }
  
  // the fit function to be used - only on total and prompt J/psi
  TF1 **fit1d = new TF1*[4];
  for(int i = 0; i < 4; i++) {
    fit1d[i] = new TF1(Form("fit_%d", i), "[0]*(1+[1]*x*x)", 0, 1);
    fit1d[i]->SetParNames("A", "l_th");
  }
  
  // get the fit range from our cosmax(pT)
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
 
  // the cycle to fit each bin and store fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.01);
  TFile *outfile = new TFile("files/finalFitRes.root", "recreate");

  double parA[4][nBinsY], eparA[4][nBinsY];
  double parL[4][nBinsY], eparL[4][nBinsY];
  double chi2[4][nBinsY], ndf[4][nBinsY], chiP[4][nBinsY];
  double cMax[nBinsY], pt[nBinsY], ept[nBinsY], zero[nBinsY];
  
  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_fit[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_fit[0]->GetYaxis()->GetBinUpEdge(i+1);
    pt[i] = (pMax+pMin)/2.;
    ept[i] = (pMax-pMin)/2.;
    zero[i] = 0;

    // get max costheta
    double cMaxVal = jumpF(cosMax->Integral(pMin, pMax)/(pMax-pMin))-0.05;
    cMax[i] = cMaxVal;

    // fit the 4 functions
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
    fit1d[3]->SetRange(0, cMaxVal);
    fit1d[3]->SetParameters(pHist[2][i]->GetBinContent(1)*1.1, 0.1);

    pHist[2][i]->Fit(fit1d[3], "R0");

    parA[3][i] = fit1d[3]->GetParameter(0);
    eparA[3][i] = fit1d[3]->GetParError(0);
    parL[3][i] = fit1d[3]->GetParameter(1);
    eparL[3][i] = fit1d[3]->GetParError(1);
    chi2[3][i] = fit1d[3]->GetChisquare();
    ndf[3][i] = fit1d[3]->GetNDF();
    chiP[3][i] = TMath::Prob(chi2[3][i], ndf[3][i]);

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
    lc.DrawLatex(0.05, pHist[0][i]->GetMaximum()*0.9, Form("#lambda_{#theta}^{total} = %.3f #pm %.3f", parL[0][i], eparL[0][i]));
    lc.DrawLatex(0.05, pHist[0][i]->GetMaximum()*0.8, Form("#lambda_{#theta}^{prompt #psi(2S)} = %.3f #pm %.3f", parL[1][i], eparL[1][i]));
    
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
 
    c->SaveAs(Form("plots/ratioFinal/fit/bin_%d.pdf", i));
    c->Clear();
  }
  outfile->Close();

  // calculating pulls - just prompt and non-prompt psi
  // NP needs to come from the corresponding input
  TFile *infile_NP = new TFile("../NP_fit/files/bkgSubRes.root");
  infile_NP->GetObject(Form("h_NPc"), h_fit[1]);
  h_fit[1]->SetDirectory(0);
  infile_NP->Close();

  // get the 1d plots
  for(int i_pt = 1; i_pt <= nBinsY; i_pt++) {
    pHist[1][i_pt-1] = h_fit[1]->ProjectionX(Form("bin%d_%d", i_pt, 2), i_pt, i_pt);
    pHist[1][i_pt-1]->SetTitle(Form("NP bin %d: [%.0f, %.0f] GeV", i_pt, yBins[i_pt-1], yBins[i_pt]));
  }
  
  double xv[nBinsX], pv_NP[nBinsX], dv_NP[nBinsX], pv_PR[nBinsX], dv_PR[nBinsX];
  // new cycle, now for PR and NP fit plots
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

    // need to restore PR fit parameters
    fit1d[1]->SetRange(0, cMax[i]);
    fit1d[1]->SetParameters(parA[1][i], parL[1][i]);
    
    for(int i_cos = 0 ; i_cos < nBinsX; i_cos++) {
      xv[i_cos] = pHist[0][i]->GetBinCenter(i_cos+1);
      
      // first non-prompt 
      double fitv = fit1d[2]->Eval(xv[i_cos]);
      double datav = pHist[1][i]->GetBinContent(i_cos+1);
      double datau = pHist[1][i]->GetBinError(i_cos+1);
      if(xv[i_cos] < cMax[i]) {
	pv_NP[i_cos] = (datav-fitv)/datau;
	dv_NP[i_cos] = (datav-fitv)/fitv * 100.;
      }
      else {
	pv_NP[i_cos] = 0;
	dv_NP[i_cos] = 0;
      }
      
      // then prompt
      fitv = fit1d[1]->Eval(xv[i_cos]);
      datav = pHist[3][i]->GetBinContent(i_cos+1);
      datau = pHist[3][i]->GetBinError(i_cos+1);
      if(xv[i_cos] < cMax[i]) {
	pv_PR[i_cos] = (datav-fitv)/datau;
	dv_PR[i_cos] = (datav-fitv)/fitv * 100.;
      }
      else {
	pv_PR[i_cos] = 0;
	dv_PR[i_cos] = 0;
      }
    }

    // now plotting the NP
    // plotting the pulls
    TH1F *fl = c->DrawFrame(0, -9, 1, 9);
    fl->SetXTitle("|cos #theta_{HX}|");
    fl->SetYTitle("pulls");
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("Non-prompt |cos #theta_{HX}| fit pulls (%.1f < p_{T} < %.1f GeV)", pMin, pMax));

    TGraph *gNP_pull = new TGraph(nBinsX, xv, pv_NP);
    for(int i_cos = nBinsX-1; i_cos > 0; i_cos--) {
      if(xv[i_cos] > cMax[i]) gNP_pull->RemovePoint(i_cos);
    }
    gNP_pull->SetLineColor(kBlack);
    gNP_pull->SetMarkerColor(kBlack);
    gNP_pull->SetMarkerStyle(20);
    gNP_pull->Draw("p");
    
    TLine *zero = new TLine(0, 0, 1, 0);
    zero->SetLineStyle(kDashed);
    zero->Draw();

    TLine *plim1 = new TLine(0, -5, 1, -5);
    plim1->SetLineStyle(kDotted);
    plim1->Draw("lsame");
    TLine *plim2 = new TLine(0, -3, 1, -3);
    plim2->SetLineStyle(kDotted);
    plim2->Draw("lsame");
    TLine *plim3 = new TLine(0, 3, 1, 3);
    plim3->SetLineStyle(kDotted);
    plim3->Draw("lsame");
    TLine *plim4 = new TLine(0, 5, 1, 5);
    plim4->SetLineStyle(kDotted);
    plim4->Draw("lsame");
    
    c->SaveAs(Form("plots/ratioFinal/fit/pullsNP_pt%d.pdf", i));
    c->Clear();

    // plotting the devs
    TH1F *fd = c->DrawFrame(0, -15, 1, 15);
    fd->SetXTitle("|cos #theta_{HX}|");
    fd->SetYTitle("relative difference (%)");
    fd->GetYaxis()->SetTitleOffset(1.3);
    fd->GetYaxis()->SetLabelOffset(0.01);
    fd->SetTitle(Form("Non-prompt |cos #theta_{HX}| rel. difference (%.1f < p_{T} < %.1f GeV)",  pMin, pMax));
  
    TGraph *gNP_dev = new TGraph(nBinsX, xv, dv_NP);
    for(int i_cos = nBinsX-1; i_cos > 0; i_cos--) {
      if(xv[i_cos] > cMax[i]) gNP_dev->RemovePoint(i_cos);
    }
    gNP_dev->SetLineColor(kBlack);		
    gNP_dev->SetMarkerColor(kBlack);
    gNP_dev->SetMarkerStyle(20);
    gNP_dev->Draw("psame");
    
    // aux lines - pull = 0 and sigma limits
    zero->Draw("lsame");

    c->SaveAs(Form("plots/ratioFinal/fit/devsNP_pt%d.pdf", i));
    c->Clear();

    // now plotting the PR
    // plotting the pulls
    TH1F *flP = c->DrawFrame(0, -9, 1, 9);
    flP->SetXTitle("|cos #theta_{HX}|");
    flP->SetYTitle("pulls");
    flP->GetYaxis()->SetTitleOffset(1.3);
    flP->GetYaxis()->SetLabelOffset(0.01);
    flP->SetTitle(Form("Prompt |cos #theta_{HX}| fit pulls (%.1f < p_{T} < %.1f GeV)", pMin, pMax));

    TGraph *gPR_pull = new TGraph(nBinsX, xv, pv_PR);
    for(int i_cos = nBinsX-1; i_cos > 0; i_cos--) {
      if(xv[i_cos] > cMax[i]) gPR_pull->RemovePoint(i_cos);
    }
    gPR_pull->SetLineColor(kBlack);
    gPR_pull->SetMarkerColor(kBlack);
    gPR_pull->SetMarkerStyle(20);
    gPR_pull->Draw("p");
    
    zero->Draw();

    plim1->Draw("lsame");
    plim2->Draw("lsame");
    plim3->Draw("lsame");
    plim4->Draw("lsame");
    
    c->SaveAs(Form("plots/ratioFinal/fit/pullsPR_pt%d.pdf", i));
    c->Clear();

    // plotting the devs
    TH1F *fdP = c->DrawFrame(0, -15, 1, 15);
    fdP->SetXTitle("|cos #theta_{HX}|");
    fdP->SetYTitle("relative difference (%)");
    fdP->GetYaxis()->SetTitleOffset(1.3);
    fdP->GetYaxis()->SetLabelOffset(0.01);
    fdP->SetTitle(Form("Prompt |cos #theta_{HX}| rel. difference (%.1f < p_{T} < %.1f GeV)",  pMin, pMax));
  
    TGraph *gPR_dev = new TGraph(nBinsX, xv, dv_PR);
    for(int i_cos = nBinsX-1; i_cos > 0; i_cos--) {
      if(xv[i_cos] > cMax[i]) gPR_dev->RemovePoint(i_cos);
    }
    gPR_dev->SetLineColor(kBlack);
    gPR_dev->SetMarkerColor(kBlack);
    gPR_dev->SetMarkerStyle(20);
    gPR_dev->Draw("psame");
    
    // aux lines - pull = 0 and sigma limits
    zero->Draw("lsame");

    c->SaveAs(Form("plots/ratioFinal/fit/devsPR_pt%d.pdf", i));
    c->Clear();

    // also plotting just the prompt and just the non-prompt J/psi
    // plotting NP
    pHist[1][i]->SetTitle(Form("Non-prompt |cos #theta_{HX}| (%.1f < p_{T} < %.1f GeV)", pMin, pMax));
    pHist[1][i]->SetStats(0);
    pHist[1][i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    pHist[1][i]->SetMinimum(0);
    pHist[1][i]->SetMaximum(pHist[1][i]->GetMaximum()*1.1);
    pHist[1][i]->SetLineColor(kRed);
    pHist[1][i]->SetMarkerColor(kRed);
    pHist[1][i]->Draw("error");
    fit1d[2]->SetLineColor(kRed);
    fit1d[2]->SetLineStyle(kDashed);
    fit1d[2]->Draw("same");

    TLine *c_limNP = new TLine(cMax[i], 0, cMax[i], pHist[1][i]->GetMaximum());
    c_limNP->SetLineStyle(kDashed);
    c_limNP->SetLineColor(kBlack);
    c_limNP->Draw();
 
    c->SaveAs(Form("plots/ratioFinal/fit/fitNP_%d.pdf", i));
    c->Clear();

    // plotting PR
    pHist[3][i]->SetTitle(Form("Prompt |cos #theta_{HX}| (%.1f < p_{T} < %.1f GeV)", pMin, pMax));
    pHist[3][i]->SetStats(0);
    pHist[3][i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    pHist[3][i]->SetMinimum(0);
    pHist[3][i]->SetMaximum(pHist[3][i]->GetMaximum()*1.1);
    pHist[3][i]->SetLineColor(kBlue);
    pHist[3][i]->SetMarkerColor(kBlue);
    pHist[3][i]->Draw("error");
    fit1d[1]->SetLineColor(kBlue);
    fit1d[1]->SetLineStyle(kDashed);
    fit1d[1]->Draw("same");

    TLine *c_limPR = new TLine(cMax[i], 0, cMax[i], pHist[3][i]->GetMaximum());
    c_limPR->SetLineStyle(kDashed);
    c_limPR->SetLineColor(kBlack);
    c_limPR->Draw();
    
    c->SaveAs(Form("plots/ratioFinal/fit/fitPR_%d.pdf", i));
    c->Clear();
  }

  string lbl_s[] = {"Data", "J", "NP", "PR"};
  TFile *outfile2 = new TFile("files/finalFitRes.root", "update");
    for(int i_t = 0; i_t < 4; i_t++) {
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
  TGraphErrors *graphCm = new TGraphErrors(nBinsY, pt, cMax, ept, zero);
  graphCm->SetName(Form("graph_cMax"));
  graphCm->Write();

  outfile2->Close();
  
  c->Destructor();
}
