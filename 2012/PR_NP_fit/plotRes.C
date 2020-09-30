// code to combine the free fit function with the ratio points
// also plots the lambda_th function

void plotRes()
{
  double M_q = 1;//3.097;
  
  // get the ratio points
  TFile *fIn = new TFile("files/ratioHist.root");
  TH2D* rHist;
  fIn->GetObject("ratioHist_ab", rHist);
  rHist->SetDirectory(0);
  fIn->Close();

  // save some useful info
  int nBinspT = rHist->GetNbinsY(), nBinsC = rHist->GetNbinsX();
  double minC = rHist->GetXaxis()->GetBinLowEdge(1), maxC = rHist->GetXaxis()->GetBinUpEdge(nBinsC);
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();

  // get the ratios in each coarse pT bin
  TH1D *hData[nBinspT];
  for(int i = 0; i < nBinspT; i++)
    hData[i] = rHist->ProjectionX(Form("hData_%d", i), i+1, i+1);

  // get the independent fit points
  // get A and lambda values for each bin
  TFile *fInd = new TFile("files/fit_res_1d.root");
  TGraphErrors* graph_A = (TGraphErrors*)fInd->Get("graph_A");
  TGraphErrors* graph_lth;
  fInd->GetObject("graph_lambda", graph_lth);
  fInd->Close();

  double *lambda = graph_lth->GetY();
  double *A = graph_A->GetY();

  // get the fit function, plot in histograms
  TF1 *wInd = new TF1("wInd", "[0]*(1+[1]*x*x)", minC, maxC);
  TH1D **hInd = new TH1D*[nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    hInd[i] = new TH1D(Form("hInd_%d", i), Form("hInd_%d", i), nBinsC, minC, maxC);
    wInd->SetParameters(A[i], lambda[i]);
    for(int j = 0; j < nBinsC; j++) {
      double min = hInd[i]->GetXaxis()->GetBinLowEdge(j+1);
      double max = hInd[i]->GetXaxis()->GetBinUpEdge(j+1);
      double fillVal = wInd->Integral(min, max)/(max-min);
      hInd[i]->SetBinContent(j+1, fillVal);
    }
  }

  // read the |costh|max(pt) function
  ifstream in;
  string dataS;
  in.open("text_output/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x/[3])", 0, 75);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2], 3.097);
  
  // plot the results
  TCanvas *c = new TCanvas("name", "title", 700, 700);
  for(int i = 0; i < nBinspT; i++) {
    // get the costh fit limit
    double cosLim = cosMax->Integral(pTBins[i], pTBins[i+1])/(pTBins[i+1]-pTBins[i]);
    double cR = floor(cosLim*10.)/10.;
    if(cosLim-cR>0.05) cR += 0.05;
    cout << cR << endl << endl;;

    // the data
    for(int j = 0; j < nBinsC; j++) {
      if(hData[i]->GetXaxis()->GetBinUpEdge(j+1)>cR+0.01) {
	cout << hData[i]->GetXaxis()->GetBinUpEdge(j+1) << endl;
	hData[i]->SetBinContent(j+1, 0);
	hData[i]->SetBinError(j+1, 0);
      }
    }
    
    hData[i]->SetStats(0);
    hData[i]->SetTitle(Form("PR/NP p_{T} bin %d: [%.0f, %.0f] GeV", i+1, pTBins[i]*M_q, pTBins[i+1]*M_q));
    hData[i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    hData[i]->SetLineColor(kBlack);
    hData[i]->SetMarkerColor(kBlack);
    hData[i]->GetYaxis()->SetRangeUser(0, A[i]*1.4);
    hData[i]->Draw();

    // the independent fit
    hInd[i]->SetLineColor(kBlue);
    hInd[i]->SetMarkerColor(kBlue);
    hInd[i]->Draw("same");

    TLine *cLim = new TLine(cR, 0, cR, A[i]*1.4);
    cLim->SetLineStyle(kDashed);
    cLim->SetLineColor(kBlack);
    cLim->Draw();

    // the legend
    TLegend *leg = new TLegend(0.2, 0.2, 0.4, 0.4);
    leg->SetTextSize(0.03);
    leg->AddEntry(hData[i], "data", "pl");
    leg->AddEntry(hInd[i], "fit", "pl");
    leg->Draw();

    // save and clear
    c->SaveAs(Form("plots/ratio_final/bin_%d.pdf", i));
    c->Clear();
  }

  // plot the results normalized
  for(int i = 0; i < nBinspT; i++) {
    // the data
    hData[i]->SetTitle(Form("PR/NP p_{T} bin %d: [%.0f, %.0f] GeV (scaled)", i+1, pTBins[i]*M_q, pTBins[i+1]*M_q));
    double norm = hData[i]->GetBinContent(1);
    hData[i]->Scale(1./norm);
    hData[i]->GetYaxis()->SetRangeUser(0.7, 1.3);
    hData[i]->Draw();

    // the independent fit
    norm = hInd[i]->GetBinContent(1);
    hInd[i]->Scale(1./norm);
    hInd[i]->Draw("hist same");

    double cosLim = cosMax->Integral(pTBins[i], pTBins[i+1])/(pTBins[i+1]-pTBins[i]);
    double cR = floor(cosLim*10.)/10.;
    if(cosLim-cR>0.05) cR += 0.05;
    TLine *cLim = new TLine(cR, 0.7, cR, 1.3);
    cLim->SetLineStyle(kDashed);
    cLim->SetLineColor(kBlack);
    cLim->Draw();

    // the legend
    TLegend *leg = new TLegend(0.2, 0.2, 0.35, 0.35);
    leg->SetTextSize(0.03);
    leg->AddEntry(hData[i], "data", "pl");
    leg->AddEntry(hInd[i], "fit", "pl");
    leg->Draw();

    // save and clear
    c->SaveAs(Form("plots/ratio_final/bin_%d_norm.pdf", i));
    c->Clear();
  }

  // draw lambda_th(pT)
  // the frame
  TH1F *func = c->DrawFrame(0, -1, 80, 1);
  func->SetXTitle("p_{T} (GeV)");
  func->SetYTitle("#lambda_{#theta}");
  func->GetYaxis()->SetTitleOffset(1.3);
  func->GetYaxis()->SetLabelOffset(0.01);
  func->SetTitle("#lambda_{#theta} (PR/NP)");

  // the three lambda_th distributions
  graph_lth->SetLineColor(kBlack);
  graph_lth->SetMarkerColor(kBlack);
  graph_lth->Draw("p");

  TLine *zero = new TLine(0, 0, 80, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  
  c->SaveAs("plots/ratio_final/lth.pdf");
  c->Clear();
  c->Destructor();
  
}

