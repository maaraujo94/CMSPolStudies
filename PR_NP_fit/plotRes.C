// code to combine all the fit functions and plot with the ratio points
// also plots the lambda_th functions all superimposed

void plotRes()
{
  double M_q = 3.097;
  
  // get the ratio points
  TFile *fIn = new TFile("files/ratioHist.root");
  TH2D* rHist;
  fIn->GetObject("cHist_ab", rHist);
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

  // get the linear/constant lambda fit points
  TFile *fDep = new TFile("files/fit_res_2d.root");
  TH1D *fitLin[nBinspT], *fitCon[nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    fitLin[i] = (TH1D*)fDep->Get(Form("Fit_linear_%d", i));
    fitCon[i] = (TH1D*)fDep->Get(Form("Fit_constant_%d", i));
  }
  // get the lambda_th
  TGraphErrors *lin_lth = (TGraphErrors*)fDep->Get("lth_linear");
  TGraphErrors *con_lth = (TGraphErrors*)fDep->Get("lth_constant");
  
  // plot the results
  TCanvas *c = new TCanvas("name", "title", 700, 700);
  for(int i = 0; i < nBinspT; i++) {
    // the data
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

    // the linear and constant fits
    fitLin[i]->SetLineColor(kGreen);
    fitLin[i]->SetMarkerColor(kGreen);
    fitLin[i]->Draw("same");

    fitCon[i]->SetLineColor(kRed);
    fitCon[i]->SetMarkerColor(kRed);
    fitCon[i]->Draw("same");

    // the legend
    TLegend *leg = new TLegend(0.7, 0.2, 0.9, 0.4);
    leg->SetTextSize(0.03);
    leg->AddEntry(hData[i], "data", "pl");
    leg->AddEntry(hInd[i], "indep", "pl");
    leg->AddEntry(fitLin[i], "linear", "pl");
    leg->AddEntry(fitCon[i], "constant", "pl");
    leg->Draw();

    // save and clear
    c->SaveAs(Form("plots/ratio_final/bin_%d.pdf", i));
    c->Clear();
  }
  fDep->Close();

  // draw lambda_th(pT)
  // the frame
  TH1F *func = c->DrawFrame(0, -1, 25, 1);
  func->SetXTitle("p_{T}/M");
  func->SetYTitle("#lambda_{#theta}");
  func->GetYaxis()->SetTitleOffset(1.3);
  func->GetYaxis()->SetLabelOffset(0.01);
  func->SetTitle("#lambda_{#theta} (p_{T}/M)");

  // the three lambda_th distributions
  graph_lth->SetLineColor(kBlue);
  graph_lth->SetMarkerColor(kBlue);
  graph_lth->Draw("p");
  lin_lth->SetLineColor(kGreen);
  lin_lth->SetMarkerColor(kGreen);
  lin_lth->SetFillColorAlpha(kGreen, 0.3);
  lin_lth->Draw("ce3");
  con_lth->SetLineColor(kRed);
  con_lth->SetMarkerColor(kRed);
  con_lth->SetFillColorAlpha(kRed, 0.3);
  con_lth->Draw("ce3");

  // the legend
  TLegend *leg = new TLegend(0.7, 0.2, 0.9, 0.4);
  leg->SetTextSize(0.03);
  leg->AddEntry(graph_lth, "indep", "pl");
  leg->AddEntry(lin_lth, "linear", "pl");
  leg->AddEntry(con_lth, "constant", "pl");
  leg->Draw();

  c->SaveAs("plots/ratio_final/lth.pdf");
  c->Clear();
  c->Destructor();
  
}

