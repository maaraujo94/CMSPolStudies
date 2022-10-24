// code to plot the fit results

void combIndSim_NP()
{
  // get the histo limits
  TFile *fIn = new TFile("../2017/PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get lambda values for each bin
  const int n_lbl = 1;
  string lbl[] = {"NP"};
  TFile *fInd7 = new TFile("../2017/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth7 = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lth7[i_t] = (TGraphErrors*)fInd7->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd7->Close();
  TFile *fInd8 = new TFile("../2018/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth8 = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lth8[i_t] = (TGraphErrors*)fInd8->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd8->Close();
  TFile *fIndS = new TFile("../Simult/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lthS = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lthS[i_t] = (TGraphErrors*)fIndS->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fIndS->Close();

  // get the systematic uncertainty envelope
  int nS = graph_lthS[0]->GetN();
  double unc_y[nS];
  for(int i = 0; i < nS; i++) {
    unc_y[i] = 0.0111;
  }
  TGraphErrors *g_unc = new TGraphErrors(nS, graph_lthS[0]->GetX(), graph_lthS[0]->GetY(), graph_lthS[0]->GetEX(), unc_y);

  // combine both sets of data for fitting
  int nC = graph_lth7[0]->GetN() * 2;
  double xvC[nC], xeC[nC], yvC[nC], yeC[nC];
  for(int i = 0; i < nC/2; i++) {
    xvC[i] = graph_lth7[0]->GetX()[i];
    yvC[i] = graph_lth7[0]->GetY()[i];
    xeC[i] = graph_lth7[0]->GetEX()[i];
    yeC[i] = graph_lth7[0]->GetEY()[i];
    
    xvC[i+nC/2] = graph_lth8[0]->GetX()[i];
    yvC[i+nC/2] = graph_lth8[0]->GetY()[i];
    xeC[i+nC/2] = graph_lth8[0]->GetEX()[i];
    yeC[i+nC/2] = graph_lth8[0]->GetEY()[i];    
  }
  TGraphErrors *g_lth = new TGraphErrors(nC, xvC, yvC, xeC, yeC);
  TGraphErrors *g_lthSF = new TGraphErrors(graph_lthS[0]->GetN(), graph_lthS[0]->GetX(), graph_lthS[0]->GetY(), graph_lthS[0]->GetEX(), graph_lthS[0]->GetEY());
  
  TF1 *flin = new TF1("flin", "[0]*x+[1]", pTBins[0], pTBins[nBinspT]);
  flin->SetParameters(0.1, 0.2);
  flin->SetLineColor(kBlack);
  flin->SetLineStyle(kDashed);
  g_lth->Fit(flin);

  TF1 *flinS = new TF1("flinS", "[0]*x+[1]", pTBins[0], pTBins[nBinspT]);
  flinS->SetParameters(0.1, 0.2);
  flinS->SetLineColor(kViolet+2);
  flinS->SetLineStyle(kDashed);
  g_lthSF->Fit(flinS);
   
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("non-prompt J/#psi #lambda_{#theta}");

  graph_lth7[0]->SetLineColor(kBlue);
  graph_lth7[0]->SetMarkerStyle(20);
  graph_lth7[0]->SetMarkerSize(.5);
  graph_lth7[0]->SetMarkerColor(kBlue);
  graph_lth7[0]->Draw("px same");
  graph_lth8[0]->SetLineColor(kRed);
  graph_lth8[0]->SetMarkerStyle(20);
  graph_lth8[0]->SetMarkerSize(.5);
  graph_lth8[0]->SetMarkerColor(kRed);
  graph_lth8[0]->Draw("px same");
  graph_lthS[0]->SetLineColor(kBlack);
  graph_lthS[0]->SetMarkerStyle(20);
  graph_lthS[0]->SetMarkerSize(.5);
  graph_lthS[0]->SetMarkerColor(kBlack);
  graph_lthS[0]->Draw("p same");
 
  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg = new TLegend(0.75, 0.75, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(graph_lth7[0], "2017", "pl");
  leg->AddEntry(graph_lth8[0], "2018", "pl");
  leg->AddEntry(graph_lthS[0], "Run2", "pl");
  leg->Draw();
  
  c->SaveAs("par_lthNP_all.pdf");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.5);
  g_unc->Draw("3");

  c->SaveAs("par_lthNP_all_band.pdf");
  
  flin->Draw("lsame");
  flinS->Draw("lsame");

  TLatex lc;
  lc.SetTextSize(0.03);
  lc.DrawLatex(30, 0.85, Form("m = ( %.2f #pm %.2f ) #times 10^{-3}", flin->GetParameter(0)*1e3, flin->GetParError(0)*1e3));
  lc.DrawLatex(30, 0.7, Form("b = %.3f #pm %.3f", flin->GetParameter(1), flin->GetParError(1)));
  lc.DrawLatex(30, 0.55, Form("#chi^{2}/ndf = %.0f/%d", flin->GetChisquare(), flin->GetNDF()));
  lc.DrawLatex(60, 0.55, Form("P(#chi^{2},ndf) = %.1f %%", TMath::Prob(flin->GetChisquare(), flin->GetNDF())*100.));

  lc.SetTextColor(kViolet+2); 
  lc.DrawLatex(30, -0.15, Form("m = ( %.2f #pm %.2f ) #times 10^{-3}", flinS->GetParameter(0)*1e3, flinS->GetParError(0)*1e3));
  lc.DrawLatex(30, -0.3, Form("b = %.3f #pm %.3f", flinS->GetParameter(1), flinS->GetParError(1)));
  lc.DrawLatex(30, -0.45, Form("#chi^{2}/ndf = %.0f/%d", flinS->GetChisquare(), flinS->GetNDF()));
  lc.DrawLatex(60, -0.45, Form("P(#chi^{2},ndf) = %.1f %%", TMath::Prob(flinS->GetChisquare(), flinS->GetNDF())*100.));

    cout << flin->GetChisquare() << "/" << flin->GetNDF() << " -> " << TMath::Prob(flin->GetChisquare(), flin->GetNDF()) << endl;

  c->SaveAs("fit_lthNP_all.pdf");
  c->Clear();
  c->Destructor();

  fIn->Close();
}
