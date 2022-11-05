void texTemp()
{
  double sig_eff[] = {0.0108, 0.0085, 0.0066, 0.0058, 0.0047, 0.0044, 0.0040, 0.0033, 0.0022, 0.0012};
  
  // PART 1 - fine-binned results (standard)
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();

  // get the fit results - baseline PR and NP
  TGraphErrors **graph_lth = new TGraphErrors*[2];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  graph_lth[1] = (TGraphErrors*)fIndB->Get("graph_lambda_NP");
  fIndB->Close();

  // tex table with sys uncerts per pt bin
  double sysp[nBinspT], sysn[nBinspT];
  ofstream ftex;
  ftex.open(Form("sys_unc.tex"));
  ftex << "\\begin{tabular}{c||c|c||c|c|c}\n";
  ftex << "$\\pt$ (GeV) & $\\sigma^{\\text{comb}}$ & $\\sigma^{\\text{eff}}$  & $\\sigma_{\\text{sys}}$ & $\\sigma_{\\text{stat}}^{\\text{PR}}$  & $\\sigma_{\\text{stat}}^{\\text{NP}}$ \\\\\n";
  ftex << "\\hline\n";

  int p_norm = 4;
  double val_d, val_eff;
  for(int i = 0; i < nBinspT; i++) {
    // pT bin
    ftex << Form("$[%.1f, %.1f]$", pTBins[i], pTBins[i+1]);
    // syst sources (start with fixed 3 decimal places
    ftex << setprecision(p_norm) << fixed;
    // combined years: same for all pT
    if(i == 0) {
      val_d = 0.012;
      ftex << Form(" & \\multirow{%d}{*}{$\\pm", nBinspT) << val_d << "$}";
    }
    else ftex << " & ";
    // single muon efficiency: already set to zero above
    if(i < 10) {
      val_eff = sig_eff[i];
      ftex << " & $\\pm" << val_eff << "$";
    }
    else {
      val_eff = 0;
      ftex << " & $-$";
    }
    // full syst (sum all sources)
    double val = sqrt(pow(val_d,2)+pow(val_eff,2));
    sysp[i] = val;
    sysn[i] = val;
    ftex << " & $\\pm" << val << "$";
    // statistical uncertainty
    ftex << " & $\\pm" << graph_lth[0]->GetEY()[i] << "$";
    ftex << " & $\\pm" << graph_lth[1]->GetEY()[i] << "$";
    ftex << "\\\\";
    ftex << "\n";
  }
  ftex << "\\end{tabular}\n";
  ftex.close();
  
  // draw uncerts, squared and stacked
  // positive conts: eff up to 46; rho factor; constant 2017-2018
  // negative cont: stat unc
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.03);
  
  // get the systs
  TH1F *f_sigEff = new TH1F("f_sigEff", "f_sigEff", nBinspT, pTBins);
  TH1F *f_sigD = new TH1F("f_sigD", "f_sigD", nBinspT, pTBins);
  
  for(int i = 0; i < nBinspT; i++) {
    // eff (symmetric)
    if(i < 10) 
      f_sigEff->SetBinContent(i+1, pow(sig_eff[i], 2));
    else 
      f_sigEff->SetBinContent(i+1, 0);
      

    // 2017-2018 (symmetric)
    f_sigD->SetBinContent(i+1, pow(0.012, 2));
  }

  // get the stats
  TH1F *f_statP = new TH1F("f_statP", "f_statP", nBinspT, pTBins);
  TH1F *f_statN = new TH1F("f_statN", "f_statN", nBinspT, pTBins);
  for(int i = 0; i < nBinspT; i++) {
    f_statP->SetBinContent(i+1, -pow(graph_lth[0]->GetEY()[i], 2));
    f_statN->SetBinContent(i+1, -pow(graph_lth[1]->GetEY()[i], 2));
  }
  
  // draw stacks
  double da_lim = 0.004;
  
  f_sigEff->SetFillColor(kGreen+1);
  f_sigEff->SetLineColor(kGreen+1);
  f_sigD->SetFillColor(kMagenta);
  f_sigD->SetLineColor(kMagenta);
  
  f_statP->SetFillColor(kBlue);
  f_statP->SetLineColor(kBlue);
  f_statN->SetFillColor(kBlue);
  f_statN->SetLineColor(kBlue);

  THStack *hstatP = new THStack("hstatP", "Title here");
  hstatP->Add(f_statP);
  hstatP->SetMinimum(-da_lim);
  hstatP->SetMaximum(da_lim);
  THStack *hstatN = new THStack("hstatN", "Title here");
  hstatN->Add(f_statN);
  hstatN->SetMinimum(-da_lim);
  hstatN->SetMaximum(da_lim);

  THStack *hsigP = new THStack("hsigP", "Stacked #sigma^{2} (prompt)");
  hsigP->SetMinimum(-da_lim);
  hsigP->Add(f_sigD);
  hsigP->Add(f_sigEff);
  hsigP->SetMaximum(da_lim);
  
  hsigP->Draw();
  hsigP->GetXaxis()->SetTitle("p_{T} (GeV)");
  hsigP->GetYaxis()->SetTitle("#sigma^{2}");
  hsigP->GetYaxis()->SetTitleOffset(2.);
  hstatP->Draw("same");

  TLegend *legF = new TLegend(0.72, 0.7, 0.97, 0.9);
  legF->SetTextSize(0.03);
  legF->AddEntry(f_sigEff, "Single #mu eff", "l");
  legF->AddEntry(f_sigD, "2017-2018", "l");
  legF->AddEntry(f_statP, "stat", "l");
  legF->Draw();

  c->SaveAs("plots/lth_uncs_pos.pdf");
  c->Clear();

  hsigP->SetTitle("Stacked #sigma^{2} (non-prompt)");
  hsigP->Draw();
  hsigP->GetXaxis()->SetTitle("p_{T} (GeV)");
  hsigP->GetYaxis()->SetTitle("#sigma^{2}");
  hsigP->GetYaxis()->SetTitleOffset(2.);
  hstatN->Draw("same");
  
  legF->Draw();
  c->SaveAs("plots/lthNP_uncs_pos.pdf");
  c->Clear();

  // plotting the lambda_theta with total (sys+stat) error
  double err_tP[nBinspT], err_tN[nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    double e_st = graph_lth[0]->GetEY()[i];
    err_tP[i] = sqrt(pow(sysp[i], 2) + pow(e_st, 2));
    e_st = graph_lth[1]->GetEY()[i];
    err_tN[i] = sqrt(pow(sysn[i], 2) + pow(e_st, 2));
  }

  // full unc only
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("Run 2 #lambda_{#theta}");
  
  TGraphErrors *lth_fP = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), graph_lth[0]->GetY(), graph_lth[0]->GetEX(), err_tP);
  lth_fP->SetMarkerColor(kBlue);
  lth_fP->SetLineColor(kBlue);
  lth_fP->Draw("p same");

  TGraphErrors *lth_fN = new TGraphErrors(nBinspT, graph_lth[1]->GetX(), graph_lth[1]->GetY(), graph_lth[1]->GetEX(), err_tN);
  lth_fN->SetMarkerColor(kRed);
  lth_fN->SetLineColor(kRed);
  lth_fN->Draw("p same");

  cout << lth_fP->GetEY()[0] << " " << graph_lth[0]->GetEY()[0] << " " << sqrt(pow(lth_fP->GetEY()[0],2) - pow(graph_lth[0]->GetEY()[0],2)) << endl; 
  cout << lth_fN->GetEY()[0] << " " << graph_lth[1]->GetEY()[0] << " " << sqrt(pow(lth_fN->GetEY()[0],2) - pow(graph_lth[1]->GetEY()[0],2)) << endl;
 
  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg2 = new TLegend(0.67, 0.7, 0.97, 0.9);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(lth_fP, "prompt J/#psi", "pl");
  leg2->AddEntry(lth_fN, "non-prompt J/#psi", "pl");
  leg2->Draw();

  c->SaveAs("plots/lth_full_unc.pdf");
  c->Clear();

  c->Destructor();
}
