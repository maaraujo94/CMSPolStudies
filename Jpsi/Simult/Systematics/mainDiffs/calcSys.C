// code to get the final set of systematics contributions: final plots + table
void calcSys()
{
  // PART 1 - fine-binned results (standard)
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - baseline vs mu effs
  TGraphErrors **graph_lth = new TGraphErrors*[6];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  graph_lth[3] = (TGraphErrors*)fIndB->Get("graph_lambda_NP");
  fIndB->Close();
  // 1 - get lambda values for each eff model
  TFile *fIndf1 = new TFile("../../../Simult_eff1/PR_fit/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fIndf1->Get(Form("graph_lambda_J"));
  fIndf1->Close();
  TFile *fIndf2 = new TFile("../../../Simult_eff2/PR_fit/files/finalFitRes.root");
  graph_lth[2] = (TGraphErrors*)fIndf2->Get(Form("graph_lambda_J"));
  fIndf2->Close();
  // 2 - get lambda values for the phi reweighings
  TFile *fIndP1 = new TFile("../../../Simult_phi2/PR_fit/files/finalFitRes.root");
  graph_lth[4] = (TGraphErrors*)fIndP1->Get(Form("graph_lambda_J"));
  fIndP1->Close();
  TFile *fIndP2 = new TFile("../../../Simult_phi3/PR_fit/files/finalFitRes.root");
  graph_lth[5] = (TGraphErrors*)fIndP2->Get(Form("graph_lambda_NP"));
  fIndP2->Close();

  // 2017 vs 2018 contribution has been fixed rather than calculated 
  double sigY = 0.012;

  // PART 3: plotting results  
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.03);
  
  // get the differences - positive and negative are the same
  // efficiency curves -> symmetrize, pos+neg
  TH1F *hs_sigEff = new TH1F("hs_sigEff", "hs_sigEff", nBinspT, pTBins);
  // 2017-2018 -> fixed value, symmetrize
  TH1F *hs_sigD = new TH1F("hs_sigD", "hs_sigD", nBinspT, pTBins);
  // lambda_phi -> asymmetric, different for PR and NP
  TH1F *hs_sigPhiPR = new TH1F("hs_sigPhiPR", "hs_sigPhiPR", nBinspT, pTBins);
  TH1F *hs_sigPhiNP = new TH1F("hs_sigPhiNP", "hs_sigPhiNP", nBinspT, pTBins);

  double aux, ct_a = 0;
  for(int i = 0; i < nBinspT; i++) {
    // muon eff curves
    if(graph_lth[0]->GetX()[i] < 50)
      aux = 0.5*(abs(graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i])+abs(graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i]));
    else aux = 0;
    hs_sigEff->SetBinContent(i+1, aux);
    // 2017 - 2018
    hs_sigD->SetBinContent(i+1, sigY);
    // phi reweighing
    aux = graph_lth[4]->GetY()[i]-graph_lth[0]->GetY()[i];
    if(graph_lth[0]->GetX()[i] < 37.5)
      hs_sigPhiPR->SetBinContent(i+1, aux);
    else hs_sigPhiPR->SetBinContent(i+1, 0);
    aux = graph_lth[5]->GetY()[i]-graph_lth[3]->GetY()[i];
    if(graph_lth[0]->GetX()[i] < 37.5)
      hs_sigPhiNP->SetBinContent(i+1, aux);
    else hs_sigPhiNP->SetBinContent(i+1, 0);
  }

  // draw the fit results

  // FIRST - set colors, styles for elements
  hs_sigEff->SetFillColor(kGreen+1);
  hs_sigEff->SetLineColor(kGreen+1);
 
  hs_sigD->SetFillColor(kMagenta);
  hs_sigD->SetLineColor(kMagenta);

  hs_sigPhiPR->SetFillColor(kOrange+1);
  hs_sigPhiPR->SetLineColor(kOrange+1);
  hs_sigPhiNP->SetFillColor(kOrange+1);
  hs_sigPhiNP->SetLineColor(kOrange+1);

  // draw uncerts, squared and stacked
  // positive conts: eff up to 50; constant 2017-2018
  // negative cont: stat unc
  
  // get the systs
  TH1F *f_sigEff = new TH1F("f_sigEff", "f_sigEff", nBinspT, pTBins);
  TH1F *f_sigD = new TH1F("f_sigD", "f_sigD", nBinspT, pTBins);
  TH1F *f_sigPhiPR = new TH1F("f_sigPhiPR", "f_sigPhiPR", nBinspT, pTBins);
  TH1F *f_sigPhiNP = new TH1F("f_sigPhiNP", "f_sigPhiNP", nBinspT, pTBins);
  
  for(int i = 0; i < nBinspT; i++) {
    // eff (symmetric)
    f_sigEff->SetBinContent(i+1, pow(hs_sigEff->GetBinContent(i+1), 2));
    // 2017-2018 (symmetric)
    f_sigD->SetBinContent(i+1, pow(hs_sigD->GetBinContent(i+1), 2));
    // phi weights (asymmetric)
    f_sigPhiPR->SetBinContent(i+1, pow(hs_sigPhiPR->GetBinContent(i+1), 2));
    f_sigPhiNP->SetBinContent(i+1, pow(hs_sigPhiNP->GetBinContent(i+1), 2));
  }

  // get the stats
  TH1F *f_statP = new TH1F("f_statP", "f_statP", nBinspT, pTBins);
  TH1F *f_statN = new TH1F("f_statN", "f_statN", nBinspT, pTBins);
  for(int i = 0; i < nBinspT; i++) {
    f_statP->SetBinContent(i+1, -pow(graph_lth[0]->GetEY()[i], 2));
    f_statN->SetBinContent(i+1, -pow(graph_lth[3]->GetEY()[i], 2));
  }
  
  // draw stacks
  double da_lim = 0.004;
  
  f_sigEff->SetFillColor(kGreen+1);
  f_sigEff->SetLineColor(kGreen+1);
  f_sigD->SetFillColor(kMagenta);
  f_sigD->SetLineColor(kMagenta);
  f_sigPhiPR->SetFillColor(kOrange+1);
  f_sigPhiPR->SetLineColor(kOrange+1);
  f_sigPhiNP->SetFillColor(kOrange+1);
  f_sigPhiNP->SetLineColor(kOrange+1);
  
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
  hsigP->Add(f_sigPhiPR);
  hsigP->SetMaximum(da_lim);
  
  hsigP->Draw();
  hsigP->GetXaxis()->SetTitle("p_{T} (GeV)");
  hsigP->GetYaxis()->SetTitle("#sigma^{2}");
  hsigP->GetYaxis()->SetTitleOffset(2.);
  hstatP->Draw("same");

  TLegend *legF = new TLegend(0.72, 0.7, 0.97, 0.9);
  legF->SetTextSize(0.03);
  legF->AddEntry(hs_sigEff, "Single #mu eff", "l");
  legF->AddEntry(hs_sigD, "2017-2018", "l");
  legF->AddEntry(f_statP, "stat", "l");
  legF->AddEntry(hs_sigD, "#lambda_{#phi}", "l");
  legF->Draw();

  c->SaveAs("plots/lth_uncs_stack.pdf");
  c->Clear();

  THStack *hsigN = new THStack("hsigN", "Stacked #sigma^{2} (non-prompt)");
  hsigN->SetMinimum(-da_lim);
  hsigN->Add(f_sigD);
  hsigN->Add(f_sigEff);
  hsigN->Add(f_sigPhiNP);
  hsigN->SetMaximum(da_lim);
  
  hsigN->Draw();
  hsigN->GetXaxis()->SetTitle("p_{T} (GeV)");
  hsigN->GetYaxis()->SetTitle("#sigma^{2}");
  hsigN->GetYaxis()->SetTitleOffset(2.);
  hstatN->Draw("same");

  legF->Draw();

  c->SaveAs("plots/lthNP_uncs_stack.pdf");
  c->Clear();

  // tex table with sys uncerts per pt bin
  double sysPR_P[nBinspT], sysPR_N[nBinspT];
  double sysNP_P[nBinspT], sysNP_N[nBinspT];
  ofstream ftexPR;
  ftexPR.open(Form("text_output/sys_unc.tex"));
  ftexPR << "\\begin{tabular}{c||c|c|c||c|c}\n";
  ftexPR << "$\\pt$ (GeV) & $\\sigma^{\\text{comb}}$ & $\\sigma^{\\text{eff}}$ & $\\sigma^{\\lambda_\\phi}$ & $\\sigma_{\\text{sys}}^{\\text{PR}}$ & $\\sigma_{\\text{stat}}^{\\text{PR}}$ \\\\\n";
  ftexPR << "\\hline\n";

  int p_norm = 3;
  double val_d, val_eff, val_phi;
  for(int i = 0; i < nBinspT; i++) {
    // pT bin
    ftexPR << Form("$[%.1f, %.1f]$", pTBins[i], pTBins[i+1]);
    // syst sources (start with fixed 2 decimal places)
    ftexPR << setprecision(p_norm) << fixed;
    // combined years: same for all pT
    if(i == 0) {
      val_d = hs_sigD->GetBinContent(i+1);
      ftexPR << Form(" & \\multirow{%d}{*}{$\\pm", nBinspT) << val_d << "$}";
    }
    else ftexPR << " & ";
    // single muon efficiency: already set to zero above
    val_eff = hs_sigEff->GetBinContent(i+1);
    if(val_eff == 0)
      ftexPR << " & $-$";
    else 
      ftexPR << " & $\\pm" << val_eff << "$";
    // lambda_phi: only negative
    val_phi = hs_sigPhiPR->GetBinContent(i+1);
    if(val_phi == 0)
      ftexPR << " & $-$";
    else
      ftexPR << " & $" << val_phi << "$";
    // full syst (sum all sources)
    double valP = sqrt(pow(val_d,2)+pow(val_eff,2));
    double valN = sqrt(pow(val_d,2)+pow(val_eff,2)+pow(val_phi,2));
    sysPR_P[i] = valP;
    sysPR_N[i] = valN;
    ftexPR << " & $^{+" << valP << "}_{-" << valN << "}$";
    // statistical uncertainty
    ftexPR << " & $\\pm" << graph_lth[0]->GetEY()[i] << "$";
    ftexPR << "\\\\";
    ftexPR << "\n";
  }
  ftexPR << "\\end{tabular}\n";
  ftexPR.close();

  ofstream ftexNP;
  ftexNP.open(Form("text_output/sysNP_unc.tex"));
  ftexNP << "\\begin{tabular}{c||c|c|c||c|c}\n";
  ftexNP << "$\\pt$ (GeV) & $\\sigma^{\\text{comb}}$ & $\\sigma^{\\text{eff}}$ & $\\sigma^{\\lambda_\\phi}$ & $\\sigma_{\\text{sys}}^{\\text{NP}}$ & $\\sigma_{\\text{stat}}^{\\text{NP}}$ \\\\\n";
  ftexNP << "\\hline\n";

  for(int i = 0; i < nBinspT; i++) {
    // pT bin
    ftexNP << Form("$[%.1f, %.1f]$", pTBins[i], pTBins[i+1]);
    // syst sources (start with fixed 2 decimal places)
    ftexNP << setprecision(p_norm) << fixed;
    // combined years: same for all pT
    if(i == 0) {
      val_d = hs_sigD->GetBinContent(i+1);
      ftexNP << Form(" & \\multirow{%d}{*}{$\\pm", nBinspT) << val_d << "$}";
    }
    else ftexNP << " & ";
    // single muon efficiency: already set to zero above
    val_eff = hs_sigEff->GetBinContent(i+1);
    if(val_eff == 0)
      ftexNP << " & $-$";
    else 
      ftexNP << " & $\\pm" << val_eff << "$";
    // lambda_phi: only negative
    val_phi = hs_sigPhiNP->GetBinContent(i+1);
    if(val_phi == 0)
      ftexNP << " & $-$";
    else
      ftexNP << " & $+" << val_phi << "$";
    // full syst (sum all sources)
    double valP = sqrt(pow(val_d,2)+pow(val_eff,2)+pow(val_phi,2));
    double valN = sqrt(pow(val_d,2)+pow(val_eff,2));
    sysNP_P[i] = valP;
    sysNP_N[i] = valN;
    ftexNP << " & $^{+" << valP << "}_{-" << valN << "}$";
    // statistical uncertainty
    ftexNP << " & $\\pm" << graph_lth[0]->GetEY()[i] << "$";
    ftexNP << "\\\\";
    ftexNP << "\n";
  }
  ftexNP << "\\end{tabular}\n";
  ftexNP.close();

  // plotting the lambda_theta with total (sys+stat) error
  // but also with just stat for comparison- graph_lth[0]
  // sysp, sysn, graph_lth[0]->GetEY() (stat)
  double errPR_P[nBinspT], errPR_N[nBinspT];
  double errNP_P[nBinspT], errNP_N[nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    double e_st = graph_lth[0]->GetEY()[i];
    errPR_P[i] = sqrt(pow(sysPR_P[i],2) + pow(e_st,2));
    errPR_N[i] = sqrt(pow(sysPR_N[i],2) + pow(e_st,2));
    e_st = graph_lth[3]->GetEY()[i];
    errNP_P[i] = sqrt(pow(sysNP_P[i],2) + pow(e_st,2));
    errNP_N[i] = sqrt(pow(sysNP_N[i],2) + pow(e_st,2));
  }

  // full unc only
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("Run 2 #lambda_{#theta}");
  
  TGraphAsymmErrors *lth_fPR = new TGraphAsymmErrors(nBinspT, graph_lth[0]->GetX(), graph_lth[0]->GetY(), graph_lth[0]->GetEX(), graph_lth[0]->GetEX(), errPR_P, errPR_N);
  lth_fPR->SetMarkerColor(kBlue);
  lth_fPR->SetLineColor(kBlue);
  lth_fPR->Draw("p same");

  TGraphAsymmErrors *lth_fNP = new TGraphAsymmErrors(nBinspT, graph_lth[3]->GetX(), graph_lth[3]->GetY(), graph_lth[3]->GetEX(), graph_lth[3]->GetEX(), errNP_P, errNP_N);
  lth_fNP->SetMarkerColor(kRed);
  lth_fNP->SetLineColor(kRed);
  lth_fNP->Draw("p same");

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg2 = new TLegend(0.67, 0.7, 0.97, 0.9);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(lth_fPR, "prompt J/#psi", "pl");
  leg2->AddEntry(lth_fNP, "non-prompt J/#psi", "pl");
  leg2->Draw();

  c->SaveAs("plots/lth_full_unc.pdf");
  c->Clear();
  c->Destructor();    
}
