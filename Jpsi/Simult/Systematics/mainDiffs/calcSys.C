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
  TGraphErrors **graph_lth = new TGraphErrors*[8];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  graph_lth[3] = (TGraphErrors*)fIndB->Get("graph_lambda_NP");
  fIndB->Close();
  // 1 - get lambda values for each eff model
  TFile *fIndf1 = new TFile("../../../Simult_eff1/PR_fit/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fIndf1->Get(Form("graph_lambda_J"));
  graph_lth[6] = (TGraphErrors*)fIndf1->Get(Form("graph_lambda_NP"));
  fIndf1->Close();
  TFile *fIndf2 = new TFile("../../../Simult_eff2/PR_fit/files/finalFitRes.root");
  graph_lth[2] = (TGraphErrors*)fIndf2->Get(Form("graph_lambda_J"));
  graph_lth[7] = (TGraphErrors*)fIndf2->Get(Form("graph_lambda_NP"));
  fIndf2->Close();
  // 2 - get lambda values for the phi reweighings
  TFile *fIndP1 = new TFile("../../../Simult_phi2/PR_fit/files/finalFitRes.root");
  graph_lth[4] = (TGraphErrors*)fIndP1->Get(Form("graph_lambda_J"));
  fIndP1->Close();
  TFile *fIndP2 = new TFile("../../../Simult_phi3/PR_fit/files/finalFitRes.root");
  graph_lth[5] = (TGraphErrors*)fIndP2->Get(Form("graph_lambda_NP"));
  fIndP2->Close();

  // PART 3: plotting results  
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);
  
  // get the differences - positive and negative are the same
  // efficiency curves -> symmetrize, pos+neg, diff for PR and NP
  TH1F *hs_sigEffPR = new TH1F("hs_sigEffPR", "hs_sigEffPR", nBinspT, pTBins);
  TH1F *hs_sigEffNP = new TH1F("hs_sigEffNP", "hs_sigEffNP", nBinspT, pTBins);
  // lambda_phi -> asymmetric, different for PR and NP
  TH1F *hs_sigPhiPR = new TH1F("hs_sigPhiPR", "hs_sigPhiPR", nBinspT, pTBins);
  TH1F *hs_sigPhiNP = new TH1F("hs_sigPhiNP", "hs_sigPhiNP", nBinspT, pTBins);
  
  double aux, ct_a = 0;
  for(int i = 0; i < nBinspT; i++) {
    // muon eff curves
    if(graph_lth[0]->GetX()[i] < 50)
      aux = 0.5*(abs(graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i])+abs(graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i]));
    else aux = 0;
    hs_sigEffPR->SetBinContent(i+1, aux);
    if(graph_lth[0]->GetX()[i] < 50)
      aux = 0.5*(abs(graph_lth[6]->GetY()[i] - graph_lth[3]->GetY()[i])+abs(graph_lth[7]->GetY()[i] - graph_lth[3]->GetY()[i]));
    else aux = 0;
    hs_sigEffNP->SetBinContent(i+1, aux);
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
  hs_sigEffPR->SetFillColor(kGreen+1);
  hs_sigEffPR->SetLineColor(kGreen+1);
  hs_sigEffNP->SetFillColor(kGreen+1);
  hs_sigEffNP->SetLineColor(kGreen+1);

  hs_sigPhiPR->SetFillColor(kBlue);
  hs_sigPhiPR->SetLineColor(kBlue);
  hs_sigPhiNP->SetFillColor(kBlue);
  hs_sigPhiNP->SetLineColor(kBlue);


  // draw uncerts, squared and stacked
  // positive conts: eff up to 50
  // negative cont: stat unc
  
  // get the systs
  TH1F *f_sigEffPR = new TH1F("f_sigEffPR", "f_sigEffPR", nBinspT, pTBins);
  TH1F *f_sigEffNP = new TH1F("f_sigEffNP", "f_sigEffNP", nBinspT, pTBins);
  TH1F *f_sigPhiPR = new TH1F("f_sigPhiPR", "f_sigPhiPR", nBinspT, pTBins);
  TH1F *f_sigPhiNP = new TH1F("f_sigPhiNP", "f_sigPhiNP", nBinspT, pTBins);
  
  for(int i = 0; i < nBinspT; i++) {
    // eff (symmetric)
    f_sigEffPR->SetBinContent(i+1, pow(hs_sigEffPR->GetBinContent(i+1), 2));
    f_sigEffNP->SetBinContent(i+1, pow(hs_sigEffNP->GetBinContent(i+1), 2));
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
  
  f_sigEffPR->SetFillColor(kGreen+1);
  f_sigEffPR->SetLineColor(kGreen+1);
  f_sigEffNP->SetFillColor(kGreen+1);
  f_sigEffNP->SetLineColor(kGreen+1);
  f_sigPhiPR->SetFillColor(kBlue);
  f_sigPhiPR->SetLineColor(kBlue);
  f_sigPhiNP->SetFillColor(kBlue);
  f_sigPhiNP->SetLineColor(kBlue);
  
  f_statP->SetFillColor(kGray);
  f_statP->SetLineColor(kGray);
  f_statN->SetFillColor(kGray);
  f_statN->SetLineColor(kGray);

  THStack *hstatP = new THStack("hstatP", "");
  hstatP->Add(f_statP);
  hstatP->SetMinimum(-da_lim);
  hstatP->SetMaximum(da_lim);
  THStack *hstatN = new THStack("hstatN", "");
  hstatN->Add(f_statN);
  hstatN->SetMinimum(-da_lim);
  hstatN->SetMaximum(da_lim);

  THStack *hsigP = new THStack("hsigP", "");
  hsigP->SetMinimum(-da_lim);
  hsigP->Add(f_sigEffPR);
  hsigP->Add(f_sigPhiPR);
  hsigP->SetMaximum(da_lim);
  
  hsigP->Draw();
  hsigP->GetXaxis()->SetTitle("p_{T} (GeV)");
  hsigP->GetYaxis()->SetTitle("#sigma^{2}");
  hsigP->GetYaxis()->SetTitleOffset(2.);
  hstatP->Draw("same");

  TLegend *legF = new TLegend(0.7, 0.785, 0.97, 0.985);
  legF->SetTextSize(0.03);
  legF->AddEntry(hs_sigPhiPR, "#beta (only negative)", "l");
  legF->AddEntry(hs_sigEffPR, "Single #mu eff", "l");
  legF->AddEntry(f_statP, "stat", "l");
  legF->Draw();

  c->SaveAs("plots/lth_uncs_stack.pdf");
  c->Clear();

  THStack *hsigN = new THStack("hsigN", "");
  hsigN->SetMinimum(-da_lim);
  hsigN->Add(f_sigEffNP);
  hsigN->Add(f_sigPhiNP);
  hsigN->SetMaximum(da_lim);
  
  hsigN->Draw();
  hsigN->GetXaxis()->SetTitle("p_{T} (GeV)");
  hsigN->GetYaxis()->SetTitle("#sigma^{2}");
  hsigN->GetYaxis()->SetTitleOffset(2.);
  hstatN->Draw("same");

  TLegend *legN = new TLegend(0.7, 0.785, 0.97, 0.985);
  legN->SetTextSize(0.03);
  legN->AddEntry(hs_sigPhiNP, "#beta (only positive)", "l");
  legN->AddEntry(hs_sigEffNP, "Single #mu eff", "l");
  legN->AddEntry(f_statP, "stat", "l");
  legN->Draw();

  c->SaveAs("plots/lthNP_uncs_stack.pdf");
  c->Clear();

  // tex table with sys uncerts per pt bin
  double sysPR_P[nBinspT], sysPR_N[nBinspT];
  double sysNP_P[nBinspT], sysNP_N[nBinspT];
  ofstream ftexPR;
  ftexPR.open(Form("text_output/sysPR_unc.tex"));
  ftexPR << "\\begin{tabular}{c|cc|c}\n";
  ftexPR << "$\\pt$ (GeV)  & $\\mu$ eff. & $\\beta$ & total \\\\\n";
  ftexPR << "\\hline\n";

  int p_norm = 3;
  double val_eff, val_phi;
  for(int i = 0; i < nBinspT; i++) {
    // pT bin
    if(pTBins[i] - (int)pTBins[i] == 0)
      ftexPR << Form("%.0f--", pTBins[i]);
    else
      ftexPR << Form("%.1f--", pTBins[i]);
    if(pTBins[i+1] - (int)pTBins[i+1] == 0)
      ftexPR << Form("%.0f", pTBins[i+1]);
    else
      ftexPR << Form("%.1f", pTBins[i+1]);
    // syst sources (start with fixed 2 decimal places)
    ftexPR << setprecision(p_norm) << fixed;
    // single muon efficiency: already set to zero above
    val_eff = hs_sigEffPR->GetBinContent(i+1);
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
    // full syst (sum all sources) - use max
    double valP = sqrt(pow(val_eff,2));
    double valN = sqrt(pow(val_eff,2)+pow(val_phi,2));
    sysPR_P[i] = valP;
    sysPR_N[i] = valN;
    double valAv = max(valP, valN);
    if(valAv > 0)
      ftexPR << " & $\\pm" << valAv << "$";
    else
      ftexPR << " & $-$";      
    ftexPR << "\\\\";
    ftexPR << "\n";
  }
  ftexPR << "\\end{tabular}\n";
  ftexPR.close();

  ofstream ftexNP;
  ftexNP.open(Form("text_output/sysNP_unc.tex"));
  ftexNP << "\\begin{tabular}{c|cc|c}\n";
  ftexNP << "$\\pt$ (GeV)  & $\\mu$ eff. & $\\beta$ & total \\\\\n";
  ftexNP << "\\hline\n";

  for(int i = 0; i < nBinspT; i++) {
    // pT bin
    if(pTBins[i] - (int)pTBins[i] == 0)
      ftexNP << Form("%.0f--", pTBins[i]);
    else
      ftexNP << Form("%.1f--", pTBins[i]);
    if(pTBins[i+1] - (int)pTBins[i+1] == 0)
      ftexNP << Form("%.0f", pTBins[i+1]);
    else
      ftexNP << Form("%.1f", pTBins[i+1]);
    // syst sources (start with fixed 2 decimal places)
    ftexNP << setprecision(p_norm) << fixed;
    // single muon efficiency: already set to zero above
    val_eff = hs_sigEffNP->GetBinContent(i+1);
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
    double valP = sqrt(pow(val_eff,2)+pow(val_phi,2));
    double valN = sqrt(pow(val_eff,2));
    sysNP_P[i] = valP;
    sysNP_N[i] = valN;
    double valAv = max(valP, valN);
    if(valAv > 0)
      ftexNP << " & $\\pm" << valAv << "$";
    else
      ftexNP << " & $-$";      
    ftexNP << "\\\\";
    ftexNP << "\n";
  }
  ftexNP << "\\end{tabular}\n";
  ftexNP.close();

  // tex tables with central values + uncerts
  ofstream ftexLPR;
  ftexLPR.open(Form("text_output/lthPR.tex"));
  ftexLPR << "\\begin{tabular}{c|cccc}\n";
  ftexLPR << "$\\pt$ (GeV) & $\\lambda_\\theta$ & $\\sigma_{\\text{stat}}$ & $\\sigma_{\\text{sys}}$ & $\\sigma_{\\text{tot}}$  \\\\\n";
  ftexLPR << "\\hline\n";
  for(int i = 0; i < nBinspT; i++) {
    // pT bin
    if(pTBins[i] - (int)pTBins[i] == 0)
      ftexLPR << Form("%.0f--", pTBins[i]);
    else
      ftexLPR << Form("%.1f--", pTBins[i]);
    if(pTBins[i+1] - (int)pTBins[i+1] == 0)
      ftexLPR << Form("%.0f", pTBins[i+1]);
    else
      ftexLPR << Form("%.1f", pTBins[i+1]);
    ftexLPR << setprecision(p_norm) << fixed;
    // central value
    ftexLPR << "& $" << graph_lth[0]->GetY()[i] << "$";
    // statistical uncertainty
    ftexLPR << " & $\\pm" << graph_lth[0]->GetEY()[i] << "$";
    // systematic uncertainty
    if(max(sysPR_P[i],sysPR_N[i]) > 0)
      ftexLPR << " & $\\pm" << max(sysPR_P[i],sysPR_N[i]) << "$";
    else
      ftexLPR << " & $-$";
    // total unc
    double valAv = max(sysPR_P[i],sysPR_N[i]);
    ftexLPR << " & $\\pm" << sqrt(pow(graph_lth[0]->GetEY()[i], 2) + pow(valAv, 2)) << "$";
    ftexLPR << "\\\\";
    ftexLPR << "\n";
  }
  ftexLPR << "\\end{tabular}\n";
  ftexLPR.close();

  // tex tables with central values + uncerts
  ofstream ftexLNP;
  ftexLNP.open(Form("text_output/lthNP.tex"));
  ftexLNP << "\\begin{tabular}{c|cccc}\n";
  ftexLNP << "$\\pt$ (GeV) & $\\lambda_\\theta$ & $\\sigma_{\\text{stat}}$ & $\\sigma_{\\text{sys}}$ & $\\sigma_{\\text{tot}}$  \\\\\n";
  ftexLNP << "\\hline\n";
  for(int i = 0; i < nBinspT; i++) {
    // pT bin
    if(pTBins[i] - (int)pTBins[i] == 0)
      ftexLNP << Form("%.0f--", pTBins[i]);
    else
      ftexLNP << Form("%.1f--", pTBins[i]);
    if(pTBins[i+1] - (int)pTBins[i+1] == 0)
      ftexLNP << Form("%.0f", pTBins[i+1]);
    else
      ftexLNP << Form("%.1f", pTBins[i+1]);
    ftexLNP << setprecision(p_norm) << fixed;
    // central value
    ftexLNP << "& $" << graph_lth[3]->GetY()[i] << "$";
    // statistical uncertainty
    ftexLNP << " & $\\pm" << graph_lth[3]->GetEY()[i] << "$";
    // systematic uncertainty
    if(max(sysNP_P[i],sysNP_N[i]) > 0)
      ftexLNP << " & $\\pm" << max(sysNP_P[i],sysNP_N[i]) << "$";
    else
      ftexLNP << " & $-$";
    // total unc
    double valAv = max(sysNP_P[i],sysNP_N[i]);
    ftexLNP << " & $\\pm" << sqrt(pow(graph_lth[3]->GetEY()[i], 2) + pow(valAv, 2)) << "$";
    ftexLNP << "\\\\";
    ftexLNP << "\n";
  }
  ftexLNP << "\\end{tabular}\n";
  ftexLNP.close();

  // extra table - all syst, total sys (pos vs neg), stat, total sum (pos vs neg)
  ofstream fexPR;
  fexPR.open(Form("text_output/uncPR_comp.tex"));
  fexPR << "\\begin{tabular}{c|cc|cc|c|cc}\n";
  fexPR << "$\\pt$ (GeV)  & $\\mu$ eff. & $\\beta$ & total sys (pos) & total sys (neg) & $\\sigma_{\\text{stat}}$ & total unc (pos) & total unc (neg) \\\\\n";
  fexPR << "\\hline\n";

  p_norm = 3;
  for(int i = 0; i < nBinspT; i++) {
    // pT bin
    if(pTBins[i] - (int)pTBins[i] == 0)
      fexPR << Form("%.0f--", pTBins[i]);
    else
      fexPR << Form("%.1f--", pTBins[i]);
    if(pTBins[i+1] - (int)pTBins[i+1] == 0)
      fexPR << Form("%.0f", pTBins[i+1]);
    else
      fexPR << Form("%.1f", pTBins[i+1]);
    // syst sources (start with fixed 2 decimal places)
    fexPR << setprecision(p_norm) << fixed;
    // single muon efficiency: already set to zero above
    val_eff = hs_sigEffPR->GetBinContent(i+1);
    if(val_eff == 0)
      fexPR << " & $-$";
    else 
      fexPR << " & $\\pm" << val_eff << "$";
    // lambda_phi: only negative
    val_phi = hs_sigPhiPR->GetBinContent(i+1);
    if(val_phi == 0)
      fexPR << " & $-$";
    else
      fexPR << " & $" << val_phi << "$";
    // full syst (sum all sources) - pos, then neg
    if(sysPR_P[i] > 0)
      fexPR << " & $+" << sysPR_P[i] << "$";
    else
      fexPR << " & $-$";      
    if(sysPR_N[i] > 0)
      fexPR << " & $-" << sysPR_N[i] << "$";
    else
      fexPR << " & $-$";
    // statistical uncertainty
    fexPR << " & $\\pm" << graph_lth[0]->GetEY()[i] << "$";
    // total unc - pos, then neg
    fexPR << " & $+" << sqrt(pow(graph_lth[0]->GetEY()[i], 2) + pow(sysPR_P[i], 2)) << "$";
    fexPR << " & $-" << sqrt(pow(graph_lth[0]->GetEY()[i], 2) + pow(sysPR_N[i], 2)) << "$";
    fexPR << "\\\\";
    fexPR << "\n";
  }
  fexPR << "\\end{tabular}\n";
  fexPR.close();

  ofstream fexNP;
  fexNP.open(Form("text_output/uncNP_comp.tex"));
  fexNP << "\\begin{tabular}{c|cc|cc|c|cc}\n";
  fexNP << "$\\pt$ (GeV)  & $\\mu$ eff. & $\\beta$ & total sys (pos) & total sys (neg) & $\\sigma_{\\text{stat}}$ & total unc (pos) & total unc (neg) \\\\\n";
  fexNP << "\\hline\n";

  p_norm = 3;
  for(int i = 0; i < nBinspT; i++) {
    // pT bin
    if(pTBins[i] - (int)pTBins[i] == 0)
      fexNP << Form("%.0f--", pTBins[i]);
    else
      fexNP << Form("%.1f--", pTBins[i]);
    if(pTBins[i+1] - (int)pTBins[i+1] == 0)
      fexNP << Form("%.0f", pTBins[i+1]);
    else
      fexNP << Form("%.1f", pTBins[i+1]);
    // syst sources (start with fixed 2 decimal places)
    fexNP << setprecision(p_norm) << fixed;
    // single muon efficiency: already set to zero above
    val_eff = hs_sigEffNP->GetBinContent(i+1);
    if(val_eff == 0)
      fexNP << " & $-$";
    else 
      fexNP << " & $\\pm" << val_eff << "$";
    // lambda_phi: only positive
    val_phi = hs_sigPhiNP->GetBinContent(i+1);
    if(val_phi == 0)
      fexNP << " & $-$";
    else
      fexNP << " & $+" << val_phi << "$";
    // full syst (sum all sources) - pos, then neg
    if(sysNP_P[i] > 0)
      fexNP << " & $+" << sysNP_P[i] << "$";
    else
      fexNP << " & $-$";      
    if(sysNP_N[i] > 0)
      fexNP << " & $-" << sysNP_N[i] << "$";
    else
      fexNP << " & $-$";
    // statistical uncertainty
    fexNP << " & $\\pm" << graph_lth[3]->GetEY()[i] << "$";
    // total unc - pos, then neg
    fexNP << " & $+" << sqrt(pow(graph_lth[3]->GetEY()[i], 2) + pow(sysNP_P[i], 2)) << "$";
    fexNP << " & $-" << sqrt(pow(graph_lth[3]->GetEY()[i], 2) + pow(sysNP_N[i], 2)) << "$";
    fexNP << "\\\\";
    fexNP << "\n";
  }
  fexNP << "\\end{tabular}\n";
  fexNP.close();

  
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
  //fl1->SetTitle("Run 2 #lambda_{#theta}");
  
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

  TLegend *leg2 = new TLegend(0.67, 0.785, 0.97, 0.985);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(lth_fPR, "prompt J/#psi", "pl");
  leg2->AddEntry(lth_fNP, "non-prompt J/#psi", "pl");
  leg2->Draw();

  c->SaveAs("plots/lth_full_unc.pdf");
  c->Clear();
  c->Destructor();

  TFile *fout = new TFile("files/finalUnc.root", "recreate");
  lth_fPR->Write("lth_fPR");
  lth_fNP->Write("lth_fNP");
  fout->Close();
}
