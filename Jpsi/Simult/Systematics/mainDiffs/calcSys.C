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
  
  // get the fit results - baseline; mu effs; 2017 vs 2018
  TGraphErrors **graph_lth = new TGraphErrors*[3];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  fIndB->Close();
  // 2 - get lambda values for each eff model
  TFile *fIndf1 = new TFile("../../../Simult_eff1/PR_fit/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fIndf1->Get(Form("graph_lambda_J"));
  fIndf1->Close();
  TFile *fIndf2 = new TFile("../../../Simult_eff2/PR_fit/files/finalFitRes.root");
  graph_lth[2] = (TGraphErrors*)fIndf2->Get(Form("graph_lambda_J"));
  fIndf2->Close();

  // 2017 vs 2018 contribution has been fixed rather than calculated 
  double sigY = 0.01;

  // PART 2 - coarse-binned results (to expand)
  // get the histo limits
  TFile *fIn_c = new TFile("../deltaR/files/chistStore.root");
  TH2D* rHist_c;
  fIn_c->GetObject("cHistB", rHist_c);
  
  int nBinspT_c = rHist_c->GetNbinsY();
  const double *pTBins_c = rHist_c->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - rho factor
  TGraphErrors **graph_lth_c = new TGraphErrors*[3];
  // 0 - get Run2 results
  TFile *fIndB_c = new TFile("../deltaR/files/finalFitRes.root");
  graph_lth_c[0] = (TGraphErrors*)fIndB_c->Get("graph_lambda_c_B");
  graph_lth_c[1] = (TGraphErrors*)fIndB_c->Get(Form("graph_lambda_c_T"));
  graph_lth_c[2] = (TGraphErrors*)fIndB_c->Get(Form("graph_lambda_c_L"));
  fIndB_c->Close();

  // expand to all pT bins
  double diff[nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    double pt = graph_lth[0]->GetX()[i];
    for(int j = 0; j < nBinspT_c; j++) {
      double pt_min = graph_lth_c[0]->GetX()[j]-graph_lth_c[0]->GetEX()[j];
      double pt_max = graph_lth_c[0]->GetX()[j]+graph_lth_c[0]->GetEX()[j];
      if(pt > pt_min && pt < pt_max) {
	if(pt < 70) diff[i] = graph_lth_c[1]->GetY()[j]-graph_lth_c[0]->GetY()[j];
	else diff[i] = graph_lth_c[2]->GetY()[j]-graph_lth_c[0]->GetY()[j];
      }
    }
  }

  // PART 3: plotting results  
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLeftMargin(0.13);  
  
  // get the differences - positive and negative are the same
  // efficiency curves -> symmetrize, pos+neg
  TH1F *hs_sigEff = new TH1F("hs_sigEff", "hs_sigEff", nBinspT, pTBins);
  // 2017-2018 -> fixed value, symmetrize
  TH1F *hs_sigD = new TH1F("hs_sigD", "hs_sigD", nBinspT, pTBins);
  // deltaR cuts -> one per pT region, symmetrize
  TH1F *hs_sigdR = new TH1F("hs_sigdR", "hs_sigdR", nBinspT, pTBins);

  double aux, ct_a = 0;
  for(int i = 0; i < nBinspT; i++) {
    // muon eff curves
    //if(graph_lth[0]->GetX()[i] < 51)
    aux = 0.5*(abs(graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i])+abs(graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i]));
    //else aux = 0;
    if(aux < 1e-3 || ct_a == 2) {
      aux = 0;
      ct_a = 1;
    }
    hs_sigEff->SetBinContent(i+1, aux);
    // 2017 - 2018
    hs_sigD->SetBinContent(i+1, sigY);
    // dR cut
    aux = abs(diff[i]);
    if(i == 8) aux = 0;
    hs_sigdR->SetBinContent(i+1, aux);
  }

  // draw the fit results

  // FIRST - set colors, styles for elements
  hs_sigEff->SetFillColor(kGreen+1);
  hs_sigEff->SetLineColor(kGreen+1);
 
  hs_sigD->SetFillColor(kMagenta);
  hs_sigD->SetLineColor(kMagenta);
 
  hs_sigdR->SetFillColor(kOrange+1);
  hs_sigdR->SetLineColor(kOrange+1);
  
  // draw uncerts, squared and stacked
  // positive conts: eff up to 46; rho factor; constant 2017-2018
  // negative cont: stat unc
  
  // get the systs
  TH1F *f_sigEff = new TH1F("f_sigEff", "f_sigEff", nBinspT, pTBins);
  TH1F *f_sigD = new TH1F("f_sigD", "f_sigD", nBinspT, pTBins);
  TH1F *f_sigdR = new TH1F("f_sigdR", "f_sigdR", nBinspT, pTBins);
  
  for(int i = 0; i < nBinspT; i++) {
    // eff (symmetric)
    f_sigEff->SetBinContent(i+1, pow(hs_sigEff->GetBinContent(i+1), 2));
    // 2017-2018 (symmetric)
    f_sigD->SetBinContent(i+1, pow(hs_sigD->GetBinContent(i+1), 2));
    // dR cut
    f_sigdR->SetBinContent(i+1, pow(hs_sigdR->GetBinContent(i+1), 2));
  }

  // get the stats
  TH1F *f_stat = new TH1F("f_stat", "f_stat", nBinspT, pTBins);
  for(int i = 0; i < nBinspT; i++) {
    f_stat->SetBinContent(i+1, -pow(graph_lth[0]->GetEY()[i], 2));
  }
  
  // draw stacks
  double da_lim = 0.004;
  
  f_sigEff->SetFillColor(kGreen+1);
  f_sigEff->SetLineColor(kGreen+1);
  f_sigD->SetFillColor(kMagenta);
  f_sigD->SetLineColor(kMagenta);
  f_sigdR->SetFillColor(kOrange+1);
  f_sigdR->SetLineColor(kOrange+1);
  
  f_stat->SetFillColor(kBlue);
  f_stat->SetLineColor(kBlue);

  THStack *hstat = new THStack("hstat", "Title here");
  hstat->Add(f_stat);
  hstat->SetMinimum(-da_lim);
  hstat->SetMaximum(da_lim);

  THStack *hsigP = new THStack("hsigP", "Stacked #sigma^{2} (prompt)");
  hsigP->SetMinimum(-da_lim);
  hsigP->Add(f_sigD);
  hsigP->Add(f_sigdR);
  hsigP->Add(f_sigEff);
  hsigP->SetMaximum(da_lim);
  
  hsigP->Draw();
  hsigP->GetXaxis()->SetTitle("p_{T} (GeV)");
  hsigP->GetYaxis()->SetTitle("#sigma^{2}");
  hsigP->GetYaxis()->SetTitleOffset(2.);
  hstat->Draw("same");

  TLegend *legF = new TLegend(0.65, 0.7, 0.9, 0.9);
  legF->SetTextSize(0.03);
  legF->AddEntry(hs_sigEff, "Single #mu eff*", "l");
  legF->AddEntry(hs_sigdR, "#rho factor", "l");
  legF->AddEntry(hs_sigD, "2017-2018", "l");
  legF->AddEntry(f_stat, "stat", "l");
  legF->Draw();

  c->SaveAs("plots/lth_uncs_pos.pdf");
  c->Clear();  
    
  // tex table with sys uncerts per pt bin
  double sysp[nBinspT], sysn[nBinspT];
  ofstream ftex;
  ftex.open(Form("sys_unc.tex"));
  ftex << "\\begin{tabular}{c||c|c|c||c|c}\n";
  ftex << "$\\pt$ (GeV) & $\\sigma^{\\text{comb}}$ & $\\sigma^{\\text{eff}}$ & $\\sigma^\\rho$ & $\\sigma_{\\text{sys}}^{\\text{PR}}$ & $\\sigma_{\\text{stat}}^{\\text{PR}}$ \\\\\n";
  ftex << "\\hline\n";

  int p_norm = 4;
  double val_d, val_eff, val_rho;
  int rho_lims[5] = {0, 8, 9, 14, 19};
  for(int i = 0; i < nBinspT; i++) {
    // pT bin
    ftex << Form("$[%.1f, %.1f]$", pTBins[i], pTBins[i+1]);
    // syst sources (start with fixed 3 decimal places
    ftex << setprecision(p_norm) << fixed;
    // combined years: same for all pT
    if(i == 0) {
      val_d = hs_sigD->GetBinContent(i+1);
      ftex << Form(" & \\multirow{%d}{*}{$\\pm", nBinspT) << val_d << "$}";
    }
    else ftex << " & ";
    // single muon efficiency: already set to zero above
    val_eff = hs_sigEff->GetBinContent(i+1);
    if(val_eff == 0)
      ftex << " & $-$";
    else 
      ftex << " & $\\pm" << val_eff << "$";
    // rho efficiency: constant in 3 pT bins
    int ct = 0;
    for(int j = 0; j < 4; j++) {
      if(i == rho_lims[j]) {
	ct++;
	val_rho = hs_sigdR->GetBinContent(i+1);
	if(val_rho > 0)
	  ftex << Form(" & \\multirow{%d}{*}{$\\pm", rho_lims[j+1]-rho_lims[j]) << val_rho << "$}";
	else
	  ftex << Form(" & \\multirow{%d}{*}{$-$}", rho_lims[j+1]-rho_lims[j]);	  
      }
    }
    if(ct == 0) ftex << " & ";
    // full syst (sum all sources)
    double val = sqrt(pow(val_d,2)+pow(val_eff,2)+pow(val_rho,2));
    sysp[i] = val;
    sysn[i] = val;
    ftex << " & $\\pm" << val << "$";
    // statistical uncertainty
    ftex << " & $\\pm" << graph_lth[0]->GetEY()[i] << "$";
    ftex << "\\\\";
    for(int j = 0; j < 4; j++) 
      if(i+1 == rho_lims[j]) 
	ftex << "\\cline{4-4}";
    ftex << "\n";
  }
  ftex << "\\end{tabular}\n";
  ftex.close();

  // plotting the lambda_theta with total (sys+stat) error
  // but also with just stat for comparison- graph_lth[0]
  // sysp, sysn, graph_lth[0]->GetEY() (stat)
  double err_tp[nBinspT], err_tn[nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    double e_st = graph_lth[0]->GetEY()[i];
    err_tp[i] = sqrt(pow(sysp[i], 2) + pow(e_st, 2));
    err_tn[i] = sqrt(pow(sysn[i], 2) + pow(e_st, 2));
  }

  // full unc only
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("Run 2 #lambda_{#theta} (prompt J/#psi)");
  
  TGraphAsymmErrors *lth_fin = new TGraphAsymmErrors(nBinspT, graph_lth[0]->GetX(), graph_lth[0]->GetY(), graph_lth[0]->GetEX(), graph_lth[0]->GetEX(), err_tn, err_tp);
  lth_fin->SetMarkerColor(kBlack);
  lth_fin->SetLineColor(kBlack);
  lth_fin->Draw("p same");

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  c->SaveAs("plots/lth_full_unc.pdf");
  c->Clear();
    
  // full in boxes, stat normal
  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("Run 2 #lambda_{#theta} (prompt J/#psi)");

  graph_lth[0]->SetMarkerColor(kBlack);
  graph_lth[0]->SetLineColor(kBlack);
  graph_lth[0]->Draw("p same");

  lth_fin->SetFillColorAlpha(kBlack,0);
  lth_fin->Draw("5 p same");

  zero->Draw();

  c->SaveAs("plots/lth_sep_unc.pdf");
  c->Destructor();

  TFile *outfile = new TFile("lthUnc.root", "recreate");
  graph_lth[0]->SetName("graph_lambda");
  graph_lth[0]->Write();
  outfile->Close();

}
