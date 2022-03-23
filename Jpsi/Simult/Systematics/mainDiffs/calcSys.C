void calcSys()
{
  // PART 1 - fine-binned results (standard)
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  TGraphErrors **graph_lth = new TGraphErrors*[6];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  fIndB->Close();
  // 1 - get results with linear lbd_4
  TFile *fIndl4 = new TFile("../Linear_l4/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fIndl4->Get(Form("graph_lambda_J"));
  fIndl4->Close();
  // 2 - get lambda values for each eff model
  TFile *fIndf1 = new TFile("../../../Simult_eff1/PR_fit/files/finalFitRes.root");
  graph_lth[2] = (TGraphErrors*)fIndf1->Get(Form("graph_lambda_J"));
  fIndf1->Close();
  TFile *fIndf2 = new TFile("../../../Simult_eff2/PR_fit/files/finalFitRes.root");
  graph_lth[3] = (TGraphErrors*)fIndf2->Get(Form("graph_lambda_J"));
  fIndf2->Close();
  // 1 - get 2017 and 2018 results
  TFile *fInd17 = new TFile("../../../2017/PR_fit/files/finalFitRes.root");
  graph_lth[4] = (TGraphErrors*)fInd17->Get(Form("graph_lambda_J"));
  fInd17->Close();
  TFile *fInd18 = new TFile("../../../2018/PR_fit/files/finalFitRes.root");
  graph_lth[5] = (TGraphErrors*)fInd18->Get(Form("graph_lambda_J"));
  fInd18->Close();

  TH1F *h_diffY = new TH1F("h_diffY", "h_diffY", nBinspT, pTBins);
  double v_int[2] = {0,0};
  for(int i = 0; i < nBinspT; i++) {
    // 2017 - 2018 case can be calculated w their uncertainties
    double val7 = graph_lth[4]->GetY()[i];
    double val8 = graph_lth[5]->GetY()[i];
    double diffY = (val7 - val8);
    double err7 = graph_lth[4]->GetEY()[i];
    double err8 = graph_lth[5]->GetEY()[i];
    double errY = sqrt(pow(err7, 2) + pow(err8, 2));

    h_diffY->SetBinContent(i+1, diffY);
    h_diffY->SetBinError(i+1, errY); 
    
    double diff78 = val7-graph_lth[0]->GetY()[i];
    v_int[0] += diff78*(pTBins[i+1]-pTBins[i])/(pTBins[nBinspT]-pTBins[0]);
    
    diff78 = val8-graph_lth[0]->GetY()[i];
    v_int[1] += diff78*(pTBins[i+1]-pTBins[i])/(pTBins[nBinspT]-pTBins[0]);
  }
  double sigY = 0.5*(abs(v_int[0])+abs(v_int[1]));
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLeftMargin(0.13);
  
  TF1 *f1 = new TF1("f1", "[0]", pTBins[0], pTBins[nBinspT]);
  f1->SetParameter(0, 0.1);
  h_diffY->Fit(f1);
  double dY_val = abs(f1->GetParameter(0)/2);
  c->Clear();
  
  // PART 2 - coarse-binned results (to expand)
  // get the histo limits
  TFile *fIn_c = new TFile("../deltaR/files/chistStore.root");
  TH2D* rHist_c;
  fIn_c->GetObject("cHistB", rHist_c);
  
  int nBinspT_c = rHist_c->GetNbinsY();
  const double *pTBins_c = rHist_c->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  TGraphErrors **graph_lth_c = new TGraphErrors*[3];
  // 0 - get Run2 results
  TFile *fIndB_c = new TFile("../deltaR/files/finalFitRes.root");
  graph_lth_c[0] = (TGraphErrors*)fIndB_c->Get("graph_lambda_B");
  graph_lth_c[1] = (TGraphErrors*)fIndB_c->Get(Form("graph_lambda_T"));
  graph_lth_c[2] = (TGraphErrors*)fIndB_c->Get(Form("graph_lambda_L"));
  fIndB_c->Close();

  for(int i = 0; i < nBinspT_c; i++) {
    cout << i << endl;
    cout << graph_lth_c[2]->GetY()[i]-graph_lth_c[0]->GetY()[i] << endl;
    cout << graph_lth_c[1]->GetY()[i]-graph_lth_c[0]->GetY()[i] << endl;
  }

  // expand to the 17 pT bins
  double diff[nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    double pt = graph_lth[0]->GetX()[i];
    for(int j = 0; j < nBinspT_c; j++) {
      double pt_min = graph_lth_c[0]->GetX()[j]-graph_lth_c[0]->GetEX()[j];
      double pt_max = graph_lth_c[0]->GetX()[j]+graph_lth_c[0]->GetEX()[j];
      if(pt > pt_min && pt < pt_max) {
	if(pt < 66) diff[i] = graph_lth_c[1]->GetY()[j]-graph_lth_c[0]->GetY()[j];
	else diff[i] = graph_lth_c[2]->GetY()[j]-graph_lth_c[0]->GetY()[j];
      }
    }
  }
  
  // get the differences - positive and negative
  // linear lambda_4 -> only positive
  TH1F *hP_sigL4 = new TH1F("h_sigL4", "h_sigL4", nBinspT, pTBins);
  // efficiency curves -> symmetrize, pos+neg
  TH1F *hP_sigEff = new TH1F("hP_sigEff", "hP_sigEff", nBinspT, pTBins);
  TH1F *hN_sigEff = new TH1F("hN_sigEff", "hN_sigEff", nBinspT, pTBins);
  // 2017-2018 -> half of the average difference, symmetrize
  TH1F *hP_sigD = new TH1F("hP_sigD", "hP_sigD", nBinspT, pTBins);
  TH1F *hN_sigD = new TH1F("hN_sigD", "hN_sigD", nBinspT, pTBins);
  // deltaR cuts -> one per pT region, symmetrize
  TH1F *hP_sigdR = new TH1F("hP_sigdR", "hP_sigdR", nBinspT, pTBins);
  TH1F *hN_sigdR = new TH1F("hN_sigdR", "hN_sigdR", nBinspT, pTBins);

  double aux;
  for(int i = 0; i < nBinspT; i++) {
    // lambda_4 lin
    hP_sigL4->SetBinContent(i+1, abs(graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i]));
    // muon eff curves
    if(graph_lth[0]->GetX()[i] < 51)
      aux = 0.5*(abs(graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i])+abs(graph_lth[3]->GetY()[i] - graph_lth[0]->GetY()[i]));
    else aux = 0;
    hP_sigEff->SetBinContent(i+1, aux);
    hN_sigEff->SetBinContent(i+1,-aux);
    // 2017 - 2018
    hP_sigD->SetBinContent(i+1, sigY);
    hN_sigD->SetBinContent(i+1,-sigY);
    // dR cut
    aux = abs(diff[i]);
    hP_sigdR->SetBinContent(i+1, aux);
    hN_sigdR->SetBinContent(i+1,-aux);
  }

  // draw the fit results

  // FIRST - draw the abs diff + Simult unc band
  double da_lim = 0.15;

  hP_sigL4->SetFillColor(kRed);
  hP_sigL4->SetLineColor(kRed);

  hP_sigEff->SetFillColor(kGreen+1);
  hP_sigEff->SetLineColor(kGreen+1);
  hN_sigEff->SetFillColor(kGreen+1);
  hN_sigEff->SetLineColor(kGreen+1);

  hP_sigD->SetFillColor(kMagenta);
  hP_sigD->SetLineColor(kMagenta);
  hN_sigD->SetFillColor(kMagenta);
  hN_sigD->SetLineColor(kMagenta);

  hP_sigdR->SetFillColor(kOrange+1);
  hP_sigdR->SetLineColor(kOrange+1);
  hN_sigdR->SetFillColor(kOrange+1);
  hN_sigdR->SetLineColor(kOrange+1);
  
  THStack *hs_pos = new THStack("hs_pos", "Stacked syst uncertainties");
  hs_pos->Add(hP_sigD);
  hs_pos->Add(hP_sigdR); 
  hs_pos->Add(hP_sigEff);
  hs_pos->Add(hP_sigL4);
  hs_pos->SetMinimum(-da_lim);
  hs_pos->SetMaximum(da_lim);
  hs_pos->Draw();

  THStack *hs_neg = new THStack("hs_neg", "Title here");
  hs_neg->Add(hN_sigD);
  hs_neg->Add(hN_sigdR);
  hs_neg->Add(hN_sigEff);
  hs_neg->SetMinimum(-da_lim);
  hs_neg->SetMaximum(da_lim);
  hs_neg->Draw("same");

  TLine *trans1A = new TLine(46, -da_lim, 46, da_lim);
  trans1A->SetLineColor(kBlack);
  trans1A->SetLineStyle(kDashed);
  trans1A->Draw();
  TLine *trans2A = new TLine(66, -da_lim, 66, da_lim);
  trans2A->SetLineColor(kBlack);
  trans2A->SetLineStyle(kDashed);
  trans2A->Draw();

  TLegend *leg = new TLegend(0.65, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(hP_sigL4, "Linear #lambda_{4}", "pl");
  leg->AddEntry(hP_sigEff, "Single #mu eff*", "pl");
  leg->AddEntry(hP_sigdR, "#rho factor", "pl");
  leg->AddEntry(hP_sigD, "2017-2018", "pl");
  leg->Draw();

  c->SaveAs("lth_sigmas.pdf");
  c->Clear();

  // same thing but now squared and stacked
  // positive conts: eff up to 46; lin lambda4; rho factor; constant 2017-2018
  // negative cont: stat unc
  
  // get the systs
  TH1F *f_sigL4 = new TH1F("f_sigL4", "f_sigL4", nBinspT, pTBins);
  TH1F *f_sigEff = new TH1F("f_sigEff", "f_sigEff", nBinspT, pTBins);
  TH1F *f_sigD = new TH1F("f_sigD", "f_sigD", nBinspT, pTBins);
  TH1F *f_sigdR = new TH1F("f_sigdR", "f_sigdR", nBinspT, pTBins);
  
  for(int i = 0; i < nBinspT; i++) {
    // lambda_4 lin
    f_sigL4->SetBinContent(i+1, pow(hP_sigL4->GetBinContent(i+1), 2));
    // eff (symmetric)
    f_sigEff->SetBinContent(i+1, pow(hP_sigEff->GetBinContent(i+1), 2));
    // 2017-2018 (symmetric)
    f_sigD->SetBinContent(i+1, pow(hP_sigD->GetBinContent(i+1), 2));
    // dR cut
    f_sigdR->SetBinContent(i+1, pow(hP_sigdR->GetBinContent(i+1), 2));
  }

  // get the stats
  TH1F *f_stat = new TH1F("f_stat", "f_stat", nBinspT, pTBins);
  for(int i = 0; i < nBinspT; i++) {
    f_stat->SetBinContent(i+1, -pow(graph_lth[0]->GetEY()[i], 2));
  }
  
  // draw stacks -> 2 cases (positive includes linear l4
  da_lim = 0.006;
  
  f_sigL4->SetFillColor(kRed);
  f_sigL4->SetLineColor(kRed);
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

  THStack *hsigN = new THStack("hsigN", "Stacked #sigma^{2}, negative");
  hsigN->SetMinimum(-da_lim);
  hsigN->Add(f_sigD);
  hsigN->Add(f_sigdR);
  hsigN->Add(f_sigEff);
  hsigN->SetMaximum(da_lim);
  THStack *hsigP = new THStack("hsigP", "Stacked #sigma^{2}, positive");
  hsigP->SetMinimum(-da_lim);
  hsigP->Add(f_sigD);
  hsigP->Add(f_sigdR);
  hsigP->Add(f_sigL4);
  hsigP->Add(f_sigEff);
  hsigP->SetMaximum(da_lim);
  
  // negative
  hsigN->Draw();
  hsigN->GetXaxis()->SetTitle("p_{T} (GeV)");
  hsigN->GetYaxis()->SetTitle("#sigma^{2}");
  hsigN->GetYaxis()->SetTitleOffset(2.);
  hstat->Draw("same");
  leg->Draw();

  c->SaveAs("lth_uncs_neg.pdf");
  c->Clear();

  // positive
  hsigP->Draw();
  hsigP->GetXaxis()->SetTitle("p_{T} (GeV)");
  hsigP->GetYaxis()->SetTitle("#sigma^{2}");
  hsigP->GetYaxis()->SetTitleOffset(2.);
  hstat->Draw("same");
  leg->Draw();
  
  c->SaveAs("lth_uncs_pos.pdf");
  c->Clear();  
  c->Destructor();
    
  // tex table with sys uncerts per pt bin
  ofstream ftex;
  ftex.open(Form("sys_unc.tex"));
  ftex << "\\begin{tabular}{c||c|c|c|c||c|c}\n";
  ftex << "$\\pt$ (GeV) & $\\sigma^{\\text{comb}}$ & $\\sigma^{\\text{fit}}$ & $\\sigma^{\\text{eff}}$ & $\\sigma^\\rho$ & $\\sigma_{\\text{sys}}$ & $\\sigma_{\\text{stat}}$ \\\\\n";
  ftex << "\\hline\n";

  int p_norm = 4;
  double val_d, val_l4, val_eff, val_rho;
  int rho_lims[4] = {0, 7, 11, 17};
  for(int i = 0; i < nBinspT; i++) {
    // pT bin
    ftex << Form("$[%.0f, %.0f]$", pTBins[i], pTBins[i+1]);
    // syst sources (start with fixed 3 decimal places
    ftex << setprecision(p_norm) << fixed;
    // combined years: same for all pT
    if(i == 0) {
      val_d = hP_sigD->GetBinContent(i+1);
      ftex << Form(" & \\multirow{%d}{*}{$\\pm", nBinspT) << val_d << "$}";
    }
    else ftex << " & ";
    // linear lambda_4: set to zero if < 0.001
    val_l4 = hP_sigL4->GetBinContent(i+1);
    if(val_l4 < 0.001) {
      val_l4 = 0;
      ftex << " & $-$";
    }
    else 
      ftex << " & $+" << val_l4 << "$";
    // single muon efficiency: already set to zero above
    val_eff = hP_sigEff->GetBinContent(i+1);
    if(val_eff == 0)
      ftex << " & $-$";
    else 
      ftex << " & $\\pm" << val_eff << "$";
    // rho efficiency: constant in 3 pT bins
    int ct = 0;
    for(int j = 0; j < 3; j++) {
      if(i == rho_lims[j]) {
	ct++;
	val_rho = hP_sigdR->GetBinContent(i+1);
	ftex << Form(" & \\multirow{%d}{*}{$\\pm", rho_lims[j+1]-rho_lims[j]) << val_rho << "$}";
      }
    }
    if(ct == 0) ftex << " & ";
    // full syst (sum all sources)
    double valP = sqrt(pow(val_d,2)+pow(val_l4,2)+pow(val_eff,2)+pow(val_rho,2));
    double valN = sqrt(pow(val_d,2)+pow(val_eff,2)+pow(val_rho,2));
    ftex << " & $^{+" << valP << "}_{-" << valN << "}$";
    ftex << " & $\\pm" << graph_lth[0]->GetEY()[i] << "$";
    ftex << "\\\\";
    for(int j = 0; j < 3; j++) 
      if(i+1 == rho_lims[j]) 
	ftex << "\\cline{5-5}";
    ftex << "\n";
  }
  ftex << "\\end{tabular}\n";
  ftex.close();

}
