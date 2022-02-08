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
  TGraphErrors **graph_lth = new TGraphErrors*[5];
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
  // 3 - get results with higher pT cut
  TFile *fIndpt = new TFile("../../../Simult_cutPt/PR_fit/files/finalFitRes.root");
  graph_lth[4] = (TGraphErrors*)fIndpt->Get(Form("graph_lambda_J"));
  fIndpt->Close();

  // PART 2 - coarse-binned results (to expand)
  // get the histo limits
  TFile *fIn_c = new TFile("../deltaR/files/chistStore.root");
  TH2D* rHist_c;
  fIn_c->GetObject("cHistB", rHist_c);
  
  int nBinspT_c = rHist_c->GetNbinsY();
  const double *pTBins_c = rHist_c->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - 3 sets
  TGraphErrors **graph_lth_c = new TGraphErrors*[3];
  // 0 - get Run2 results
  TFile *fIndB_c = new TFile("../deltaR/files/finalFitRes.root");
  graph_lth_c[0] = (TGraphErrors*)fIndB_c->Get("graph_lambda_B");
  graph_lth_c[1] = (TGraphErrors*)fIndB_c->Get(Form("graph_lambda_T"));
  graph_lth_c[2] = (TGraphErrors*)fIndB_c->Get(Form("graph_lambda_L"));
  fIndB_c->Close();

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
  
  // get the differences
  TH1F *g_sigL4 = new TH1F("h_sigL4", "h_sigL4", nBinspT, pTBins);
  TH1F *g_sigF = new TH1F("h_sigF", "h_sigF", nBinspT, pTBins);
  TH1F *g_sigFI = new TH1F("h_sigFI", "h_sigFI", nBinspT, pTBins);
  TH1F *g_sigpT = new TH1F("h_sigpT", "h_sigpT", nBinspT, pTBins);
  TH1F *g_sigR = new TH1F("h_sigR", "h_sigR", nBinspT, pTBins);

  TH1F *g2_sigL4 = new TH1F("h2_sigL4", "h2_sigL4", nBinspT, pTBins);
  TH1F *g2_sigF = new TH1F("h2_sigF", "h2_sigF", nBinspT, pTBins);
  TH1F *g2_sigFI = new TH1F("h2_sigFI", "h2_sigFI", nBinspT, pTBins);
  TH1F *g2_sigpT = new TH1F("h2_sigpT", "h2_sigpT", nBinspT, pTBins);
  TH1F *g2_sigR = new TH1F("h2_sigR", "h2_sigR", nBinspT, pTBins);
  for(int i = 0; i < nBinspT; i++) {
    // lambda_4 lin
    if(graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i]>0)
      g_sigL4->SetBinContent(i+1, graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i]);
    else
      g2_sigL4->SetBinContent(i+1, graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i]);
    // weigh by f(pT)
    if(graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i]>0)
      g_sigF->SetBinContent(i+1, graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i]);
    else
      g2_sigF->SetBinContent(i+1, graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i]);
    // weigh by 1/f(pT)
    if(graph_lth[3]->GetY()[i] - graph_lth[0]->GetY()[i]>0)
      g_sigFI->SetBinContent(i+1, graph_lth[3]->GetY()[i] - graph_lth[0]->GetY()[i]);
    else
      g2_sigFI->SetBinContent(i+1, graph_lth[3]->GetY()[i] - graph_lth[0]->GetY()[i]);
    // pT cut
    if(graph_lth[4]->GetY()[i] - graph_lth[0]->GetY()[i]>0)
      g_sigpT->SetBinContent(i+1, graph_lth[4]->GetY()[i] - graph_lth[0]->GetY()[i]);
    else
      g2_sigpT->SetBinContent(i+1, graph_lth[4]->GetY()[i] - graph_lth[0]->GetY()[i]);
    // dR cut
    if(diff[i] > 0)
      g_sigR->SetBinContent(i+1, diff[i]);
    else
      g2_sigR->SetBinContent(i+1, diff[i]);
  }

  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);

  // FIRST - draw the abs diff + Simult unc band
  double da_lim = 0.06;
  
  g_sigF->SetFillColor(kBlack);
  g_sigFI->SetFillColor(kBlue);
  g_sigpT->SetFillColor(kGreen+1);
  g_sigL4->SetFillColor(kRed);
  g_sigR->SetFillColor(kOrange+2);
  g_sigF->SetLineColor(kBlack);
  g_sigFI->SetLineColor(kBlue);
  g_sigpT->SetLineColor(kGreen+1);
  g_sigL4->SetLineColor(kRed);
  g_sigR->SetLineColor(kOrange+2);

  g2_sigF->SetFillColor(kBlack);
  g2_sigFI->SetFillColor(kBlue);
  g2_sigpT->SetFillColor(kGreen+1);
  g2_sigL4->SetFillColor(kRed);
  g2_sigR->SetFillColor(kOrange+2);
  g2_sigF->SetLineColor(kBlack);
  g2_sigFI->SetLineColor(kBlue);
  g2_sigpT->SetLineColor(kGreen+1);
  g2_sigL4->SetLineColor(kRed);
  g2_sigR->SetLineColor(kOrange+2);
  
  THStack *hs_pos = new THStack("hs_pos", "Title here");
  hs_pos->Add(g_sigF);
  hs_pos->Add(g_sigFI);
  hs_pos->Add(g_sigL4);
  hs_pos->Add(g_sigpT);
  hs_pos->Add(g_sigR);
  hs_pos->SetMinimum(-da_lim);
  hs_pos->SetMaximum(da_lim);
  hs_pos->Draw();

  THStack *hs_neg = new THStack("hs_neg", "Title here");
  hs_neg->Add(g2_sigF);
  hs_neg->Add(g2_sigFI);
  hs_neg->Add(g2_sigL4);
  hs_neg->Add(g2_sigpT);
  hs_neg->Add(g2_sigR);
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
  leg->AddEntry(g_sigF, "f(p_{T})", "pl");
  leg->AddEntry(g_sigFI, "1/f(p_{T})", "pl");
  leg->AddEntry(g_sigpT, "p_{T} cut", "pl");
  leg->AddEntry(g_sigL4, "Lin #lambda_{4}", "pl");
  leg->AddEntry(g_sigR, "#rho factor", "pl");
  leg->Draw();

  c->SaveAs("lth_sigmas.pdf");
  c->Clear();

  // same thing but now squared and stacked
  // positive conts: pT cut up to 46; lin lambda4; rho factor from 46
  // negative cont: stat unc
  
  // get the systs
  TH1F *f_sigL4 = new TH1F("f_sigL4", "f_sigL4", nBinspT, pTBins);
  TH1F *f_sigpT = new TH1F("f_sigpT", "f_sigpT", nBinspT, pTBins);
  TH1F *f_sigR = new TH1F("f_sigR", "f_sigR", nBinspT, pTBins);

  for(int i = 0; i < nBinspT; i++) {
    // lambda_4 lin
    f_sigL4->SetBinContent(i+1, pow(graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i], 2));
    // pT cut (up to 46 GeV)
    if(graph_lth[0]->GetX()[i] < 46)
      f_sigpT->SetBinContent(i+1, pow(graph_lth[4]->GetY()[i] - graph_lth[0]->GetY()[i], 2));
    else
      f_sigpT->SetBinContent(i+1, 0);
    // dR cut
    if(graph_lth[0]->GetX()[i] > 46)
      f_sigR->SetBinContent(i+1, pow(diff[i], 2));
    else
      f_sigR->SetBinContent(i+1, 0);
  }

  // get the stats
  TH1F *f_stat = new TH1F("f_stat", "f_stat", nBinspT, pTBins);
  for(int i = 0; i < nBinspT; i++) {
    f_stat->SetBinContent(i+1, -pow(graph_lth[0]->GetEY()[i], 2));
  }
  
  // draw stacks
  da_lim = 0.005;
  
  f_sigpT->SetFillColor(kGreen+1);
  f_sigL4->SetFillColor(kRed);
  f_sigR->SetFillColor(kOrange+2);
  f_sigpT->SetLineColor(kGreen+1);
  f_sigL4->SetLineColor(kRed);
  f_sigR->SetLineColor(kOrange+2);

  f_stat->SetFillColor(kBlue);
  f_stat->SetLineColor(kBlue);

  THStack *hsig = new THStack("hsig", "Title here");
  hsig->Add(f_sigpT);
  hsig->Add(f_sigR);
  hsig->Add(f_sigL4);
  hsig->SetMinimum(-da_lim);
  hsig->SetMaximum(da_lim);
  hsig->Draw();

  THStack *hstat = new THStack("hstat", "Title here");
  hstat->Add(f_stat);
  hstat->SetMinimum(-da_lim);
  hstat->SetMaximum(da_lim);
  hstat->Draw("same");

  TLine *trans1F = new TLine(46, -da_lim, 46, da_lim);
  trans1F->SetLineColor(kBlack);
  trans1F->SetLineStyle(kDashed);
  trans1F->Draw();
  TLine *trans2F = new TLine(66, -da_lim, 66, da_lim);
  trans2F->SetLineColor(kBlack);
  trans2F->SetLineStyle(kDashed);
  trans2F->Draw();

  TLegend *legF = new TLegend(0.65, 0.7, 0.9, 0.9);
  legF->SetTextSize(0.03);
  legF->AddEntry(f_sigpT, "p_{T} cut", "pl");
  legF->AddEntry(f_sigL4, "Lin #lambda_{4}", "pl");
  legF->AddEntry(f_sigR, "#rho factor", "pl");
  legF->AddEntry(f_stat, "stat", "pl");
  legF->Draw();

  c->SaveAs("lth_uncs.pdf");
  c->Clear();

  
  c->Destructor();

}
