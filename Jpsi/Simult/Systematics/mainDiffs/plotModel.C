// code to compare the deviations for several scenarios
// 1) free alpha - Run2
// 2) linear mass bkg - Run2
// 3) linear lambda_4 - Run2

void plotModel()
{
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - 4 sets
  TGraphErrors **graph_lth = new TGraphErrors*[4];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  fIndB->Close();
  // 1 - get free alpha results
  TFile *fInda = new TFile("../Alpha_free/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fInda->Get(Form("graph_lambda_J"));
  fInda->Close();
  // 2 - get lambda values for linear mass bkg
  TFile *fIndmb = new TFile("../Mass_bkg_alt/files/finalFitRes.root");
  graph_lth[2] = (TGraphErrors*)fIndmb->Get(Form("graph_lambda_J"));
  fIndmb->Close();
  // 3 - get results with linear lbd_4
  TFile *fIndl4 = new TFile("../Linear_l4/files/finalFitRes.root");
  graph_lth[3] = (TGraphErrors*)fIndl4->Get(Form("graph_lambda_J"));
  fIndl4->Close();

  // get the differences
  double diff[3][nBinspT], za[nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    diff[0][i] = (graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[1][i] = (graph_lth[2]->GetY()[i] - graph_lth[0]->GetY()[i]);
    diff[2][i] = (graph_lth[3]->GetY()[i] - graph_lth[0]->GetY()[i]);
    
    za[i] = 0;
  }
  TGraphErrors *g_lthA = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[0], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthMB = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[1], graph_lth[0]->GetEX(), za);
  TGraphErrors *g_lthL4 = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[2], graph_lth[0]->GetEX(), za);
  
  TGraphErrors *g_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_lth[0]->GetEY());
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  double d_lim = 60;

  // FIRST - draw the abs diff + Simult unc band
  double da_lim = 0.3;
  
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#Delta#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("#Delta#lambda_{#theta} (prompt J/#psi)");
  
  g_lthA->SetLineColor(kBlue);
  g_lthA->SetMarkerColor(kBlue);
  g_lthA->SetMarkerStyle(20);
  g_lthA->SetMarkerSize(.75);
  g_lthA->Draw("p same");

  g_lthMB->SetLineColor(kRed);
  g_lthMB->SetMarkerColor(kRed);
  g_lthMB->SetMarkerStyle(20);
  g_lthMB->SetMarkerSize(.75);
  g_lthMB->Draw("p same");

  g_lthL4->SetLineColor(kGreen+1);
  g_lthL4->SetMarkerColor(kGreen+1);
  g_lthL4->SetMarkerStyle(20);
  g_lthL4->SetMarkerSize(.75);
  g_lthL4->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.3);
  g_unc->Draw("ce3");

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  //zero->Draw();
  TLine *trans1A = new TLine(46, -da_lim, 46, da_lim);
  trans1A->SetLineColor(kBlack);
  trans1A->SetLineStyle(kDashed);
  trans1A->Draw();
  TLine *trans2A = new TLine(66, -da_lim, 66, da_lim);
  trans2A->SetLineColor(kBlack);
  trans2A->SetLineStyle(kDashed);
  trans2A->Draw();

  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(g_lthA, "Free #alpha", "pl");
  leg->AddEntry(g_lthMB, "Linear mass bkg", "pl");
  leg->AddEntry(g_lthL4, "Linear #lambda_{4}", "pl");
  leg->Draw();

  c->SaveAs("lth_absDiff_M.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();


}
