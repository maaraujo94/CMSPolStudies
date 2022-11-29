// code to compare the deviations for several scenarios
// 1) free alpha - Run2
// 2) linear mass bkg - Run2

void plotModel()
{
  // get the histo limits
  TFile *fIn = new TFile("../../PR_fit/files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results - 2 sets
  TGraphErrors **graph_lth = new TGraphErrors*[2];
  // 0 - get Run2 results
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  graph_lth[0] = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  fIndB->Close();
  // 2 - get lambda values for linear mass bkg
  TFile *fIndmb = new TFile("../Mass_bkg_alt/files/finalFitRes.root");
  graph_lth[1] = (TGraphErrors*)fIndmb->Get(Form("graph_lambda_J"));
  fIndmb->Close();

  // get the differences
  double diff[1][nBinspT], za[nBinspT];
  double err[1][nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    diff[0][i] = (graph_lth[1]->GetY()[i] - graph_lth[0]->GetY()[i]);
 
    double unc1 = graph_lth[0]->GetEY()[i];
    double unc2 = graph_lth[1]->GetEY()[i];
    err[0][i] = sqrt(abs(pow(unc1,2)-pow(unc2,2)));

    za[i] = 0;
  }
  TGraphErrors *g_lthMB = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), diff[0], graph_lth[0]->GetEX(), za);
  
  TGraphErrors *g_unc = new TGraphErrors(nBinspT, graph_lth[0]->GetX(), za, graph_lth[0]->GetEX(), graph_lth[0]->GetEY());
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  
  double d_lim = 60;

  // FIRST - draw the abs diff + Simult unc band
  double da_lim = 0.3;
  
  TH1F *fl1 = c->DrawFrame(pTBins[0]-5, -da_lim, pTBins[nBinspT], da_lim);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#Delta#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("#Delta#lambda_{#theta} (mass bkg subtraction)");
  
  g_lthMB->SetLineColor(kRed);
  g_lthMB->SetMarkerColor(kRed);
  g_lthMB->SetMarkerStyle(20);
  g_lthMB->SetMarkerSize(.75);
  g_lthMB->Draw("p same");

  g_unc->SetLineColor(kBlack);
  g_unc->SetFillColorAlpha(kBlack, 0.1);
  g_unc->Draw("ce3");

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  //zero->Draw();
 
  c->SaveAs("plots/lth_absDiff_M.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();


}
