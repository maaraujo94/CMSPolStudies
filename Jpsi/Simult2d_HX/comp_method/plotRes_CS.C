// code to plot the fit results - PR case
// compares the 2d and 1d fits for CS (w/ and w/o ltp)

#import "../ptbins.C"

void plotRes_CS()
{
  // get the fit results - 2d
  TFile *fInd = new TFile("../../Simult2d_CS/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lPR = new TGraphErrors*[3];
  string lbl[3] = {"theta", "phi", "thph"};
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lPR[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();
  
  // get the fit results - 2d w/o ltp
  TFile *fIndz = new TFile("../../Simult2d_CS/PR_fit/files/finalFitRes_ltpz.root");
  TGraphErrors **graph_lPRz = new TGraphErrors*[3];
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lPRz[i_t] = (TGraphErrors*)fIndz->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fIndz->Close();

  //get the fit results - 1d
  // lth is immediate
  TFile *find_th = new TFile("../../Simult_CS/PR_fit/files/finalFitRes.root");
  TGraphErrors *graph_lth = (TGraphErrors*)find_th->Get(Form("graph_lambda_J"));
  find_th->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.02);
  int col[] = {kBlack, kRed, kGreen+1};

  // draw lambda_th(pT)
  TH1F *flth = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
  flth->SetXTitle("p_{T} (GeV)");
  flth->SetYTitle("#lambda_{#theta}^{PR}");
  flth->GetYaxis()->SetTitleOffset(1.3);
  flth->GetYaxis()->SetLabelOffset(0.01);
  flth->SetTitle("");
    
  graph_lPR[0]->SetLineColor(col[0]);
  graph_lPR[0]->SetMarkerStyle(20);
  graph_lPR[0]->SetMarkerSize(.5);
  graph_lPR[0]->SetMarkerColor(col[0]);
  graph_lPR[0]->Draw("p same");

  graph_lPRz[0]->SetLineColor(col[0]);
  graph_lPRz[0]->SetMarkerStyle(24);
  graph_lPRz[0]->SetMarkerSize(.5);
  graph_lPRz[0]->SetMarkerColor(col[0]);
  graph_lPRz[0]->SetLineStyle(kDashed);
  graph_lPRz[0]->Draw("p same");

  graph_lth->SetLineColor(col[1]);
  graph_lth->SetMarkerStyle(25);
  graph_lth->SetMarkerSize(.5);
  graph_lth->SetMarkerColor(col[1]);
  graph_lth->SetLineStyle(kDashed);
  graph_lth->Draw("p same");
  
  TLine *zero = new TLine(ptBins[0]-5, 0, ptBins[nPtBins], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *legth = new TLegend(0.2, 0.15, 0.5, 0.3);
  legth->SetTextSize(0.03);
  legth->SetBorderSize(0);
  legth->SetFillColorAlpha(kWhite, 0);
  legth->AddEntry(graph_lPR[0], "#lambda_{#theta}^{HX} (2d)", "pl");
  legth->AddEntry(graph_lPRz[0], "#lambda_{#theta}^{HX} (2d w/o #lambda_{#theta#phi})", "pl");
  legth->AddEntry(graph_lth, "#lambda_{#theta}^{HX} (1d)", "pl");
  legth->Draw();
  
  c->SaveAs("lth_CS_comp.pdf");
  c->Clear();

  c->Destructor();
 
}
