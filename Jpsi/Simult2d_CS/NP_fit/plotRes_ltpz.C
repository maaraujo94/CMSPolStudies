// code to plot the fit results

#import "../ptbins.C"

void plotRes_ltpz()
{
  // get the fit results
  // get A, lambda, chiProb values for each bin
  TFile *fInd = new TFile("files/finalFitRes_ltpz.root");
  TGraphErrors *graph_A = (TGraphErrors*)fInd->Get(Form("graph_A"));
  TGraph *graph_chi = (TGraph*)fInd->Get(Form("graph_chiP"));

  TGraphErrors **graph_lth = new TGraphErrors*[3];
  string lbl[3] = {"theta", "phi", "thph"};
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lth[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.02);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda^{NP}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("");
    
  int col[] = {kBlack, kRed, kGreen+1};
  for(int i = 0; i < 3; i++) {
    graph_lth[i]->SetLineColor(col[i]);
    graph_lth[i]->SetMarkerStyle(20);
    graph_lth[i]->SetMarkerSize(.5);
    graph_lth[i]->SetMarkerColor(col[i]);
    graph_lth[i]->Draw("p same");
  }

  TLine *zero = new TLine(ptBins[0]-5, 0, ptBins[nPtBins], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg = new TLegend(0.65, 0.8, 0.95, 0.95);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite, 0);
  leg->AddEntry(graph_lth[0], "#lambda_{#theta}^{CS}", "pl");
  leg->AddEntry(graph_lth[1], "#lambda_{#phi}^{CS}", "pl");
  leg->AddEntry(graph_lth[2], "#lambda_{#theta#phi}^{CS}", "pl");
  leg->Draw();
  
  c->SaveAs("plots/ratioFinal_ltpz/par_lth.pdf");
  c->Clear();

  // draw A(pT)
  c->SetTopMargin(0.1);
    
  c->SetLogy();
  TH1F *fa = c->DrawFrame(ptBins[0], 1, ptBins[nPtBins], 1e2);
  fa->SetXTitle("p_{T} (GeV)");
  fa->SetYTitle("A");
  fa->GetYaxis()->SetTitleOffset(1.3);
  fa->GetYaxis()->SetLabelOffset(0.01);
  fa->SetTitle("Run 2 A");
 
  graph_A->SetLineColor(col[0]);
  graph_A->SetMarkerColor(col[0]);
  graph_A->Draw("p same");

  c->SaveAs("plots/ratioFinal_ltpz/par_A.pdf");
  c->Clear();

  // draw chiProb(pT)
  c->SetLogy(0);
  TH1F *fc = c->DrawFrame(ptBins[0], 0, ptBins[nPtBins], 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("P(#chi^{2}, ndf)");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);
  fc->SetTitle("Run 2 P(#chi^{2}, ndf)");

  graph_chi->SetLineColor(col[0]);
  graph_chi->SetMarkerColor(col[0]);
  graph_chi->SetMarkerStyle(20);
  graph_chi->SetMarkerSize(.75);
  graph_chi->Draw("p same");

  c->SaveAs("plots/ratioFinal_ltpz/par_chiP.pdf");
  c->Clear();
  c->Destructor();
  


}
