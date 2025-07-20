// code to plot the fit results for both states simultaneously
// includes both stat and sys uncertainties

void compStates_full()
{
  // get the fit results
  // get lambda values for each bin
  const int n_lbl = 2;
  string lbl[] = {"PR", "NP"}; 
  string nm_lbl[] = {"prompt", "non-prompt"};
  
  TFile *fIndJ = new TFile("../Jpsi/Simult/Systematics/mainDiffs/files/finalUnc.root");
  TGraphErrors **graph_lthJ = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lthJ[i_t] = (TGraphErrors*)fIndJ->Get(Form("lth_f%s", lbl[i_t].c_str()));
  }    
  fIndJ->Close();
  TFile *fIndP = new TFile("../Psi2/Simult/Systematics/mainDiffs/files/finalUnc.root");
  TGraphErrors **graph_lthP = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lthP[i_t] = (TGraphErrors*)fIndP->Get(Form("lth_f%s", lbl[i_t].c_str()));
  }    
  fIndP->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(15, -1, 125, 1);
  fl->SetXTitle("#it{p}_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.4);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->GetXaxis()->SetTitleOffset(1.1);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->GetXaxis()->CenterTitle(true);
  
  int col_j[2] = {kBlue, kRed};
  int col_p[2] = {kViolet+1, kPink+7};

  TLegend *legPR = new TLegend(0.65, 0.85, 0.95, 0.95);
  legPR->SetTextSize(0.03);
  legPR->SetBorderSize(0);
  legPR->SetFillColorAlpha(kWhite,0);
  TLegend *legNP = new TLegend(0.65, 0.2, 0.95, 0.3);
  legNP->SetTextSize(0.03);
  legNP->SetBorderSize(0);
  legNP->SetFillColorAlpha(kWhite,0);

  for(int i_l = 0; i_l < 2; i_l++) {
    
    graph_lthJ[i_l]->SetMarkerSize(.75);
    graph_lthJ[i_l]->SetMarkerStyle(20);
    graph_lthJ[i_l]->SetMarkerColor(col_j[i_l]);
    graph_lthJ[i_l]->SetLineColor(col_j[i_l]);
    graph_lthJ[i_l]->Draw("p same");
    
    graph_lthP[i_l]->SetMarkerStyle(25);
    graph_lthP[i_l]->SetMarkerSize();
    graph_lthP[i_l]->SetMarkerColor(col_p[i_l]);
    graph_lthP[i_l]->SetLineColor(col_p[i_l]);
    graph_lthP[i_l]->Draw("p same");
  }
  legPR->AddEntry(graph_lthJ[0], "prompt J/#psi", "pe");
  legPR->AddEntry(graph_lthP[0], "prompt #psi(2S)", "pe");
  legPR->Draw();  
  legNP->AddEntry(graph_lthJ[1], "non-prompt J/#psi", "pe");
  legNP->AddEntry(graph_lthP[1], "non-prompt #psi(2S)", "pe");
  legNP->Draw();  

  TLine *zero = new TLine(15, 0, 125, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  c->SaveAs("par_lthF.pdf");
  c->Clear();

  // now plot states separately
  // PR J/psi
  TH1F *flJ_PR = c->DrawFrame(15, -1, 125, 1);
  flJ_PR->SetXTitle("p_{T} (GeV)");
  flJ_PR->SetYTitle("#lambda_{#theta}");
  flJ_PR->GetYaxis()->SetTitleOffset(1.3);
  flJ_PR->GetYaxis()->SetLabelOffset(0.01);

  graph_lthJ[0]->SetMarkerSize(.75);
  graph_lthJ[0]->SetMarkerStyle(20);
  graph_lthJ[0]->SetMarkerColor(col_j[0]);
  graph_lthJ[0]->SetLineColor(col_j[0]);
  graph_lthJ[0]->Draw("p same");
  
  zero->Draw();

  c->SaveAs("par_lthPR_jpsi.pdf");
  c->Clear();

  // NP J/psi
  TH1F *flJ_NP = c->DrawFrame(15, -1, 125, 1);
  flJ_NP->SetXTitle("p_{T} (GeV)");
  flJ_NP->SetYTitle("#lambda_{#theta}");
  flJ_NP->GetYaxis()->SetTitleOffset(1.3);
  flJ_NP->GetYaxis()->SetLabelOffset(0.01);

  graph_lthJ[1]->SetMarkerSize(.75);
  graph_lthJ[1]->SetMarkerStyle(20);
  graph_lthJ[1]->SetMarkerColor(col_j[1]);
  graph_lthJ[1]->SetLineColor(col_j[1]);
  graph_lthJ[1]->Draw("p same");
  
  zero->Draw();

  c->SaveAs("par_lthNP_jpsi.pdf");
  c->Clear();

  // PR psi(2S)
  TH1F *flP_PR = c->DrawFrame(15, -1, 125, 1);
  flP_PR->SetXTitle("p_{T} (GeV)");
  flP_PR->SetYTitle("#lambda_{#theta}");
  flP_PR->GetYaxis()->SetTitleOffset(1.3);
  flP_PR->GetYaxis()->SetLabelOffset(0.01);

  graph_lthP[0]->SetMarkerSize(.75);
  graph_lthP[0]->SetMarkerStyle(20);
  graph_lthP[0]->SetMarkerColor(col_p[0]);
  graph_lthP[0]->SetLineColor(col_p[0]);
  graph_lthP[0]->Draw("p same");
  
  zero->Draw();

  c->SaveAs("par_lthPR_psip.pdf");
  c->Clear();

  // NP psi(2S)
  TH1F *flP_NP = c->DrawFrame(15, -1, 125, 1);
  flP_NP->SetXTitle("p_{T} (GeV)");
  flP_NP->SetYTitle("#lambda_{#theta}");
  flP_NP->GetYaxis()->SetTitleOffset(1.3);
  flP_NP->GetYaxis()->SetLabelOffset(0.01);

  graph_lthP[1]->SetMarkerSize(.75);
  graph_lthP[1]->SetMarkerStyle(20);
  graph_lthP[1]->SetMarkerColor(col_p[1]);
  graph_lthP[1]->SetLineColor(col_p[1]);
  graph_lthP[1]->Draw("p same");
  
  zero->Draw();

  c->SaveAs("par_lthNP_psip.pdf");
  c->Clear();
c->Destructor();
}
