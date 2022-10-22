// code to plot the results with full uncertainty

void plotLth()
{
  // get the fit results
  TGraphErrors **graph_lth = new TGraphErrors*[2];
  TFile *fInd = new TFile("lthUnc.root");
  graph_lth[0] = (TGraphErrors*)fInd->Get("graph_lambda");
  fInd->Close();

  TFile *fIndN = new TFile("lthUnc_NP.root");
  graph_lth[1] = (TGraphErrors*)fIndN->Get("graph_lambda");
  fIndN->Close();
  
  int n_pt = graph_lth[0]->GetN();
  double ptmin = graph_lth[0]->GetX()[0]-graph_lth[0]->GetEX()[0];
  double ptmax = graph_lth[0]->GetX()[n_pt-1]+graph_lth[0]->GetEX()[n_pt-1];
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);

  // now draw just prompt and non-prompt J/psi results
  TH1F *fl2 = c->DrawFrame(ptmin-5, -1, ptmax, 1);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("Run 2 #lambda_{#theta}");

  // prompt J/psi
  graph_lth[0]->SetLineColor(kBlue);
  graph_lth[0]->SetMarkerColor(kBlue);
  graph_lth[0]->Draw("p same");

  // non-prompt J/psi
  graph_lth[1]->SetLineColor(kRed);
  graph_lth[1]->SetMarkerColor(kRed);
  graph_lth[1]->Draw("p same");

  TLine *zero = new TLine(ptmin-5, 0, ptmax, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  
  TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(graph_lth[0], "prompt J/#psi", "pl");
  leg2->AddEntry(graph_lth[1], "non-prompt J/#psi", "pl");
  leg2->Draw();

  c->SaveAs("plots/lth_full_unc_both.pdf");
  c->Destructor();


}
