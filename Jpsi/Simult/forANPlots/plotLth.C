// code to plot the fit results

#import "../ptbins.C"

void plotLth()
{
  // get the fit results
  // get lambda values for each bin
  string lbl[] = {"Data", "NP", "PR", "J"};
  TFile *fInd = new TFile("../PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth = new TGraphErrors*[4];
  for(int i_t = 0; i_t < 4; i_t++) {
    graph_lth[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();

  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  int cols[] = {kViolet-1, kRed, kBlack, kBlue, kGreen};

  for(int i = 0; i < 4; i++) {
    graph_lth[i]->SetLineColor(cols[i]);
    graph_lth[i]->SetMarkerColor(cols[i]);
  }
  TLegend *leg = new TLegend(0.7, 0.7, 0.97, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(graph_lth[0], "Peak", "pl");
  leg->AddEntry(graph_lth[1], "non-prompt J/#psi", "pl");
  leg->AddEntry(graph_lth[2], "PR", "pl");
  leg->AddEntry(graph_lth[3], "prompt J/#psi", "pl");
  
  // draw lambda_th(pT) - just peak
  TH1F *fl1 = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#lambda_{#theta}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  fl1->SetTitle("Run 2 #lambda_{#theta}");
  
  TLine *zero = new TLine(ptBins[0]-5, 0, ptBins[nPtBins], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  graph_lth[0]->Draw("p same");

  leg->Draw(); 
  
  c->SaveAs("plots/ratioFinal/lth1_S.pdf");

  // add pr lambda_th
  graph_lth[1]->Draw("p same");
  graph_lth[2]->Draw("p same");
  
  c->SaveAs("plots/ratioFinal/lth2_S.pdf");
  
  // add prompt jpsi lambda_th
  graph_lth[3]->Draw("p same");
  
  c->SaveAs("plots/ratioFinal/lth3_S.pdf");
  c->Clear();

  // now draw just prompt and non-prompt J/psi results
  TH1F *fl2 = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("Run 2 #lambda_{#theta}");

  // prompt J/psi
  graph_lth[3]->Draw("p same");

  // non-prompt J/psi
  graph_lth[1]->Draw("p same");

  zero->Draw();
  
  TLegend *leg2 = new TLegend(0.67, 0.7, 0.97, 0.9);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(graph_lth[3], "prompt J/#psi", "pl");
  leg2->AddEntry(graph_lth[1], "non-prompt J/#psi", "pl");
  leg2->Draw();

  c->SaveAs("plots/ratioFinal/lth_main.pdf");
  c->Clear();
  c->Destructor();
}
