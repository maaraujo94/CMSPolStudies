void compLth()
{
  TFile *fin1 = new TFile("files/finalFitRes.root");
  TGraphErrors *g_lth1 = (TGraphErrors*)fin1->Get("graph_lambda_J");
  fin1->Close();

  TFile *fin2 = new TFile("files/finalFitRes_full.root");
  TGraphErrors *g_lth2 = (TGraphErrors*)fin2->Get("graph_lambda_J");
  fin2->Close();

  TCanvas *c = new TCanvas("", "", 700, 700);

  // draw lambda_th(pT)
  double xmin = g_lth1->GetX()[0]-g_lth1->GetEX()[0];
  double xmax = g_lth1->GetX()[g_lth1->GetN()-1]+g_lth1->GetEX()[g_lth1->GetN()-1];
  
  TH1F *fl = c->DrawFrame(xmin, -1, xmax, 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("Run 2 #lambda_{#theta}");

  g_lth1->SetLineColor(kBlack);
  g_lth1->SetMarkerColor(kBlack);
  g_lth1->Draw("p same");

  g_lth2->SetLineColor(kBlue);
  g_lth2->SetMarkerColor(kBlue);
  g_lth2->Draw("p same");

  TLine *zero = new TLine(xmin, 0, xmax, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(g_lth1, "with model", "pl");
  leg->AddEntry(g_lth2, "without model", "pl");
  leg->Draw();

  c->SaveAs("plots/lth_comp.pdf");
  c->Destructor();
}
