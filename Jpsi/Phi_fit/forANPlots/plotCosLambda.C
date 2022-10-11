void plotCosLambda()
{
  // get the LSB and RSB lamda results with everything free
  TGraphErrors **g_l2 = new TGraphErrors*[4];
  TGraphErrors **g_l4 = new TGraphErrors*[4];
  
  TFile *finL = new TFile("../PR_fit/files/LSB_fitres.root");
  finL->GetObject("fit_0_l2", g_l2[0]);
  finL->GetObject("fit_0_l4", g_l4[0]);
  finL->GetObject("fit_1_l2", g_l2[2]);
  finL->GetObject("fit_1_l4", g_l4[2]);
  finL->Close();

  TFile *finR = new TFile("../PR_fit/files/RSB_fitres.root");
  finR->GetObject("fit_0_l2", g_l2[1]);
  finR->GetObject("fit_0_l4", g_l4[1]);
  finR->GetObject("fit_1_l2", g_l2[3]);
  finR->GetObject("fit_1_l4", g_l4[3]);
  finR->Close();

  int n = g_l2[0]->GetN();
  double xmin = g_l2[0]->GetX()[0]-g_l2[0]->GetEX()[0]-5;
  double xmax = g_l2[0]->GetX()[n-1]+g_l2[0]->GetEX()[n-1]+5;
  
  // plot the results for each SB and each lambda
  TCanvas *c = new TCanvas("", "", 900, 900);
  
  TH1F *fl1 = c->DrawFrame(xmin, -2.99, xmax, 2.99);
  fl1->SetXTitle("p_{T} (GeV)");
  fl1->SetYTitle("#lambda_{2}");
  fl1->GetYaxis()->SetTitleOffset(1.3);
  fl1->GetYaxis()->SetLabelOffset(0.01);
  //fl1->SetTitle("Run 2 #lambda_{#theta}");

  g_l2[0]->SetMarkerStyle(20);
  g_l2[0]->SetMarkerSize(0.75);
  g_l2[0]->Draw("psame");

  TLatex cLl2;
  cLl2.DrawLatex(xmin+5, 2.5, "LSB #lambda_{2}");

  c->SaveAs("plots/bkgCosth/LSB_l2.pdf");
  c->Clear();

  TH1F *fl2 = c->DrawFrame(xmin, -2.99, xmax, 2.99);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#lambda_{4}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  //fl2->SetTitle("Run 2 #lambda_{#theta}");

  g_l4[0]->SetMarkerStyle(20);
  g_l4[0]->SetMarkerSize(0.75);
  g_l4[0]->Draw("psame");

  TLatex cLl4;
  cLl4.DrawLatex(xmin+5, 2.5, "LSB #lambda_{4}");

  c->SaveAs("plots/bkgCosth/LSB_l4.pdf");
  c-> Clear();
  
  TH1F *fl3 = c->DrawFrame(xmin, -2.99, xmax, 2.99);
  fl3->SetXTitle("p_{T} (GeV)");
  fl3->SetYTitle("#lambda_{2}");
  fl3->GetYaxis()->SetTitleOffset(1.3);
  fl3->GetYaxis()->SetLabelOffset(0.01);
  //fl3->SetTitle("Run 2 #lambda_{#theta}");

  g_l2[1]->SetMarkerStyle(20);
  g_l2[1]->SetMarkerSize(0.75);
  g_l2[1]->Draw("psame");

  TLatex cRl2;
  cRl2.DrawLatex(xmin+5, 2.5, "RSB #lambda_{2}");

  c->SaveAs("plots/bkgCosth/RSB_l2.pdf");
  c->Clear();

  TH1F *fl4 = c->DrawFrame(xmin, -2.99, xmax, 2.99);
  fl4->SetXTitle("p_{T} (GeV)");
  fl4->SetYTitle("#lambda_{4}");
  fl4->GetYaxis()->SetTitleOffset(1.3);
  fl4->GetYaxis()->SetLabelOffset(0.01);
  //fl2->SetTitle("Run 2 #lambda_{#theta}");

  g_l4[1]->SetMarkerStyle(20);
  g_l4[1]->SetMarkerSize(0.75);
  g_l4[1]->Draw("psame");

  TLatex cRl4;
  cRl4.DrawLatex(xmin+5, 2.5, "RSB #lambda_{4}");

  c->SaveAs("plots/bkgCosth/RSB_l4.pdf");
  c->Clear();

  // now plot the results with lambda_2 fixed
  TH1F *flf1 = c->DrawFrame(xmin, -2.99, xmax, 2.99);
  flf1->SetXTitle("p_{T} (GeV)");
  flf1->SetYTitle("#lambda_{2}");
  flf1->GetYaxis()->SetTitleOffset(1.3);
  flf1->GetYaxis()->SetLabelOffset(0.01);
  //flf1->SetTitle("Run 2 #lambda_{#theta}");

  g_l2[2]->SetMarkerStyle(20);
  g_l2[2]->SetMarkerSize(0.75);
  g_l2[2]->SetMarkerColor(kBlue);
  g_l2[2]->SetLineColor(kBlue);
  g_l2[2]->Draw("psame");

  TLatex cLlf2;
  cLlf2.DrawLatex(xmin+5, 2.5, "LSB #lambda_{2}");

  c->SaveAs("plots/bkgCosth/LSB_fixed_l2.pdf");
  c->Clear();

  TH1F *flf2 = c->DrawFrame(xmin, -2.99, xmax, 2.99);
  flf2->SetXTitle("p_{T} (GeV)");
  flf2->SetYTitle("#lambda_{4}");
  flf2->GetYaxis()->SetTitleOffset(1.3);
  flf2->GetYaxis()->SetLabelOffset(0.01);
  //flf2->SetTitle("Run 2 #lambda_{#theta}");

  g_l4[2]->SetMarkerStyle(20);
  g_l4[2]->SetMarkerSize(0.75);
  g_l4[2]->SetMarkerColor(kBlue);
  g_l4[2]->SetLineColor(kBlue);
  g_l4[2]->Draw("psame");

  TLatex cLlf4;
  cLlf4.DrawLatex(xmin+5, 2.5, "LSB #lambda_{4}");

  c->SaveAs("plots/bkgCosth/LSB_fixed_l4.pdf");
  c-> Clear();
  
  TH1F *flf3 = c->DrawFrame(xmin, -2.99, xmax, 2.99);
  flf3->SetXTitle("p_{T} (GeV)");
  flf3->SetYTitle("#lambda_{2}");
  flf3->GetYaxis()->SetTitleOffset(1.3);
  flf3->GetYaxis()->SetLabelOffset(0.01);
  //flf3->SetTitle("Run 2 #lambda_{#theta}");

  g_l2[3]->SetMarkerStyle(20);
  g_l2[3]->SetMarkerSize(0.75);
  g_l2[3]->SetMarkerColor(kBlue);
  g_l2[3]->SetLineColor(kBlue);
  g_l2[3]->Draw("psame");

  TLatex cRlf2;
  cRlf2.DrawLatex(xmin+5, 2.5, "RSB #lambda_{2}");

  c->SaveAs("plots/bkgCosth/RSB_fixed_l2.pdf");
  c->Clear();

  TH1F *flf4 = c->DrawFrame(xmin, -2.99, xmax, 2.99);
  flf4->SetXTitle("p_{T} (GeV)");
  flf4->SetYTitle("#lambda_{4}");
  flf4->GetYaxis()->SetTitleOffset(1.3);
  flf4->GetYaxis()->SetLabelOffset(0.01);
  //flf2->SetTitle("Run 2 #lambda_{#theta}");

  g_l4[3]->SetMarkerStyle(20);
  g_l4[3]->SetMarkerSize(0.75);
  g_l4[3]->SetMarkerColor(kBlue);
  g_l4[3]->SetLineColor(kBlue);
  g_l4[3]->Draw("psame");

  TLatex cRlf4;
  cRlf4.DrawLatex(xmin+5, 2.5, "RSB #lambda_{4}");

  c->SaveAs("plots/bkgCosth/RSB_fixed_l4.pdf");

  c->Destructor();
}
