void plot_Clim()
{
  // get clim values from Jpsi
  TFile *fin_J = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/PR_fit/files/finalFitRes.root");
  TGraph *g_CJ = (TGraph*)fin_J->Get("graph_cMax");
  fin_J->Close();

  // get clim values from Psi(2S)
  TFile *fin_P = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Simult/PR_fit/files/finalFitRes.root");
  TGraph *g_CP = (TGraph*)fin_P->Get("graph_cMax");
  fin_P->Close();

  double pmin = g_CP->GetX()[0]-g_CP->GetEX()[0];
  int nmax = g_CJ->GetN();
  double pmax = g_CJ->GetX()[nmax-1]+g_CJ->GetEX()[nmax-1];

  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);

  TH1F *fc = c->DrawFrame(pmin-5, 0, pmax+5, 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("|cos#theta|_{max}");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);

  g_CJ->SetMarkerColor(kBlack);
  g_CJ->SetLineColor(kBlack);
  g_CJ->SetMarkerStyle(20);
  g_CJ->SetMarkerSize(.5);
  g_CJ->Draw("psame");

  g_CP->SetMarkerColor(kBlue);
  g_CP->SetLineColor(kBlue);
  g_CP->SetMarkerStyle(25);
  g_CP->SetMarkerSize(.5);
  g_CP->Draw("psame");

  TLegend *leg = new TLegend(0.65, 0.5, 0.95, 0.6);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  leg->AddEntry(g_CJ, "J/#psi", "pl");
  leg->AddEntry(g_CP, "#psi(2S)", "pl");
  leg->Draw();

  c->SaveAs("cosMax.pdf");
  c->Clear();
  c->Destructor();
}
