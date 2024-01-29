// macro to plot fitted costh dists for slides
void plotCosDist_NP()
{
  // read the histos from subtraction
  TFile *infile = new TFile("../NP_fit/files/bkgSubRes.root");
  const int n_inp = 3;
  TH2D **h_rat = new TH2D*[n_inp];
  string lbl[] = {"NP", "SB", "NPc"};
  for(int i = 0; i < n_inp; i++) {
    infile->GetObject(Form("h_%sB", lbl[i].c_str()), h_rat[i]);
    h_rat[i]->SetDirectory(0);
  }
  infile->Close();

  // get the binning
  int nBinsX = h_rat[0]->GetNbinsX(), nBinsY = h_rat[0]->GetNbinsY();
  const double *yBins = h_rat[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *h_rat1d[n_inp][nBinsY];
  for(int i_t = 0; i_t < n_inp; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      h_rat1d[i_t][i-1] = h_rat[i_t]->ProjectionX(Form("bin%d_%d_r", i, i_t+1), i, i);
    }
  }

  // the cycle to plot each bin
  TCanvas *c = new TCanvas("", "", 700, 700);    
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  int cols[] = {kRed+3, kGreen+1, kRed};
    
  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_rat[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_rat[0]->GetYaxis()->GetBinUpEdge(i+1);

    // draw all 3 fit results
    h_rat1d[0][i]->SetTitle("");
    h_rat1d[0][i]->SetStats(0);
    h_rat1d[0][i]->SetLineColor(cols[0]);
    h_rat1d[0][i]->SetMarkerColor(cols[0]);
    h_rat1d[0][i]->SetMinimum(0);
    h_rat1d[0][i]->SetMaximum(h_rat1d[0][i]->GetBinContent(1)*1.3);
    h_rat1d[0][i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    h_rat1d[0][i]->Draw("error");

    h_rat1d[1][i]->SetLineColor(cols[1]);
    h_rat1d[1][i]->SetMarkerColor(cols[1]);
    h_rat1d[1][i]->Draw("error same");
    
    h_rat1d[2][i]->SetLineColor(cols[2]);
    h_rat1d[2][i]->SetMarkerColor(cols[2]);
    h_rat1d[2][i]->Draw("error same");

    TLatex lcr1;
    lcr1.SetTextSize(0.04);
    lcr1.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.9, "Run 2");
    lcr1.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));
    lcr1.SetTextColor(cols[0]);
    lcr1.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.825, "NPS");
    lcr1.SetTextColor(cols[1]);
    lcr1.DrawLatex(0.15, h_rat1d[1][i]->GetMaximum()*1.1, "f_{B_{NP}}^{NPS}*B_{NP}");
    lcr1.SetTextColor(cols[2]);
    lcr1.DrawLatex(0.15, h_rat1d[2][i]->GetMaximum()*0.85, "non-prompt #psi(2S)");
    
    c->SaveAs(Form("plots/ratioFinal_NP/dists/bin1B_%d.pdf", i));
    c->Clear();
  }
  
  c->Destructor();
}
