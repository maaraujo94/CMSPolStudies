// macro to plot costh dists for slides
void plotCosDist()
{
  // read the histos from subtraction - normalized by f_NP/f_bkg
  TFile *infile = new TFile("../PR_fit/files/bkgSubRes.root");
  TH2D **h_base = new TH2D*[5]; 
  string lbl[] = {"Data", "NP", "PR", "J", "SB"};
  for(int i = 0; i < 5; i++) {
    infile->GetObject(Form("h_%sB", lbl[i].c_str()), h_base[i]);
    h_base[i]->SetDirectory(0);
  }
  infile->Close();

  // get the binning
  int nBinsX = h_base[0]->GetNbinsX(), nBinsY = h_base[0]->GetNbinsY();
  const double *yBins = h_base[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *h_base1d[5][nBinsY];
  for(int i_t = 0; i_t < 5; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      h_base1d[i_t][i-1] = h_base[i_t]->ProjectionX(Form("bin%d_%d", i, i_t+1), i, i);
    }
  }
  
  // the cycle to plot each bin
  TCanvas *c = new TCanvas("", "", 700, 700);    
  int cols[] = {kViolet-1, kRed, kBlack, kBlue, kGreen};
  
  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_base[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_base[0]->GetYaxis()->GetBinUpEdge(i+1);

    // just peak costh
    h_base1d[0][i]->SetTitle("");
    h_base1d[0][i]->SetStats(0);
    h_base1d[0][i]->SetLineColor(cols[0]);
    h_base1d[0][i]->SetMarkerColor(cols[0]);
    h_base1d[0][i]->SetMinimum(0);
    h_base1d[0][i]->SetMaximum(h_base1d[0][i]->GetBinContent(1)*1.5);
    h_base1d[0][i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    h_base1d[0][i]->Draw("error");
    
    TLatex lcb1;
    lcb1.SetTextSize(0.04);
    lcb1.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.9, "2017");
    lcb1.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));
    lcb1.SetTextColor(cols[0]);
    lcb1.DrawLatex(0.15, h_base1d[0][i]->GetMaximum()*0.8, "Peak");
    
    c->SaveAs(Form("plots/ratioFinal/dists/bin1B_%d.pdf", i));

    // peak + NP (norm) + (peak-NP) costh
    h_base1d[1][i]->SetStats(0);
    h_base1d[1][i]->SetLineColor(cols[1]);
    h_base1d[1][i]->SetMarkerColor(cols[1]);
    h_base1d[1][i]->Draw("error same");

    lcb1.SetTextColor(cols[1]);
    lcb1.DrawLatex(0.15, h_base1d[1][i]->GetMaximum()*1.1, "NP (scaled)");

    h_base1d[2][i]->SetStats(0);
    h_base1d[2][i]->SetLineColor(cols[2]);
    h_base1d[2][i]->SetMarkerColor(cols[2]);
    h_base1d[2][i]->Draw("error same");

    lcb1.SetTextColor(cols[2]);
    lcb1.DrawLatex(0.15, h_base1d[2][i]->GetMaximum()*0.7, "PR");

    c->SaveAs(Form("plots/ratioFinal/dists/bin2B_%d.pdf", i));    
    c->Clear();

    // peak + NP (norm) + (peak-NP) + SB (norm) + pr J/psi costh
    h_base1d[0][i]->Draw("error");
    
    TLatex lcb3;
    lcb3.SetTextSize(0.04);
    lcb3.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.9, "2017");
    lcb3.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));
    lcb3.SetTextColor(cols[0]);
    lcb3.DrawLatex(0.15, h_base1d[0][i]->GetMaximum()*0.8, "Peak");

    h_base1d[1][i]->Draw("error same");

    lcb3.SetTextColor(cols[1]);
    lcb3.DrawLatex(0.15, h_base1d[1][i]->GetMaximum()*1.1, "NP (scaled)");

    h_base1d[2][i]->Draw("error same");

    lcb3.SetTextColor(cols[2]);
    lcb3.DrawLatex(0.15, h_base1d[2][i]->GetMaximum()*1.0, "PR");

    h_base1d[3][i]->SetStats(0);
    h_base1d[3][i]->SetLineColor(cols[3]);
    h_base1d[3][i]->SetMarkerColor(cols[3]);
    h_base1d[3][i]->Draw("error same");

    lcb3.SetTextColor(cols[3]);
    lcb3.DrawLatex(0.15, h_base1d[3][i]->GetMaximum()*0.7, "J/#psi");

    h_base1d[4][i]->SetStats(0);
    h_base1d[4][i]->SetLineColor(cols[4]);
    h_base1d[4][i]->SetMarkerColor(cols[4]);
    h_base1d[4][i]->Draw("error same");

    lcb3.SetTextColor(cols[4]);
    lcb3.DrawLatex(0.15, h_base1d[4][i]->GetMaximum()*1.1, "bkg (scaled)");

    c->SaveAs(Form("plots/ratioFinal/dists/bin3B_%d.pdf", i));    
    c->Clear();
  }
  c->Destructor();
}
