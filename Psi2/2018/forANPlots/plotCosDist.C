// macro to plot costh dists for slides
void plotCosDist()
{
  // read the histos from subtraction - before normalization
  TFile *infile = new TFile("../PR_fit/files/histoStore.root");
  const int n_inp = 2;
  TH2D **h_base = new TH2D*[n_inp]; 
  string lbl[] = {"PR", "NP"};
  for(int i = 0; i < n_inp; i++) {
    infile->GetObject(Form("%sH", lbl[i].c_str()), h_base[i]);
    h_base[i]->SetDirectory(0);
  }
  // read the RATIO histos from subtraction
  TH2D **h_rat = new TH2D*[n_inp]; 
  for(int i = 0; i < n_inp; i++) {
    infile->GetObject(Form("r%sH", lbl[i].c_str()), h_rat[i]);
    h_rat[i]->SetDirectory(0);
  }
  infile->Close();

  // get the binning
  int nBinsX = h_base[0]->GetNbinsX(), nBinsY = h_base[0]->GetNbinsY();
  const double *yBins = h_base[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *h_base1d[n_inp][nBinsY];
  for(int i_t = 0; i_t < n_inp; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      h_base1d[i_t][i-1] = h_base[i_t]->ProjectionX(Form("bin%d_%d", i, i_t+1), i, i);
    }
  }
  TH1D *h_rat1d[n_inp][nBinsY];
  for(int i_t = 0; i_t < n_inp; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      h_rat1d[i_t][i-1] = h_rat[i_t]->ProjectionX(Form("binR%d_%d", i, i_t+1), i, i);
    }
  }
  
  // the cycle to plot each bin
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  int cols[] = {kViolet-1, kRed};
  
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
    lcb1.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.9, "2018");
    lcb1.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.85, Form("%.1f-%.1f GeV", pMin, pMax));
    lcb1.SetTextColor(cols[0]);
    lcb1.DrawLatex(0.15, h_base1d[0][i]->GetMaximum()*0.75, "Peak");
    
    c->SaveAs(Form("plots/ratioFinal/dists/bin1B_%d.pdf", i));

    // peak + NP costh
    h_base1d[1][i]->SetStats(0);
    h_base1d[1][i]->SetLineColor(cols[1]);
    h_base1d[1][i]->SetMarkerColor(cols[1]);
    h_base1d[1][i]->Draw("error same");

    lcb1.SetTextColor(cols[1]);
    lcb1.DrawLatex(0.15, h_base1d[1][i]->GetMaximum()*0.85, "NP");

    c->SaveAs(Form("plots/ratioFinal/dists/bin2B_%d.pdf", i));    
    c->Clear();

    // just peak costh - ratio
    h_rat1d[0][i]->SetTitle("");
    h_rat1d[0][i]->SetStats(0);
    h_rat1d[0][i]->SetLineColor(cols[0]);
    h_rat1d[0][i]->SetMarkerColor(cols[0]);
    h_rat1d[0][i]->SetMinimum(0);
    h_rat1d[0][i]->SetMaximum(h_rat1d[0][i]->GetBinContent(1)*1.5);
    h_rat1d[0][i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    h_rat1d[0][i]->Draw("error");
    
    TLatex lcr1;
    lcr1.SetTextSize(0.04);
    lcr1.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.9, "2018");
    lcr1.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.85, Form("%.1f-%.1f GeV", pMin, pMax));
    lcr1.SetTextColor(cols[0]);
    lcr1.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.75, "Peak/MC");
    
    c->SaveAs(Form("plots/ratioFinal/dists/bin1_%d.pdf", i));

    // peak + NP costh - ratio
    h_rat1d[1][i]->SetStats(0);
    h_rat1d[1][i]->SetLineColor(cols[1]);
    h_rat1d[1][i]->SetMarkerColor(cols[1]);
    h_rat1d[1][i]->Draw("error same");

    lcr1.SetTextColor(cols[1]);
    lcr1.DrawLatex(0.15, h_rat1d[1][i]->GetMaximum()*0.85, "NP/MC");

    c->SaveAs(Form("plots/ratioFinal/dists/bin2_%d.pdf", i));    
    c->Clear();

  }
  c->Destructor();
}
