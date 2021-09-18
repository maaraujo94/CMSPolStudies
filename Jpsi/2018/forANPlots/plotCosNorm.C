// macro to plot costh dists for slides
void plotCosNorm()
{
  // read the histos from subtraction - normalized by f_NP/f_bkg
  TFile *infile = new TFile("../PR_fit/files/bkgSubRes.root");
  TH2D **h_rat = new TH2D*[5];
  TH2D **h_base = new TH2D*[5]; 
  string lbl[] = {"Data", "NP", "PR", "J", "SB"};
  for(int i = 0; i < 5; i++) {
    infile->GetObject(Form("h_%s", lbl[i].c_str()), h_rat[i]);
    infile->GetObject(Form("h_%sB", lbl[i].c_str()), h_base[i]);
    h_rat[i]->SetDirectory(0);
    h_base[i]->SetDirectory(0);
  }
  infile->Close();

  // read the base histos - non-normalized
  TFile *infile2 = new TFile("../PR_fit/files/histoStore.root");
  TH2D **h_nnr = new TH2D*[2];
  TH2D **h_nnb = new TH2D*[2];
  string lbl2[] = {"data", "NP"};
  for(int i = 0; i < 2; i++) {
    infile2->GetObject(Form("%sH_ab", lbl2[i].c_str()), h_nnb[i]);
    h_nnb[i]->SetDirectory(0);
  }
  infile2->GetObject("ratioH_ab", h_nnr[0]);
  h_nnr[0]->SetDirectory(0);
  infile2->GetObject("ratNPH_ab", h_nnr[1]);
  h_nnr[1]->SetDirectory(0);
  infile2->Close();

  // get the binning
  int nBinsX = h_rat[0]->GetNbinsX(), nBinsY = h_rat[0]->GetNbinsY();
  const double *yBins = h_rat[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *h_rat1d[5][nBinsY];
  TH1D *h_base1d[5][nBinsY];
  TH1D *h_nnr1d[2][nBinsY];
  TH1D *h_nnb1d[2][nBinsY];
  for(int i_t = 0; i_t < 5; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      h_rat1d[i_t][i-1] = h_rat[i_t]->ProjectionX(Form("bin%d_%d_r", i, i_t+1), i, i);
      h_base1d[i_t][i-1] = h_base[i_t]->ProjectionX(Form("bin%d_%d_b", i, i_t+1), i, i);
    }
  }
  for(int i_t = 0; i_t < 2; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      h_nnr1d[i_t][i-1] = h_nnr[i_t]->ProjectionX(Form("bin%d_%d_nr", i, i_t+1), i, i);
      h_nnb1d[i_t][i-1] = h_nnb[i_t]->ProjectionX(Form("bin%d_%d_nb", i, i_t+1), i, i);
    }
  }
  
 
  // the cycle to plot each bin
  TCanvas *c = new TCanvas("", "", 700, 700);    
  int cols[] = {kViolet-1, kRed, kBlack, kBlue, kGreen};
  
  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_rat[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_rat[0]->GetYaxis()->GetBinUpEdge(i+1);

    // PART 1: plot the distributions before /MC
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
    lcb1.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));
    lcb1.SetTextColor(cols[0]);
    lcb1.DrawLatex(0.15, h_base1d[0][i]->GetMaximum()*0.8, "Peak");
    
    c->SaveAs(Form("plots/ratioFinal/bin1B_%d.pdf", i));
    c->Clear();

    // peak + NP (non-norm) costh
    h_base1d[0][i]->Draw("error");
    
    TLatex lcb2;
    lcb2.SetTextSize(0.04);
    lcb2.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.9, "2018");
    lcb2.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));
    lcb2.SetTextColor(cols[0]);
    lcb2.DrawLatex(0.15, h_base1d[0][i]->GetMaximum()*0.8, "Peak");
    
    h_nnb1d[1][i]->SetStats(0);
    h_nnb1d[1][i]->SetLineColor(cols[1]);
    h_nnb1d[1][i]->SetMarkerColor(cols[1]);
    h_nnb1d[1][i]->Draw("error same");

    lcb2.SetTextColor(cols[1]);
    lcb2.DrawLatex(0.15, h_nnb1d[1][i]->GetMaximum()*0.8, "NP");

    c->SaveAs(Form("plots/ratioFinal/bin2B_%d.pdf", i));    
    c->Clear();

    // peak + NP (norm) + (peak-NP) costh
    h_base1d[0][i]->Draw("error");
    
    TLatex lcb3;
    lcb3.SetTextSize(0.04);
    lcb3.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.9, "2018");
    lcb3.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));
    lcb3.SetTextColor(cols[0]);
    lcb3.DrawLatex(0.15, h_base1d[0][i]->GetMaximum()*0.8, "Peak");

    h_base1d[1][i]->SetStats(0);
    h_base1d[1][i]->SetLineColor(cols[1]);
    h_base1d[1][i]->SetMarkerColor(cols[1]);
    h_base1d[1][i]->Draw("error same");

    lcb3.SetTextColor(cols[1]);
    lcb3.DrawLatex(0.15, h_base1d[1][i]->GetMaximum()*1.1, "NP (scaled)");

    h_base1d[2][i]->SetStats(0);
    h_base1d[2][i]->SetLineColor(cols[2]);
    h_base1d[2][i]->SetMarkerColor(cols[2]);
    h_base1d[2][i]->Draw("error same");

    lcb3.SetTextColor(cols[2]);
    lcb3.DrawLatex(0.15, h_base1d[2][i]->GetMaximum()*0.7, "PR");

    c->SaveAs(Form("plots/ratioFinal/bin3B_%d.pdf", i));    
    c->Clear();

    // peak + NP (norm) + (peak-NP) + SB (norm) + pr J/psi costh
    h_base1d[0][i]->Draw("error");
    
    TLatex lcb4;
    lcb4.SetTextSize(0.04);
    lcb4.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.9, "2018");
    lcb4.DrawLatex(0.7, h_base1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));
    lcb4.SetTextColor(cols[0]);
    lcb4.DrawLatex(0.15, h_base1d[0][i]->GetMaximum()*0.8, "Peak");

    h_base1d[1][i]->Draw("error same");

    lcb4.SetTextColor(cols[1]);
    lcb4.DrawLatex(0.15, h_base1d[1][i]->GetMaximum()*1.1, "NP (scaled)");

    h_base1d[2][i]->Draw("error same");

    lcb4.SetTextColor(cols[2]);
    lcb4.DrawLatex(0.15, h_base1d[2][i]->GetMaximum()*1.0, "PR");

    h_base1d[3][i]->SetStats(0);
    h_base1d[3][i]->SetLineColor(cols[3]);
    h_base1d[3][i]->SetMarkerColor(cols[3]);
    h_base1d[3][i]->Draw("error same");

    lcb4.SetTextColor(cols[3]);
    lcb4.DrawLatex(0.15, h_base1d[3][i]->GetMaximum()*0.7, "J/#psi");

    h_base1d[4][i]->SetStats(0);
    h_base1d[4][i]->SetLineColor(cols[4]);
    h_base1d[4][i]->SetMarkerColor(cols[4]);
    h_base1d[4][i]->Draw("error same");

    lcb4.SetTextColor(cols[4]);
    lcb4.DrawLatex(0.15, h_base1d[4][i]->GetMaximum()*1.1, "bkg (scaled)");

    c->SaveAs(Form("plots/ratioFinal/bin4B_%d.pdf", i));    
    c->Clear();

    // PART 2: plotting ratio distributions
    // just peak/MC costh
    h_rat1d[0][i]->SetTitle("");
    h_rat1d[0][i]->SetStats(0);
    h_rat1d[0][i]->SetLineColor(kBlack);
    h_rat1d[0][i]->SetMarkerColor(kBlack);
    h_rat1d[0][i]->SetMinimum(0);
    h_rat1d[0][i]->SetMaximum(h_rat1d[0][i]->GetBinContent(1)*1.5);
    h_rat1d[0][i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    h_rat1d[0][i]->Draw("error");
    
    TLatex lcr1;
    lcr1.SetTextSize(0.04);
    lcr1.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.8, "Peak/MC");
    lcr1.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.9, "2018");
    lcr1.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));
    
    c->SaveAs(Form("plots/ratioFinal/bin1_%d.pdf", i));
    c->Clear();

    // peak/MC + NP/MC (non-norm)
    h_rat1d[0][i]->Draw("error");
    
    TLatex lcr2;
    lcr2.SetTextSize(0.04);
    lcr2.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.8, "Peak/MC");
    lcr2.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.9, "2018");
    lcr2.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));

    h_nnr1d[1][i]->SetStats(0);
    h_nnr1d[1][i]->SetLineColor(kBlue);
    h_nnr1d[1][i]->SetMarkerColor(kBlue);
    h_nnr1d[1][i]->Draw("error same");

    lcr2.SetTextColor(kBlue);
    lcr2.DrawLatex(0.15, h_nnr1d[1][i]->GetMaximum()*0.8, "NP/MC");

    c->SaveAs(Form("plots/ratioFinal/bin2_%d.pdf", i));
    c->Clear();

    // peak/MC + NP/MC (norm) + (peak-NP)/MC 
    h_rat1d[0][i]->Draw("error");
    
    TLatex lcr3;
    lcr3.SetTextSize(0.04);
    lcr3.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.8, "Peak/MC");
    lcr3.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.9, "2018");
    lcr3.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));

    h_rat1d[1][i]->SetStats(0);
    h_rat1d[1][i]->SetLineColor(kBlue);
    h_rat1d[1][i]->SetMarkerColor(kBlue);
    h_rat1d[1][i]->Draw("error same");

    lcr3.SetTextColor(kBlue);
    lcr3.DrawLatex(0.15, h_rat1d[1][i]->GetMaximum()*1.1, "NP/MC (scaled)");
    
    h_rat1d[2][i]->SetStats(0);
    h_rat1d[2][i]->SetLineColor(kViolet);
    h_rat1d[2][i]->SetMarkerColor(kViolet);
    h_rat1d[2][i]->Draw("error same");

    lcr3.SetTextColor(kViolet);
    lcr3.DrawLatex(0.15, h_rat1d[2][i]->GetMaximum()*0.7, "PR/MC");

    c->SaveAs(Form("plots/ratioFinal/bin3_%d.pdf", i));
    c->Clear();

    // peak/MC + NP/MC (norm) + (peak-NP)/MC + SB/MC (norm) + pr J/psi/MC costh
    h_rat1d[0][i]->Draw("error");
    
    TLatex lcr4;
    lcr4.SetTextSize(0.04);
    lcr4.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.8, "Peak/MC");
    lcr4.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.9, "2018");
    lcr4.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));

    h_rat1d[1][i]->Draw("error same");

    lcr4.SetTextColor(kBlue);
    lcr4.DrawLatex(0.15, h_rat1d[1][i]->GetMaximum()*1.1, "NP/MC (scaled)");

    h_rat1d[2][i]->Draw("error same");

    lcr4.SetTextColor(kViolet);
    lcr4.DrawLatex(0.15, h_rat1d[2][i]->GetMaximum()*0.95, "PR/MC");

    h_rat1d[3][i]->SetStats(0);
    h_rat1d[3][i]->SetLineColor(kRed);
    h_rat1d[3][i]->SetMarkerColor(kRed);
    h_rat1d[3][i]->Draw("error same");

    lcr4.SetTextColor(kRed);
    lcr4.DrawLatex(0.15, h_rat1d[3][i]->GetMaximum()*0.7, "J/#psi/MC");

    h_rat1d[4][i]->SetStats(0);
    h_rat1d[4][i]->SetLineColor(kGreen);
    h_rat1d[4][i]->SetMarkerColor(kGreen);
    h_rat1d[4][i]->Draw("error same");

    lcr4.SetTextColor(kGreen);
    lcr4.DrawLatex(0.15, h_rat1d[4][i]->GetMaximum()*1.1, "bkg (scaled)");

    c->SaveAs(Form("plots/ratioFinal/bin4_%d.pdf", i));    
    c->Clear();


  }
  
  c->Destructor();
}
