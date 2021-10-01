// macro to get the f_bkg uncertainty and compare to direct fit results
void fbkgUnc()
{
  // get the fractions obtained directly from the fits
  TFile *inSBo = new TFile("../../PR_fit/files/mfit.root");
  TGraphErrors *fit_fBG = (TGraphErrors*)inSBo->Get("fit_fBG");
  inSBo->Close();

  // fix f_SB to span 0 to 100
  int nbg = fit_fBG->GetN();
  double *fy = fit_fBG->GetY();
  double *fey = fit_fBG->GetEY();
  for(int i = 0; i < nbg; i++) {
    fy[i] *= 100.;
    fey[i] *= 100.;
  }

  // get the histos from the generation and obtain mean, std dev
  TH1F **h_fbg = new TH1F*[nbg];
  double bgy[nbg], bgey[nbg];
  TFile *inSB = new TFile("files/fbgDists.root");
  for(int i = 0; i < nbg; i++) {
    inSB->GetObject(Form("h_fbg_%d", i), h_fbg[i]);
    bgy[i] = h_fbg[i]->GetMean() * 100.;
    bgey[i] = h_fbg[i]->GetStdDev() * 100.;
  }
  inSB->Close();
  
  // define the new TGraphs
  TGraphErrors *gen_fBG = new TGraphErrors(nbg, fit_fBG->GetX(), bgy, fit_fBG->GetEX(), bgey);

  // plotting the comparisons
  TCanvas *c = new TCanvas("", "", 900, 900);

  // f_bg first
  double x_min = fit_fBG->GetX()[0]-fit_fBG->GetEX()[0]-5;
  double x_max = fit_fBG->GetX()[nbg-1]+fit_fBG->GetEX()[nbg-1]+5;
  TH1F *ffbg = c->DrawFrame(x_min, 0, x_max, 15);
  ffbg->SetXTitle("p_{T} (GeV)");
  ffbg->SetYTitle("f_{BG} (%)");
  ffbg->GetYaxis()->SetTitleOffset(1.3);
  ffbg->GetYaxis()->SetLabelOffset(0.01);
  ffbg->SetTitle("2018 f_{BG} vs p_{T}");
  
  fit_fBG->SetMarkerStyle(20);
  fit_fBG->SetMarkerSize(.75);
  fit_fBG->SetMarkerColor(kBlack);
  fit_fBG->SetLineColor(kBlack);
  fit_fBG->Draw("psame");

  gen_fBG->SetMarkerStyle(24);
  gen_fBG->SetMarkerSize(.75);
  gen_fBG->SetMarkerColor(kBlue);
  gen_fBG->SetLineColor(kBlue);
  gen_fBG->Draw("psame");

  c->SaveAs(Form("fBG_comp.pdf"));
  c->Clear();

  TFile *fout = new TFile("files/fbkgUnc.root", "recreate");
  gen_fBG->SetName("f_bkg");
  gen_fBG->Write();
  fout->Close();
}
