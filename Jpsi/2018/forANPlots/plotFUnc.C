void plotFUnc()
{
  // get the fractions obtained directly from the fits
  TFile *inSBo = new TFile("../PR_fit/files/mfit.root");
  TGraphErrors *fit_fBGf = (TGraphErrors*)inSBo->Get("fit_fBG");
  inSBo->Close();
  TFile *inNPo = new TFile("../PR_fit/files/ltfit.root");
  TGraphErrors *fit_fNP = (TGraphErrors*)inNPo->Get("fit_b_fNP");
  inNPo->Close();

  // fix f_SB to span 0 to 100
  int nbg = fit_fBGf->GetN(), nnp = fit_fNP->GetN();
  double *fx = fit_fBGf->GetX();
  double *fy = fit_fBGf->GetY();
  double *fex = fit_fBGf->GetEX();
  double *fey = fit_fBGf->GetEY();
  for(int i = 0; i < nbg; i++) {
    fy[i] *= 100.;
    fey[i] *= 100.;
  }
  TGraphErrors *fit_fBG = new TGraphErrors(nbg, fx, fy, fex, fey);

  // get the histos from the generation and obtain mean, std dev
  TH1F **h_fbg = new TH1F*[nbg];
  double bgy[nbg], bgey[nbg];
  TFile *inSB = new TFile("fbg_gen.root");
  for(int i = 0; i < nbg; i++) {
    inSB->GetObject(Form("h_fbg_%d", i), h_fbg[i]);
    bgy[i] = h_fbg[i]->GetMean() * 100.;
    bgey[i] = h_fbg[i]->GetStdDev() * 100.;
  }
  inSB->Close();
  
  TH1F **h_fnp = new TH1F*[nnp];
  double npy[nnp], npey[nnp];
  TFile *inNP = new TFile("fnp_gen.root");
  for(int i = 0; i < nnp; i++) {
    inNP->GetObject(Form("h_fnp_%d", i), h_fnp[i]);
    npy[i] = h_fnp[i]->GetMean() * 100.;
    npey[i] = h_fnp[i]->GetStdDev() * 100.;
  }
  inNP->Close();

  // define the new TGraphs
  TGraphErrors *gen_fBG = new TGraphErrors(nbg, fx, bgy, fex, bgey);
  double *fnx = fit_fNP->GetX();
  double *fnex = fit_fNP->GetEX();
  TGraphErrors *gen_fNP = new TGraphErrors(nnp, fnx, npy, fnex, npey);

  // plotting the comparisons
  TCanvas *c = new TCanvas("", "", 900, 900);

  // f_bg first
  TH1F *ffbg = c->DrawFrame(fx[0]-fex[0]-5, 0, fx[nbg-1]+fex[nbg-1]+5, 15);
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

  // f_np second
  TH1F *ffnp = c->DrawFrame(fnx[0]-fnex[0]-5, 0, fnx[nnp-1]+fnex[nnp-1]+5, 50);
  ffnp->SetXTitle("p_{T} (GeV)");
  ffnp->SetYTitle("f_{NP} (%)");
  ffnp->GetYaxis()->SetTitleOffset(1.3);
  ffnp->GetYaxis()->SetLabelOffset(0.01);
  ffnp->SetTitle("2018 f_{NP} vs p_{T}");
  
  fit_fNP->SetMarkerStyle(20);
  fit_fNP->SetMarkerSize(.75);
  fit_fNP->SetMarkerColor(kBlack);
  fit_fNP->SetLineColor(kBlack);
  fit_fNP->Draw("psame");

  gen_fNP->SetMarkerStyle(24);
  gen_fNP->SetMarkerSize(.75);
  gen_fNP->SetMarkerColor(kBlue);
  gen_fNP->SetLineColor(kBlue);
  gen_fNP->Draw("psame");

  c->SaveAs(Form("fNP_comp.pdf"));
  c->Clear();

  
}
