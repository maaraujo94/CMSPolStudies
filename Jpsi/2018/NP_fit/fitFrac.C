// macro to fit the background fractions fNP and fSB
void fitFrac()
{
  // get background fractions
  TFile *fin1 = new TFile("files/mfit.root");
  TGraphErrors *fSB = (TGraphErrors*)fin1->Get("fit_fBG");
  fin1->Close();

  int nSB = fSB->GetN();
  double *yeSB = fSB->GetEY();
  for(int j = 0; j < nSB; j++) {
    yeSB[j] = fSB->GetY()[j]/100.;
  }
  
  // fit and plot f_SB
  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1F *fr1 = c->DrawFrame(20., 0.0, 125, 0.1);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f_{bkg}");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("2018 f_{bkg}");

  fSB->SetLineColor(kBlack);
  fSB->SetMarkerColor(kBlack);
  fSB->SetMarkerStyle(20);
  fSB->Draw("p");

  // fit function - parameters M(a,mu), a, mu
  TF1 *f_fit1 = new TF1("fit_SB", "[0]", 0, 125);
  f_fit1->SetParNames("fbkg");
  f_fit1->SetParameter(0, 0.07);
  f_fit1->SetLineColor(kBlue);
  fSB->Fit("fit_SB");
  f_fit1->SetRange(0, 125);
  f_fit1->Draw("same");

  c->SaveAs("plots/fSB_fit.pdf");
  c->Clear();

  // plot f_bkg pretty in %
  double *xvSB = fSB->GetX();
  double *xeSB = fSB->GetEX();
  double yvSB[nSB];
  // scale for fitting
  for(int j = 0; j < nSB; j++) {
    yvSB[j] = fSB->GetY()[j]*100;
    yeSB[j] = 0;
  }
  TGraphErrors *fSB_p = new TGraphErrors(nSB, xvSB, yvSB, xeSB, yeSB);

  TH1F *fp = c->DrawFrame(20, 0, 125, 15);
  fp->SetXTitle("p_{T} (GeV)");
  fp->SetYTitle("f_{bkg} (%)");
  fp->GetYaxis()->SetTitleOffset(1.3);
  fp->GetYaxis()->SetLabelOffset(0.01);
  fp->SetTitle(Form("2018 f_{bkg}"));
  
  fSB_p->SetMarkerStyle(20);
  fSB_p->SetMarkerSize(.75);
  fSB_p->SetLineColor(kBlue);
  fSB_p->SetMarkerColor(kBlue);
  fSB_p->Draw("p");	

  f_fit1->SetParameter(0, f_fit1->GetParameter(0)*100.);
  f_fit1->SetLineColor(kBlack);
  f_fit1->SetLineStyle(kDashed);
  f_fit1->Draw("lsame");

  c->SaveAs("plots/fSB_interp.pdf");
  c->Clear();
  c->Destructor();

  f_fit1->SetParameter(0, f_fit1->GetParameter(0)/100.);
  
  TFile *fout = new TFile("files/bkgFrac.root", "recreate");
  fSB->Write();
  f_fit1->Write();
  fout->Close();

}
