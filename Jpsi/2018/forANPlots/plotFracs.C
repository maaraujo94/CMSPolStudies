// macro to plot the bkg fracs and prompt J/psi frac over pT
void plotFracs()
{
  // get fit fBG and fNP
  TFile *fin1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/2018/PR_fit/files/mfit.root");
  TGraphErrors *fSB = (TGraphErrors*)fin1->Get("fit_fBG");
  fin1->Close();

  TFile *fin2 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/2018/PR_fit/files/ltfit.root");
  TGraphErrors *fNP_p = (TGraphErrors*)fin2->Get("fit_b_fNP");
  fin2->Close();

  // get the plateau fit for fBG
  TF1 *f_fit1 = new TF1("f_fit1", "0.029*(1-exp(-[0]*(x-[1])))/(1-exp(-[0]*(20-[1])))", 0, 125);
  f_fit1->SetParameters(0.02, -10);
  f_fit1->SetLineColor(kBlack);
  f_fit1->SetLineStyle(kDashed);
  fSB->Fit("f_fit1", "N");
  f_fit1->SetRange(25, 120);

  // define NP (it's as % currently) and prompt J/psi
  int n = fNP_p->GetN();
  double *x = fNP_p->GetX();
  double *ex = fNP_p->GetEX();
  double *ey = fNP_p->GetEY();
  double yJ[n], yN[n];
  for(int i = 0; i < n; i++) {
    yN[i] = fNP_p->GetY()[i]/100;
    double f_sb = f_fit1->Eval(x[i]);
    yJ[i] = 1-yN[i]-f_sb;
  }
  TGraphErrors *fPRJ = new TGraphErrors(n, x, yJ, ex, ey);
  TGraphErrors *fNP = new TGraphErrors(n, x, yN, ex, ey);

  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1F *fr1 = c->DrawFrame(20, 0.0, 125, 1.0);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("f comparison");

  fNP->SetLineColor(kRed);
  fNP->SetMarkerColor(kRed);
  fNP->SetMarkerStyle(20);
  fNP->Draw("p");

  fSB->SetLineColor(kGreen);
  fSB->SetMarkerColor(kGreen);
  fSB->SetMarkerStyle(20);
  fSB->Draw("p");

  f_fit1->Draw("same");
  
  fPRJ->SetLineColor(kBlue);
  fPRJ->SetMarkerColor(kBlue);
  fPRJ->SetMarkerStyle(20);
  fPRJ->Draw("p");

  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(fPRJ, "prompt J/#psi", "pl");
  leg->AddEntry(fNP, "NP", "pl");
  leg->AddEntry(fSB, "SB", "pl");
  leg->Draw();

  c->SaveAs("plots/f_comp.pdf");
  c->Clear();
  c->Destructor();
}
