void plotFracs()
{
  TFile *fin1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/2018/PR_fit/files/mfit.root");
  TGraphErrors *fSB = (TGraphErrors*)fin1->Get("fit_fBG");
  fin1->Close();

  TFile *fin2 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/2018/PR_fit/files/tfit.root");
  TGraphErrors *fNP = (TGraphErrors*)fin2->Get("fit_fNP");
  fin2->Close();

  int n = fSB->GetN();
  double *x = fSB->GetX();
  double *ex = fSB->GetEX();
  double *ey = fSB->GetEY();
  double y[n];
  for(int i = 0; i < n; i++) {
    y[i] = 1-fSB->GetY()[i]-fNP->GetY()[i];
  }
  TGraphErrors *fPRJ = new TGraphErrors(n, x, y, ex, ey);

  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1F *fr1 = c->DrawFrame(20, 0.0, 125, 1.0);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("f comparison");

  TF1 *f_fit2 = new TF1("f_fit2", "[0]*(1-exp(-[1]*(x-[2])))", 0, 125);
  f_fit2->SetParameters(0.4, 0.05, 0);
  f_fit2->FixParameter(2,0);
  f_fit2->SetLineColor(kBlack);
  f_fit2->SetLineStyle(kDashed);
  fNP->Fit("f_fit2", "N");
  f_fit2->SetRange(25, 120);
  f_fit2->Draw("same");

  fNP->SetLineColor(kRed);
  fNP->SetMarkerColor(kRed);
  fNP->SetMarkerStyle(20);
  fNP->Draw("p");

  TF1 *f_fit1 = new TF1("f_fit1", "0.029*(1-exp(-[0]*(x-[1])))/(1-exp(-[0]*(20-[1])))", 0, 125);
  f_fit1->SetParameters(0.02, -10);
  f_fit1->SetLineColor(kBlack);
  f_fit1->SetLineStyle(kDashed);
  fSB->Fit("f_fit1", "N");
  f_fit1->SetRange(25, 120);
  f_fit1->Draw("same");

  fSB->SetLineColor(kGreen);
  fSB->SetMarkerColor(kGreen);
  fSB->SetMarkerStyle(20);
  fSB->Draw("p");

  TF1 *f_sum = new TF1("f_sum", "1-0.029*(1-exp(-[0]*(x-[1])))/(1-exp(-[0]*(20-[1])))-[2]*(1-exp(-[3]*(x-[4])))", 25, 120);
  f_sum->SetParameters(f_fit1->GetParameter(0), f_fit1->GetParameter(1), f_fit2->GetParameter(0), f_fit2->GetParameter(1), f_fit2->GetParameter(2));
  f_sum->SetLineColor(kBlack);
  f_sum->SetLineStyle(kDashed);
  f_sum->SetRange(25, 120);
  f_sum->Draw("same");

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
