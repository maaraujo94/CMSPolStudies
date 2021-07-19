// macro to fit the background fractions fNP and fSB
void fitFrac()
{
  // get background fractions
  TFile *fin1 = new TFile("files/mfit.root");
  TGraphErrors *fSB = (TGraphErrors*)fin1->Get("fit_fBG");
  fin1->Close();

  // NP fraction is in %, need to fix that
  TFile *fin2 = new TFile("files/ltfit.root");
  TGraphErrors *fNP_p = (TGraphErrors*)fin2->Get("fit_b_fNP");
  fin2->Close();

  int n = fNP_p->GetN();
  double *xv = fNP_p->GetX();
  double *xe = fNP_p->GetEX();
  double *yv = fNP_p->GetY();
  double *ye = fNP_p->GetEY();
  // scale for fitting
  for(int j = 0; j < n; j++) {
    yv[j] /= 100.;
    ye[j] = xv[j]*1e-4;
  }
  TGraphErrors *fNP = new TGraphErrors(n, xv, yv, xe, ye);

  // fit and plot f_SB
  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1F *fr1 = c->DrawFrame(20., 0.0, 125, 0.08);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f_{bkg}");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("2017 f_{bkg}");

  fSB->SetLineColor(kBlack);
  fSB->SetMarkerColor(kBlack);
  fSB->SetMarkerStyle(20);
  fSB->Draw("p");

  // fit function - parameters M(a,mu), a, mu
  TF1 *f_fit1 = new TF1("fit_SB", "0.022*(1-exp(-[0]*(x-[1])))/(1-exp(-[0]*(20-[1])))", 0, 125);
  f_fit1->SetParNames("a", "mu");
  f_fit1->SetParameters(0.02, -10);
  f_fit1->SetLineColor(kBlue);
  fSB->Fit("fit_SB");
  f_fit1->SetRange(0, 125);
  f_fit1->Draw("same");

  c->SaveAs("plots/fSB_fit.pdf");
  c->Clear();

  // fit and plot f_NP
  TH1F *fr2 = c->DrawFrame(0., 0., 125, 0.4);
  fr2->SetXTitle("p_{T} (GeV)");
  fr2->SetYTitle("f_{NP}");
  fr2->GetYaxis()->SetTitleOffset(1.3);
  fr2->GetYaxis()->SetLabelOffset(0.01);
  fr2->SetTitle("f_{NP}");

  fNP->SetLineColor(kBlack);
  fNP->SetMarkerColor(kBlack);
  fNP->SetMarkerStyle(20);
  fNP->Draw("p");

  // fit function - parameters M, a (mu = 0)
  TF1 *f_fit2 = new TF1("fit_NP", "[0]*(1-exp(-[1]*(x-[2])))", 0, 125);
  f_fit2->SetParNames("M", "a", "mu");
  f_fit2->SetParameters(0.3, 0.05, 10);
  f_fit2->FixParameter(2,0);
  f_fit2->SetLineColor(kBlue);
  fNP->Fit("fit_NP");
  f_fit2->SetRange(0, 125);
  f_fit2->Draw("same");

  c->SaveAs("plots/fNP_fit.pdf");
  c->Clear();
  c->Destructor();

  TFile *fout = new TFile("files/bkgFrac.root", "recreate");
  fSB->Write();
  f_fit1->Write();
  fNP->Write();
  f_fit2->Write();
  fout->Close();

  ofstream fout_t;
  fout_t.open("text_output/fit_frac.tex");
  fout_t << "\\begin{tabular}{c||c|c|c}\n";
  fout_t << " & $M$ & $a$ & $\\mu$  \\\\\n";
  fout_t << "\\hline\n";
  fout_t << "$f_{bkg}$ ";
  double val = 0.022/(1-exp(-f_fit1->GetParameter(0)*(20-f_fit1->GetParameter(1))));
  int p_norm = 1;
  if(val > 0 && val < 1 ) 
    p_norm = ceil(-log10(val))+1;
  fout_t << " & " <<  setprecision(p_norm) << fixed << val;
  for(int i = 0; i < 2; i++) {
    val = f_fit1->GetParameter(i);
    p_norm = 1.;
    if(val > 0 && val < 1 ) 
      p_norm = ceil(-log10(val))+1;
    fout_t << " & " <<  setprecision(p_norm) << fixed << val;
  }
  fout_t << "\\\\\n";
  fout_t << "\\end{tabular}\n";
  fout_t.close();
}
