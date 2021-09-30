// macro to fit the background fractions fNP and fSB
void fitFrac()
{
  // get background fractions
  TFile *fin1 = new TFile("files/mfit.root");
  TGraphErrors *fSB = (TGraphErrors*)fin1->Get("fit_fBG");
  fin1->Close();

  // fit and plot f_SB
  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1F *fr1 = c->DrawFrame(20., 0.0, 125, 0.08);
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
  TF1 *f_fit1 = new TF1("fit_SB", "[2]*0.020*(1-exp(-[0]*(x-[1])))/(1-exp(-[0]*(20-[1])))", 0, 125);
  f_fit1->SetParNames("a", "mu", "ph");
  f_fit1->SetParameters(0.02, -10, 1);
  f_fit1->FixParameter(2, 1);
  f_fit1->SetLineColor(kBlue);
  fSB->Fit("fit_SB");
  f_fit1->SetRange(0, 125);
  f_fit1->Draw("same");

  c->SaveAs("plots/fSB_fit.pdf");
  c->Clear();

  // plot f_bkg pretty in %
  int nSB = fSB->GetN();
  double *xvSB = fSB->GetX();
  double *xeSB = fSB->GetEX();
  double yvSB[nSB], yeSB[nSB];
  // scale for fitting
  for(int j = 0; j < nSB; j++) {
    yvSB[j] = fSB->GetY()[j]*100;
    yeSB[j] = 0;
  }
  TGraphErrors *fSB_p = new TGraphErrors(nSB, xvSB, yvSB, xeSB, yeSB);

  TH1F *fp = c->DrawFrame(20, 0, 125, 15);
  fp->SetXTitle("p_{T} (GeV)");
  fp->SetYTitle("f_{bkg} (%)");
  fp->GetYaxis()->SetTitleOffset(1.5);
  fp->GetYaxis()->SetLabelOffset(0.01);
  fp->SetTitle(Form("2018 f_{bkg}"));
  
  fSB_p->SetMarkerStyle(20);
  fSB_p->SetMarkerSize(.75);
  fSB_p->SetLineColor(kBlue);
  fSB_p->SetMarkerColor(kBlue);
  fSB_p->Draw("p");	

  f_fit1->SetParameter(2, 100);
  f_fit1->SetLineColor(kBlack);
  f_fit1->SetLineStyle(kDashed);
  f_fit1->Draw("lsame");

  c->SaveAs("plots/fSB_interp.pdf");
  c->Clear();
  c->Destructor();

  f_fit1->SetParameter(2,1);
  
  TFile *fout = new TFile("files/bkgFrac.root", "recreate");
  fSB->Write();
  f_fit1->Write();
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
