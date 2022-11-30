#import "../ptbins.C"
void plotfNP()
{
  // get fNP from 1d fit, with all pars free (just points)
  TFile *fin1 = new TFile("../PR_fit/files/ltfit.root");
  TGraphErrors *g_f1d = (TGraphErrors*)fin1->Get("fit_fNP");
  fin1->Close();
  // get fNP from 2d fit, with uncertainties
  TFile *fin2 = new TFile("../PR_fit/files/NPFrac.root");
  TH1D *h_f2d = (TH1D*)fin2->Get("fnp_unc");
  h_f2d->SetDirectory(0);
  fin2->Close();
  
  // prep the plotting
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);

  TH1F *fp = c->DrawFrame(ptBins[0]-5, 0, ptBins[nPtBins], 50);
  fp->SetXTitle("p_{T} (GeV)");
  fp->SetYTitle("f_{NP} (%)");
  fp->GetYaxis()->SetTitleOffset(1.4);
  fp->GetYaxis()->SetLabelOffset(0.01);
  fp->SetTitle(Form("f_{NP} vs p_{T}"));

  h_f2d->Scale(100);
  h_f2d->SetMarkerColor(kBlack);
  h_f2d->SetLineColor(kBlack);
  h_f2d->Draw("error same");

  g_f1d->SetMarkerColor(kRed);
  g_f1d->SetLineColor(kRed);
  g_f1d->Draw("psame");

  c->SaveAs("plots/fNP.pdf");
  c->Destructor();
}
