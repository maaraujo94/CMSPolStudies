// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

void plotFrac()
{
  TFile *finB = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/bkgFits/files/bkgFrac_NP.root");
  TH1D* h_fB = (TH1D*)finB->Get("fbkg_unc");
  h_fB->SetDirectory(0);
  finB->Close();

  TFile *finN = new TFile("files/bkgFrac.root");
  TH1D* h_fN = (TH1D*)finN->Get("fbkg_unc");
  h_fN->SetDirectory(0);
  finN->Close();

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);

  // plot all fracs (bkg corrected NP)
  double ptmin = 20, ptmax = 125;
  
  TH1F *fr2 = c->DrawFrame(ptmin, 0.0, ptmax, 15);
  fr2->SetXTitle("#it{p}_{T} (GeV)");
  fr2->SetYTitle("fraction");
  fr2->GetYaxis()->SetTitleOffset(1.4);
  fr2->GetYaxis()->SetLabelOffset(0.01);
  fr2->GetXaxis()->SetLabelOffset(0.015);
  fr2->GetXaxis()->SetTitleOffset(1.3);

  h_fN->Scale(100.);
  h_fN->SetLineColor(kRed);
  h_fN->SetMarkerColor(kRed);
  h_fN->SetMarkerStyle(21);
  h_fN->SetMarkerSize(.75);
  h_fN->Draw("error same");

  h_fB->Scale(100.);
  h_fB->SetMarkerColor(kBlue);
  h_fB->SetLineColor(kBlue);
  h_fB->SetMarkerStyle(34);
  h_fB->SetMarkerSize(1.);
  h_fB->Draw("error same");
  
  TLatex lc;
  lc.SetTextSize(0.03);

  // draw CMS text    
  double xp = getPos(ptmin, ptmax, 0.75, 0);
  double yp = getPos(0, 100, 0.95, 0);
  lc.DrawLatex(xp, yp, "CMS");
  // draw L
  yp = getPos(0, 100, 0.9, 0);
  lc.DrawLatex(xp, yp, "#bf{L = 103.3 fb^{-1}}");
  // draw sqrt(s)
  yp = getPos(0, 100, 0.85, 0);
  lc.DrawLatex(xp, yp, "#bf{#sqrt{s} = 13 TeV}");
  // draw y
  yp = getPos(0, 100, 0.75, 0);
  lc.DrawLatex(xp, yp, "#bf{|#it{y}| < 1.2}");

  TLegend *leg = new TLegend(0.15, 0.7, 0.65, 0.85);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  leg->AddEntry(h_fB, "Base f_{bkg}", "pl");
  leg->AddEntry(h_fN, "New f_{bkg}", "pl");
  leg->Draw();


  c->SaveAs("f_comp.pdf");
  c->Clear();
  c->Destructor();

}
