void fNPcorr()
{
  // get the base f_NP
  TFile *inNP = new TFile("files/NPFrac.root");
  TH2D *h_fnp = (TH2D*)inNP->Get("h_fnp");
  h_fnp->SetDirectory(0);
  inNP->Close();
  
  // get the f_bkg^NP from NP_fit folder
  TFile *inSB = new TFile("../bkgFits/files/bkgFrac_NP.root");
  TH2D *h_fbkg = (TH2D*)inSB->Get("h_fbkg");
  h_fbkg->SetDirectory(0);
  inSB->Close();

  // get the binning
  int nBinsX = h_fnp->GetNbinsX(), nBinsY = h_fnp->GetNbinsY();
  const double *yBins = h_fnp->GetYaxis()->GetXbins()->GetArray();
  double minX = h_fnp->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_fnp->GetXaxis()->GetBinUpEdge(nBinsX);

  // get 1-f_bkg^NP
  TH2D *h_fc = new TH2D("h_fc", "h_fc", nBinsX, minX, maxX, nBinsY, yBins);
  for(int i_x = 0; i_x < nBinsX; i_x++) {
    for(int i_y = 0; i_y < nBinsY; i_y++) {
      h_fc->SetBinContent(i_x+1, i_y+1, 1.-h_fbkg->GetBinContent(i_x+1, i_y+1));
      h_fc->SetBinError(i_x+1, i_y+1, 0);
    }
  }

  // get corrected f_NP
  TH2D *h_fNPc = (TH2D*)h_fnp->Clone("h_fNPc");
  h_fNPc->Multiply(h_fc);

  // plot the comparison
  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1D *h_fnp1d = h_fnp->ProjectionY("h_fnp1d", 1, 1);
  TH1D *h_fbkg1d = h_fbkg->ProjectionY("h_fbkg1d", 1, 1);
  TH1D *h_fnpc1d = h_fNPc->ProjectionY("h_fnpc1d", 1, 1);

  h_fnp1d->Scale(100.);
  h_fbkg1d->Scale(100.);
  h_fnpc1d->Scale(100.);

  cout << h_fnp1d->GetBinError(1) << " " << h_fbkg1d->GetBinError(1) << " "<< h_fnpc1d->GetBinError(1) << endl;
  
  h_fnp1d->SetStats(0);
  h_fnp1d->SetMinimum(0);
  h_fnp1d->SetMaximum(50);
  h_fnp1d->GetXaxis()->SetTitle("p_{T} (GeV)");
  h_fnp1d->GetYaxis()->SetTitle("f (%)");
  h_fnp1d->GetYaxis()->SetTitleOffset(1.3);
  h_fnp1d->GetYaxis()->SetLabelOffset(0.01);
  h_fnp1d->SetTitle("2017 f_{NP}^{corr}");
  h_fnp1d->SetLineColor(kRed);
  h_fnp1d->SetMarkerColor(kRed);
  h_fnp1d->SetLineStyle(kDashed);
  h_fnp1d->SetMarkerStyle(20);
  h_fnp1d->SetMarkerSize(.75);
  h_fnp1d->Draw("error");
  
  h_fbkg1d->SetLineColor(kGreen+1);
  h_fbkg1d->SetMarkerColor(kGreen+1);
  h_fbkg1d->SetMarkerStyle(20);
  h_fbkg1d->SetMarkerSize(.75);
  h_fbkg1d->Draw("error same");
    
  h_fnpc1d->SetLineColor(kRed);
  h_fnpc1d->SetMarkerColor(kRed);
  h_fnpc1d->SetMarkerStyle(20);
  h_fnpc1d->SetMarkerSize(.75);
  h_fnpc1d->Draw("error same");

  TLegend *leg = new TLegend(0.7, 0.6, 0.9, 0.9);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_fnp1d, "f_{NP}", "pl");
  leg->AddEntry(h_fbkg1d, "f_{bkg}^{NP}", "pl");
  leg->AddEntry(h_fnpc1d, "f_{NP}^{c}", "pl");
  leg->Draw();
  
  c->SaveAs("plots/f_NP_corr.pdf");
  c->Destructor();
  
  TFile *fout = new TFile("files/NPFrac.root", "update");
  h_fNPc->Write(0, TObject::kOverwrite);
  fout->Close();
}
