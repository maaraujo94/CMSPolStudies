void plot_EffAcc()
{
  // get MC dist from histoStore
  TFile *fin2 = new TFile("../PR_fit/files/histoStore.root");
  TH2D *h_MC2 = (TH2D*)fin2->Get("MCH");

  // get pT bins
  const double *yBins = h_MC2->GetYaxis()->GetXbins()->GetArray();

  TH1D **h1_MC = new TH1D*[3];
  int nb[] = {1,4,7};
  for(int i = 0; i< 3; i++) {
    h1_MC[i] = (TH1D*)h_MC2->ProjectionX(Form("h1_MC_%d",i), nb[i], nb[i]);
  }

  TCanvas *c = new TCanvas("","",900,900);
  c->SetLeftMargin(0.11);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.02);

  int colpt[] = {kBlack, kRed+1, kBlue};

  for(int i = 0; i < 3; i++) {
    h1_MC[i]->Scale(1/h1_MC[i]->GetBinContent(1));
    h1_MC[i]->SetLineColor(colpt[i]);
    h1_MC[i]->SetMarkerColor(colpt[i]);
    h1_MC[i]->SetMarkerStyle(20);
    h1_MC[i]->SetMarkerSize(.75);
    h1_MC[i]->SetStats(0);
    h1_MC[i]->SetTitle("");
    h1_MC[i]->GetXaxis()->SetTitle("|cos #theta_{HX}|");
    //    h1_MC[i]->GetYaxis()->SetTitle("Det acc * eff");
    h1_MC[i]->GetYaxis()->SetTitleOffset(1.6);
    h1_MC[i]->GetYaxis()->SetLabelOffset(0.01);
    h1_MC[i]->GetXaxis()->SetLabelOffset(0.01);
    h1_MC[i]->GetXaxis()->SetTitleOffset(1.25);
    h1_MC[i]->GetXaxis()->CenterTitle(true);
    h1_MC[i]->SetMinimum(0);
    h1_MC[i]->SetMaximum(1.5);

    if(i == 0)
      h1_MC[i]->Draw("error");
    else
      h1_MC[i]->Draw("same");
  }

  TLatex lcpt;
  double x = 0.05, y = 1.3;
  lcpt.SetTextSize(0.04);
  lcpt.DrawLatex(x, y, "#bf{2018 #psi(2S)}");

  lcpt.SetTextSize(0.03);
  y = 0.8;
  lcpt.SetTextColor(colpt[0]);
  y -= 0.09;
  lcpt.DrawLatex(x, y, Form("#bf{[%.0f, %.0f] GeV}", yBins[nb[0]], yBins[nb[0]+1]));
  lcpt.SetTextColor(colpt[1]);
  y -= 0.09;
  lcpt.DrawLatex(x, y, Form("#bf{[%.0f, %.0f] GeV}", yBins[nb[1]], yBins[nb[1]+1]));
  lcpt.SetTextColor(colpt[2]);
  y -= 0.09;
  lcpt.DrawLatex(x, y, Form("#bf{[%.0f, %.0f] GeV}", yBins[nb[2]], yBins[nb[2]+1]));
  


  c->SaveAs("plots/eff_acc.pdf");
  c->Clear();
  c->Destructor();

  fin2->Close();
}
