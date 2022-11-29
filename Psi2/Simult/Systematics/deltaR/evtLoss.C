void evtLoss()
{
  // get baseline evts
  TFile *finB = new TFile("../../PR_fit/files/histoStore.root");
  TH2D *h_dataB = (TH2D*)finB->Get("PRH");
  h_dataB->SetDirectory(0);
  TH2D *h_mcB = (TH2D*)finB->Get("MCH");
  h_mcB->SetDirectory(0);
  finB->Close();

  // get tight cut events
  TFile *finT = new TFile("../../../Simult_dR1/PR_fit/files/histoStore.root");
  TH2D *h_dataT = (TH2D*)finT->Get("PRH");
  h_dataT->SetDirectory(0);
  TH2D *h_mcT = (TH2D*)finT->Get("MCH");
  h_mcT->SetDirectory(0);
  finT->Close();

  // get loose cut events
  TFile *finL = new TFile("../../../Simult_dR2/PR_fit/files/histoStore.root");
  TH2D *h_dataL = (TH2D*)finL->Get("PRH");
  h_dataL->SetDirectory(0);
  TH2D *h_mcL = (TH2D*)finL->Get("MCH");
  h_mcL->SetDirectory(0);
  finL->Close();

  //get binning
  int nBinsX = h_dataB->GetNbinsX(), nBinsY = h_dataB->GetNbinsY();
  const double *yBins = h_dataB->GetYaxis()->GetXbins()->GetArray();
  double minX = h_dataB->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_dataB->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/nBinsX;

  // calculate the event loss (in percentage)
  double diffL[nBinsY], diffT[nBinsY];
  double mdiffL[nBinsY], mdiffT[nBinsY];
  double za[nBinsY], eptBins[nBinsY], ptBins[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    ptBins[i] = 0.5*(yBins[i+1]+yBins[i]);
    eptBins[i] = 0.5*(yBins[i+1]-yBins[i]);
    za[i] = 0;

    double vB = h_dataB->Integral(1, nBinsX, i+1, i+1);
    double vL = h_dataL->Integral(1, nBinsX, i+1, i+1);
    double vT = h_dataT->Integral(1, nBinsX, i+1, i+1);

    diffL[i] = 1-vL/vB;
    diffT[i] = 1-vT/vB;

    vB = h_mcB->Integral(1, nBinsX, i+1, i+1);
    vL = h_mcL->Integral(1, nBinsX, i+1, i+1);
    vT = h_mcT->Integral(1, nBinsX, i+1, i+1);

    mdiffL[i] = 1-vL/vB;
    mdiffT[i] = 1-vT/vB;
  }

  // define the tgraphs for the two sets of evts
  TGraphErrors *g_lossL = new TGraphErrors(nBinsY, ptBins, diffL, eptBins, za);
  TGraphErrors *g_lossT = new TGraphErrors(nBinsY, ptBins, diffT, eptBins, za);
  TGraphErrors *g_mlossL = new TGraphErrors(nBinsY, ptBins, mdiffL, eptBins, za);
  TGraphErrors *g_mlossT = new TGraphErrors(nBinsY, ptBins, mdiffT, eptBins, za);

  TCanvas *c = new TCanvas("", "", 700, 700);
  TH1F *fh = c->DrawFrame(yBins[0], 0, yBins[nBinsY], 1.);
  fh->SetXTitle("p_{T} (GeV)");
  fh->SetYTitle("fraction of lost events");
  fh->GetYaxis()->SetTitleOffset(1.3);
  fh->GetYaxis()->SetLabelOffset(0.01);
  fh->SetTitle("Lost events after #DeltaR cuts");

  g_lossL->SetMarkerStyle(20);
  g_lossL->SetMarkerSize(.75);
  g_lossL->SetMarkerColor(kBlue);
  g_lossL->SetLineColor(kBlue);
  g_lossL->Draw("psame");

  g_mlossL->SetMarkerStyle(20);
  g_mlossL->SetMarkerSize(.75);
  g_mlossL->SetMarkerColor(kBlue);
  g_mlossL->SetLineStyle(kDashed);
  g_mlossL->SetLineColor(kBlue);
  g_mlossL->Draw("psame");

  g_lossT->SetMarkerStyle(20);
  g_lossT->SetMarkerSize(.75);
  g_lossT->SetMarkerColor(kRed);
  g_lossT->SetLineColor(kRed);
  g_lossT->Draw("psame");

  g_mlossT->SetMarkerStyle(20);
  g_mlossT->SetMarkerSize(.75);
  g_mlossT->SetMarkerColor(kRed);
  g_mlossT->SetLineColor(kRed);
  g_mlossT->SetLineStyle(kDashed);
  g_mlossT->Draw("psame");

  TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(g_lossL, "#DeltaR>0.15", "pl");
  leg->AddEntry(g_lossT, "#DeltaR>0.17", "pl");
  leg->Draw();

  TLine *cOff = new TLine(50, 0, 50, 1);
  cOff->SetLineStyle(kDashed);
  cOff->Draw();

  c->SaveAs("plots/evtLoss.pdf");
  c->Clear();
  c->Destructor();

}
