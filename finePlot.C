// code to plot the fine bins of all 3 ratios

void finePlot()
{
  // getting PR hist
  TFile *fPR = new TFile("2d_fit/files/ratioHist.root");
  TH2D *hPR = (TH2D*)fPR->Get("ratioHist_ab");
  hPR->SetDirectory(0);
  fPR->Close();
  // getting NP hist
  TFile *fNP = new TFile("NP_fit/files/ratioHist.root");
  TH2D *hNP = (TH2D*)fNP->Get("ratioHist_ab");
  hNP->SetDirectory(0);
  fNP->Close();
  // getting PR hist
  TFile *fPN = new TFile("PR_NP_fit/files/ratioHist.root");
  TH2D *hPN = (TH2D*)fPN->Get("ratioHist_ab");
  hPN->SetDirectory(0);
  fPN->Close();

  int nBinsY = hPR->GetNbinsY();
  
  // get the 1D dists
  TH1D *pPR[nBinsY], *pNP[nBinsY], *pPN[nBinsY];
  for(int i = 1; i <= nBinsY; i++) {
    pPR[i-1] = hPR->ProjectionX(Form("PR_bin%d", i), i, i);
    pNP[i-1] = hNP->ProjectionX(Form("NP_bin%d", i), i, i);
    pPN[i-1] = hPN->ProjectionX(Form("PR_NP_bin%d", i), i, i);
  }

  // plot the dists
  TCanvas *c = new TCanvas("", "", 700, 700);
  double norm;
  for(int i = 0; i < nBinsY; i++) {
    // scale dists
    norm = pPR[i]->GetBinContent(1);
    pPR[i]->Scale(1./norm);
    norm = pNP[i]->GetBinContent(1);
    pNP[i]->Scale(1./norm);
    norm = pPN[i]->GetBinContent(1);
    pPN[i]->Scale(1./norm);

    double mL = max(pPR[i]->GetMaximum(), pNP[i]->GetMaximum())+0.1;
    pPR[i]->SetMaximum(mL);
    pNP[i]->SetMaximum(mL);

    pPR[i]->Draw();
    c->SaveAs(Form("plots/PR_bin%d.pdf", i+1));
    c->Clear();

    pNP[i]->Draw();
    c->SaveAs(Form("plots/NP_bin%d.pdf", i+1));
    c->Clear();

    pPN[i]->Draw();
    c->SaveAs(Form("plots/PN_bin%d.pdf", i+1));
    c->Clear();
  }
  
  // getting PR/NP hist
  TFile *fPNc = new TFile("PR_NP_fit/files/ratioHist.root");
  TH2D *hPNc = (TH2D*)fPNc->Get("cHist_ab");
  hPNc->SetDirectory(0);
  fPNc->Close();

  nBinsY = hPNc->GetNbinsY();
  
  // get the 1D dists
  TH1D *pPNc[nBinsY];
  for(int i = 1; i <= nBinsY; i++) {
    pPNc[i-1] = hPNc->ProjectionX(Form("PR_NP_c_bin%d", i), i, i);
  }

  for(int i = 0; i < nBinsY; i++) {
    pPNc[i]->SetStats(0);
    pPNc[i]->SetMaximum(pPNc[i]->GetBinContent(1)*1.4);
    pPNc[i]->SetTitle(Form("PR/NP coarse bin %d", i+1));
    pPNc[i]->Draw();
    c->SaveAs(Form("plots/PN_c_bin%d.pdf", i+1));
    c->Clear();
  }
  c->Destructor();
  
}


