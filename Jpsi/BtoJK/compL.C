void compL()
{
  // get the graphs from B->J/psi K and from inclusive NP J/psi
  TFile *fB = new TFile("files/lambdas.root");
  TGraphAsymmErrors* lambdaN_0 = (TGraphAsymmErrors*)fB->Get("lambdaN_0");
  TGraphAsymmErrors* lambdaN_1 = (TGraphAsymmErrors*)fB->Get("lambdaN_1");
  TGraphAsymmErrors* lambdaN_2 = (TGraphAsymmErrors*)fB->Get("lambdaN_2");
  fB->Close();

  TFile *fNP = new TFile("../2018/NP_fit/files/fit_res_1d.root");
  TGraphErrors* graph_lth = (TGraphErrors*)fNP->Get("graph_lambda"); 
  TGraphErrors* graph_lth_hpt = (TGraphErrors*)fNP->Get("graph_lambda_hpt");
  fNP->Close();

  // get the binning from the ratio hist
  TFile *fIn = new TFile("../2018/NP_fit/files/ratioHist.root");
  TH2D* rHist;
  TH2D* rHist_hpt;
  fIn->GetObject("ratioHist_ab", rHist);
  fIn->GetObject("ratioHist_ab_hpt", rHist_hpt);

  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  int nBinspT_hpt = rHist_hpt->GetNbinsY();
  const double *pTBins_hpt = rHist_hpt->GetYaxis()->GetXbins()->GetArray();

  int ntot = nBinspT+nBinspT_hpt;
  double xvals[ntot], ex[ntot];
  double yvals[ntot], ey[ntot];

  for(int i = 0; i < ntot; i++)
    {
      cout << i << endl;
      if(i < nBinspT) {
	xvals[i] = graph_lth->GetX()[i];
	yvals[i] = graph_lth->GetY()[i];
	ex[i] = graph_lth->GetEX()[i];
	ey[i] = graph_lth->GetEY()[i];
      }
      else {
	xvals[i] = graph_lth_hpt->GetX()[i-nBinspT];
	yvals[i] = graph_lth_hpt->GetY()[i-nBinspT];
	ex[i] = graph_lth_hpt->GetEX()[i-nBinspT];
	ey[i] = graph_lth_hpt->GetEY()[i-nBinspT];
      }
    }
  TGraphErrors *graph_lth_tot = new TGraphErrors(ntot, xvals, yvals, ex, ey);
  
  TCanvas *c = new TCanvas("", "", 700, 700);

  TH1F *fl = c->DrawFrame(pTBins[0], -0.5, pTBins_hpt[nBinspT_hpt], 0.5);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("#lambda_{#theta} (NP)");

  // combine both lambda_th distributions
  graph_lth_tot->SetLineColor(kBlack);
  graph_lth_tot->SetMarkerColor(kBlack);
  graph_lth_tot->SetFillColorAlpha(kBlack, 0.5);
  graph_lth_tot->Draw("le3");
  
  lambdaN_0->SetLineColor(kBlue);
  lambdaN_0->SetMarkerColor(kBlue);
  lambdaN_0->SetFillColorAlpha(kBlue, 0.5);
  lambdaN_0->Draw("le3");
  lambdaN_1->SetLineColor(kBlue+1);
  lambdaN_1->SetMarkerColor(kBlue+1);
  lambdaN_1->SetFillColorAlpha(kBlue+1, 0.5);
  lambdaN_1->Draw("le3");
  lambdaN_2->SetLineColor(kBlue+2);
  lambdaN_2->SetMarkerColor(kBlue+2);
  lambdaN_2->SetFillColorAlpha(kBlue+2, 0.5);
  lambdaN_2->Draw("le3");

  TLegend *leg_l = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg_l->SetTextSize(0.03);
  leg_l->AddEntry(graph_lth, "NP/MC fit result", "pl");
  leg_l->AddEntry(lambdaN_0, lambdaN_0->GetTitle(), "pl");
  leg_l->AddEntry(lambdaN_1, lambdaN_1->GetTitle(), "pl");
  leg_l->AddEntry(lambdaN_2, lambdaN_2->GetTitle(), "pl");
  leg_l->Draw();

  c->SaveAs("plots/lth_comp.pdf");
  
  fIn->Close();
}
