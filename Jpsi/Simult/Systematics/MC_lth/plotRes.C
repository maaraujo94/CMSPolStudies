// code to plot the fit results

const int n_colors = 50;
int st_col;

int do_round(double val)
{
  int valR = (int)val;
  if (val-valR > 0.5) return valR+1;
  else return valR;
}

void init_color()
{
  const int n_stops = 5;

  float stops[n_stops] = {0, 0.25, 0.5, 0.75, 1.0};
  float red[n_stops] = {1,1,0,0,0.9};
  float green[n_stops] = {0,0.6,0.7,0.5,0.1};
  float blue[n_stops] = {0,0,0,1,0};

  // For each defined gradient...
  for (int g = 1; g < n_stops; g++) {
      // create the colors...
      int nColorsGradient = (Int_t) (floor(n_colors*stops[g]) - floor(n_colors*stops[g-1]));
      for (int c = 0; c < nColorsGradient; c++) {
	TColor *col = new TColor( Float_t(red[g-1]   + c * (red[g]   - red[g-1])  / nColorsGradient),
				 Float_t(green[g-1] + c * (green[g] - green[g-1])/ nColorsGradient),
				 Float_t(blue[g-1]  + c * (blue[g]  - blue[g-1]) / nColorsGradient),
				  1);
	if(c == 0 && g == 1) st_col = col->GetNumber();
      }
  }

}

void plotRes()
{
  // get the lth values
  double lth[21];
  for(int i = 0; i < 21; i++) {
    lth[i] = -1.+i*0.1;
  }
  
  // get the histo limits
  TFile *fIn = new TFile("files/histoStore.root");
  TH2D* rHist;
  fIn->GetObject("MCH_10", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get A, lambda, chiProb values for each bin
  TFile *fInd = new TFile("files/finalFitRes.root");
  TGraphErrors **graph_A = new TGraphErrors*[21];
  TGraphErrors **graph_lth = new TGraphErrors*[21];
  TGraph **graph_chi = new TGraph*[21];
  for(int i_t = 0; i_t < 21; i_t++) {
    graph_A[i_t] = (TGraphErrors*)fInd->Get(Form("graph_A_%d", i_t));
    graph_lth[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%d", i_t));
    graph_chi[i_t] = (TGraph*)fInd->Get(Form("graph_chiP_%d", i_t));
  }    
  fInd->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);
  c->SetLeftMargin(0.1);
  
  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0]-5, -1.1, pTBins[nBinspT]+5, 1.35);
  fl->SetXTitle("#it{p}_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->GetXaxis()->SetLabelOffset(0.015);
  fl->GetXaxis()->SetTitleOffset(1.3);
  fl->GetXaxis()->CenterTitle(true);

  init_color();
  double n_max = 21;

  int st[] = {20, 21, 24, 25};
  for(int i = 0; i < 21; i++) {
    float cv = (n_colors-1)/n_max * i + st_col;
    int gc = do_round(cv);
    graph_lth[i]->SetLineColor(gc);
    graph_lth[i]->SetMarkerColor(gc);
    // graph_lth[i]->SetMarkerStyle(st[i%4]);
    graph_lth[i]->SetMarkerStyle(20);
    graph_lth[i]->SetMarkerSize(.75);
    graph_lth[i]->Draw("p same");

    TLine *cons = new TLine(pTBins[0]-5, -1+0.1*i, pTBins[nBinspT]+5, -1+0.1*i);
    cons->SetLineColor(gc);
    cons->SetLineStyle(kDashed);
    cons->Draw();
  }

  // need a all-black graph for plotting
  TGraphErrors *g_leg = (TGraphErrors*)graph_lth[0]->Clone("g_leg");
  g_leg->SetMarkerColor(kBlack);
  g_leg->SetLineColor(kBlack);
  TF1 *flin = new TF1("flin", "2", pTBins[0]-5, pTBins[nBinspT]+5);
  flin->SetLineColor(kBlack);
  flin->SetLineStyle(kDashed);
  flin->Draw("same");
  
  TLegend *legAll = new TLegend(0.7, 0.875, 1., 0.975);
  legAll->SetTextSize(0.03);
  legAll->SetBorderSize(0);
  legAll->SetFillColorAlpha(kWhite,0);
  legAll->AddEntry(flin, "Generated", "l");
  legAll->AddEntry(g_leg, "Measured", "pl");
  legAll->Draw();


  c->SaveAs("plots/ratioFinal/par_lth.pdf");
  c->Clear();

  // also draw the deviation from expected
  double dv[21][nBinspT];
  for(int i = 0; i < nBinspT; i++) {
    for(int j = 0; j < 21; j++) {
      dv[j][i] = graph_lth[j]->GetY()[i]-lth[j];
    }
  }
  TGraphErrors **graph_dlth = new TGraphErrors*[21];
  for(int i = 0; i < 21; i++) {
    graph_dlth[i] = new TGraphErrors(nBinspT, graph_lth[i]->GetX(), dv[i], graph_lth[i]->GetEX(), graph_lth[i]->GetEY());
  }
  
  TH1F *fdl = c->DrawFrame(pTBins[0]-5, -0.029, pTBins[nBinspT]+5, 0.029);
  fdl->SetXTitle("#it{p}_{T} (GeV)");
  fdl->SetYTitle("#Delta#lambda_{#theta}");
  fdl->GetYaxis()->SetTitleOffset(1.4);
  fdl->GetYaxis()->SetLabelOffset(0.01);
  fdl->GetXaxis()->SetLabelOffset(0.015);
  fdl->GetXaxis()->SetTitleOffset(1.3);
  fdl->GetXaxis()->CenterTitle(true);

  int n_plot[] = {8, 12, 14};
  for(int i = 0; i < 21; i++) {
    for(int i_p = 0; i_p < 3; i_p++) {
      if(i == n_plot[i_p]) {
	float cv = (n_colors-1)/n_max * i + st_col;
	int gc = do_round(cv);
	graph_dlth[i]->SetLineColor(gc);
	graph_dlth[i]->SetMarkerColor(gc);
	graph_dlth[i]->SetMarkerStyle(st[i%4]);
	graph_dlth[i]->SetMarkerSize(.75);
	graph_dlth[i]->Draw("p same");
      }
    }
  }

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT]+5, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg = new TLegend(0.7, 0.8, 1., 0.95);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  for(int i = 0; i < 21; i++) {
    for(int i_p = 0; i_p < 3; i_p++) {
      if(i == n_plot[i_p]) {
	leg->AddEntry(graph_dlth[i], Form("#lambda_{#theta} = %.1f", lth[i]), "pl");
      }
    }
  }
  leg->Draw();
    
  c->SaveAs("plots/ratioFinal/par_dlth.pdf");
  c->Clear();

  // draw A(pT)
  /* c->SetLogy();
  TH1F *fa = c->DrawFrame(pTBins[0], 1e-2, pTBins[nBinspT], 6e-1);
  fa->SetXTitle("p_{T} (GeV)");
  fa->SetYTitle("A");
  fa->GetYaxis()->SetTitleOffset(1.3);
  fa->GetYaxis()->SetLabelOffset(0.01);
  fa->SetTitle("Run 2 A");

  // combine both lambda_th distributions
  for(int i = 0; i < 21; i++) {
    graph_A[i]->SetLineColor(col[i]);
    graph_A[i]->SetMarkerColor(col[i]);
    graph_A[i]->Draw("p same");
  }

  c->SaveAs("plots/ratioFinal/par_A.pdf");
  c->Clear();

  // draw chiProb(pT)
  c->SetLogy(0);
  TH1F *fc = c->DrawFrame(pTBins[0], 0, pTBins[nBinspT], 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("P(#chi^{2}, ndf)");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);
  fc->SetTitle("Run 2 P(#chi^{2}, ndf)");

  // combine both lambda_th distributions
  for(int i = 0; i < 21; i++) {
    graph_chi[i]->SetLineColor(col[i]);
    graph_chi[i]->SetMarkerColor(col[i]);
    graph_chi[i]->SetMarkerStyle(20);
    graph_chi[i]->SetMarkerSize(.75);
    graph_chi[i]->Draw("p same");
  }

  c->SaveAs("plots/ratioFinal/par_chiP.pdf");
  c->Clear();*/
  c->Destructor();
  
  fIn->Close();


}
