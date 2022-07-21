// code to plot the fit results

void plotRes()
{
  // get the histo limits
  TFile *fIn = new TFile("files/bkgSubRes.root");
  TH2D* rHist;
  fIn->GetObject("h_Data", rHist);
  
  int nBinspT = rHist->GetNbinsY();
  const double *pTBins = rHist->GetYaxis()->GetXbins()->GetArray();
  
  // get the fit results
  // get A, lambda, chiProb values for each bin
  string lbl[] = {"Data", "NP", "PR", "J"};
  TFile *fInd = new TFile("files/finalFitRes.root");
  TGraphErrors **graph_A = new TGraphErrors*[4];
  TGraphErrors **graph_lth = new TGraphErrors*[4];
  TGraph **graph_chi = new TGraph*[4];
  for(int i_t = 0; i_t < 4; i_t++) {
    graph_A[i_t] = (TGraphErrors*)fInd->Get(Form("graph_A_%s", lbl[i_t].c_str()));
    graph_lth[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
    graph_chi[i_t] = (TGraph*)fInd->Get(Form("graph_chiP_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();

    // get the final lth from the base SB/MC fit
  TFile *fIndB = new TFile("../../PR_fit/files/finalFitRes.root");
  TGraphErrors *graph_lthBase = (TGraphErrors*)fIndB->Get("graph_lambda_J");
  fIndB->Close();
  
  for(int i = 0; i < graph_lthBase->GetN(); i++) {
    cout << graph_lthBase->GetX()[i] << " " << graph_lthBase->GetY()[i] << " " << graph_lth[3]->GetY()[i] << endl;
  }

  // slightly shifting the central x values
  int nB = graph_lthBase->GetN();
  double xB[nB], exlB[nB], exhB[nB];
  for(int i = 0; i < nB; i++) {
    xB[i] = graph_lthBase->GetX()[i]+0.5;
    exlB[i] = graph_lthBase->GetEX()[i]+0.5;
    exhB[i] = graph_lthBase->GetEX()[i]-0.5;
  }
  TGraphAsymmErrors *graph_lthB = new TGraphAsymmErrors(nB, xB, graph_lthBase->GetY(), exlB, exhB, graph_lthBase->GetEY(), graph_lthBase->GetEY());

  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(pTBins[0]-5, -1, pTBins[nBinspT], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("2017 #lambda_{#theta}");

  int col[] = {kViolet, kRed, kBlack, kBlue};
  for(int i = 0; i < 4; i++) {
    graph_lth[i]->SetLineColor(col[i]);
    graph_lth[i]->SetMarkerColor(col[i]);
    graph_lth[i]->Draw("p same");
  }

  TLine *zero = new TLine(pTBins[0]-5, 0, pTBins[nBinspT], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  TLine *trans1 = new TLine(46, -1, 46, 1);
  trans1->SetLineColor(kBlack);
  trans1->SetLineStyle(kDashed);
  trans1->Draw();

  TLegend *leg = new TLegend(0.65, 0.12, 0.9, 0.32);
  leg->SetTextSize(0.03);
  leg->AddEntry(graph_lth[0], "total", "pl");
  leg->AddEntry(graph_lth[1], "NP contrib", "pl");
  leg->AddEntry(graph_lth[2], "prompt", "pl");
  leg->AddEntry(graph_lth[3], "prompt #psi(2S)", "pl");
  leg->Draw();
  
  c->SaveAs("plots/ratioFinal/par_lth.pdf");
  c->Clear();

  // draw just final lambda_th(pT) - comp btw std, alt
  double val[nBinspT];
  for(int i = 0; i < nBinspT; i++) { 
    val[i] = graph_lth[3]->GetY()[i] - graph_lthBase->GetY()[i];
  }
  TGraphErrors *g_lthD = new TGraphErrors(nBinspT, graph_lthBase->GetX(), val, graph_lthBase->GetEX(), graph_lthBase->GetEY());

  double d_lim = 0.4;

  TH1F *fl2 = c->DrawFrame(pTBins[0]-5, -d_lim, pTBins[nBinspT], d_lim);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#delta#lambda_{#theta}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("2017 #delta#lambda_{#theta} (prompt #psi(2S))");

  g_lthD->SetLineColor(kBlack);
  g_lthD->SetMarkerColor(kBlack);
  g_lthD->SetMarkerStyle(20);
  g_lthD->SetMarkerSize(.5);
  g_lthD->Draw("p same");

  zero->Draw();
  TLine *trans1D = new TLine(46, -d_lim, 46, d_lim);
  trans1D->SetLineColor(kBlack);
  trans1D->SetLineStyle(kDashed);
  trans1D->Draw();
  
  c->SaveAs("par_lth_F.pdf");
  c->Clear();

  // draw A(pT)
  c->SetLogy();
  TH1F *fa = c->DrawFrame(pTBins[0], 1e-2, pTBins[nBinspT], 4e-1);
  fa->SetXTitle("p_{T} (GeV)");
  fa->SetYTitle("A");
  fa->GetYaxis()->SetTitleOffset(1.3);
  fa->GetYaxis()->SetLabelOffset(0.01);
  fa->SetTitle("2017 A");

  // combine both lambda_th distributions
  for(int i = 0; i < 4; i++) {
    graph_A[i]->SetLineColor(col[i]);
    graph_A[i]->SetMarkerColor(col[i]);
    graph_A[i]->Draw("p same");
  }

  TLine *trans1_A = new TLine(46, 1e-2, 46, 6e-1);
  trans1_A->SetLineColor(kBlack);
  trans1_A->SetLineStyle(kDashed);
  trans1_A->Draw();

  c->SaveAs("plots/ratioFinal/par_A.pdf");
  c->Clear();

  // draw chiProb(pT)
  c->SetLogy(0);
  TH1F *fc = c->DrawFrame(pTBins[0], 0, pTBins[nBinspT], 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("P(#chi^{2}, ndf)");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);
  fc->SetTitle("2017 P(#chi^{2}, ndf)");

  // combine both lambda_th distributions
  for(int i = 0; i < 4; i++) {
    graph_chi[i]->SetLineColor(col[i]);
    graph_chi[i]->SetMarkerColor(col[i]);
    graph_chi[i]->SetMarkerStyle(20);
    graph_chi[i]->SetMarkerSize(.75);
    graph_chi[i]->Draw("p same");
  }

  TLine *trans1_C = new TLine(46, 0, 46, 1);
  trans1_C->SetLineColor(kBlack);
  trans1_C->SetLineStyle(kDashed);
  trans1_C->Draw();

  c->SaveAs("plots/ratioFinal/par_chiP.pdf");
  c->Clear();
  c->Destructor();
  
  fIn->Close();


}
