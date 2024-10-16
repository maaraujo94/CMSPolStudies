// code to plot the fit results - NP case
// compares the 2d and 1d fits for HX (w/ and w/o ltp)

#import "../ptbins.C"

void plotRes_NP()
{
  // get the fit results - 2d
  TFile *fInd = new TFile("../NP_fit/files/finalFitRes.root");
  TGraphErrors **graph_lNP = new TGraphErrors*[3];
  string lbl[3] = {"theta", "phi", "thph"};
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lNP[i_t] = (TGraphErrors*)fInd->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd->Close();

  // get the fit results - 2d w/o ltp
  TFile *fIndz = new TFile("../NP_fit/files/finalFitRes_ltpz.root");
  TGraphErrors **graph_lNPz = new TGraphErrors*[3];
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lNPz[i_t] = (TGraphErrors*)fIndz->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fIndz->Close();
  
  //get the fit results - 1d
  // lth is immediate
  TFile *find_th = new TFile("../../Simult/NP_fit/files/finalFitRes.root");
  TGraphErrors *graph_lth = (TGraphErrors*)find_th->Get(Form("graph_lambda_NPc"));
  find_th->Close();
  // lphi need to get from beta
  TFile *find_ph = new TFile("../../Phi_fit/NP_fit/files/finalFitRes.root");
  TGraphErrors *graph_beta = (TGraphErrors*)find_ph->Get(Form("graph_lambda_NPc"));
  find_ph->Close();

  int nv = graph_beta->GetN();
  double v_lph[nv], e_lph[nv];
  for(int i = 0; i < nv; i++) {
    double v_c = graph_lth->GetY()[i];
    double v_b = graph_beta->GetY()[i];
    double e_c = graph_lth->GetEY()[i];
    double e_b = graph_beta->GetEY()[i];
    v_lph[i] = 0.5*v_b*(3.+v_c);
    e_lph[i] = sqrt(pow((3.+v_c)/2., 2)*pow(e_b,2) + pow(v_b/2.,2)*pow(e_c,2));
  }
  TGraphErrors *graph_lph = new TGraphErrors(nv, graph_beta->GetX(), v_lph, graph_beta->GetEX(), e_lph);
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.02);
  int col[] = {kBlack, kRed, kGreen+1};

  // draw lambda_th(pT)
  TH1F *flth = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
  flth->SetXTitle("p_{T} (GeV)");
  flth->SetYTitle("#lambda_{#theta}^{NP}");
  flth->GetYaxis()->SetTitleOffset(1.3);
  flth->GetYaxis()->SetLabelOffset(0.01);
  flth->SetTitle("");
    
  graph_lNP[0]->SetLineColor(col[0]);
  graph_lNP[0]->SetMarkerStyle(20);
  graph_lNP[0]->SetMarkerSize(.5);
  graph_lNP[0]->SetMarkerColor(col[0]);
  graph_lNP[0]->Draw("p same");

  graph_lNPz[0]->SetLineColor(col[0]);
  graph_lNPz[0]->SetMarkerStyle(24);
  graph_lNPz[0]->SetMarkerSize(.5);
  graph_lNPz[0]->SetMarkerColor(col[0]);
  graph_lNPz[0]->SetLineStyle(kDashed);
  graph_lNPz[0]->Draw("p same");

  graph_lth->SetLineColor(col[1]);
  graph_lth->SetMarkerStyle(25);
  graph_lth->SetMarkerSize(.5);
  graph_lth->SetMarkerColor(col[1]);
  graph_lth->SetLineStyle(kDashed);
  graph_lth->Draw("p same");

  TLine *zero = new TLine(ptBins[0]-5, 0, ptBins[nPtBins], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *legth = new TLegend(0.2, 0.15, 0.5, 0.3);
  legth->SetTextSize(0.03);
  legth->SetBorderSize(0);
  legth->SetFillColorAlpha(kWhite, 0);
  legth->AddEntry(graph_lNP[0], "#lambda_{#theta}^{HX} (2d)", "pl");
  legth->AddEntry(graph_lNPz[0], "#lambda_{#theta}^{HX} (2d w/o #lambda_{#theta#phi})", "pl");
  legth->AddEntry(graph_lth, "#lambda_{#theta}^{HX} (1d)", "pl");
  legth->Draw();
  
  c->SaveAs("lthNP_HX_comp.pdf");
  c->Clear();

  // now lambda_phi
  TH1F *flph = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
  flph->SetXTitle("p_{T} (GeV)");
  flph->SetYTitle("#lambda_{#phi}^{NP}");
  flph->GetYaxis()->SetTitleOffset(1.3);
  flph->GetYaxis()->SetLabelOffset(0.01);
  flph->SetTitle("");
    
  graph_lNP[1]->SetLineColor(col[0]);
  graph_lNP[1]->SetMarkerStyle(20);
  graph_lNP[1]->SetMarkerSize(.5);
  graph_lNP[1]->SetMarkerColor(col[0]);
  graph_lNP[1]->Draw("p same");

  graph_lNPz[1]->SetLineColor(col[0]);
  graph_lNPz[1]->SetMarkerStyle(24);
  graph_lNPz[1]->SetMarkerSize(.5);
  graph_lNPz[1]->SetMarkerColor(col[0]);
  graph_lNPz[1]->SetLineStyle(kDashed);
  graph_lNPz[1]->Draw("p same");

  graph_lph->SetLineColor(col[1]);
  graph_lph->SetMarkerStyle(25);
  graph_lph->SetMarkerSize(.5);
  graph_lph->SetMarkerColor(col[1]);
  graph_lph->SetLineStyle(kDashed);
  graph_lph->Draw("p same");

  zero->Draw();

  TLegend *legph = new TLegend(0.2, 0.15, 0.5, 0.3);
  legph->SetTextSize(0.03);
  legph->SetBorderSize(0);
  legph->SetFillColorAlpha(kWhite, 0);
  legph->AddEntry(graph_lNP[1], "#lambda_{#phi}^{HX} (2d)", "pl");
  legph->AddEntry(graph_lNPz[1], "#lambda_{#theta}^{HX} (2d w/o #lambda_{#theta#phi})", "pl");
  legph->AddEntry(graph_lph, "#lambda_{#phi}^{HX} (1d)", "pl");
  legph->Draw();

  c->SaveAs("lphNP_HX_comp.pdf");
  c->Clear();

  c->Destructor();
 
}
