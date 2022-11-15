// code to plot the fit results

#import "../Simult/ptbins.C"

void compYears()
{
  // get the fit results
  // get lambda values for each bin
  const int n_lbl = 3;
  string lbl[] = {"J", "NP", "Data"};
  string sv_lbl[] = {"", "NP", "Peak"};
  string nm_lbl[] = {"prompt", "non-prompt", "peak"};
  
  TFile *fInd7 = new TFile("../2017/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth7 = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lth7[i_t] = (TGraphErrors*)fInd7->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd7->Close();
  TFile *fInd8 = new TFile("../2018/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth8 = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lth8[i_t] = (TGraphErrors*)fInd8->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd8->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);

  for(int i_l = 0; i_l < 3; i_l++) {

    // get the unc band
    double unc_v[nPtBins], unc_y[nPtBins];
    for(int i = 0; i < nPtBins; i++) {
      unc_v[i] = 0.012;
      unc_y[i] = 0.5*(graph_lth8[i_l]->GetY()[i]+graph_lth7[i_l]->GetY()[i]);
      cout << i << " " << (graph_lth7[i_l]->GetY()[i]-graph_lth8[i_l]->GetY()[i])/(graph_lth8[i_l]->GetEY()[i]+graph_lth7[i_l]->GetEY()[i]) << endl;
    }
    TGraphErrors *g_unc = new TGraphErrors(nPtBins, graph_lth7[i_l]->GetX(), unc_y, graph_lth7[i_l]->GetEX(), unc_v);
    
    // draw lambda_th(pT)
    TH1F *fl = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle("#lambda_{#theta}");
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("%s #psi(2S) #lambda_{#theta}", nm_lbl[i_l].c_str()));

    graph_lth7[i_l]->SetLineColor(kBlue);
    graph_lth7[i_l]->SetMarkerSize(.5);
    graph_lth7[i_l]->SetMarkerStyle(20);
    graph_lth7[i_l]->SetMarkerColor(kBlue);
    graph_lth7[i_l]->SetMarkerColor(kBlue);
    graph_lth7[i_l]->Draw("p same");
    graph_lth8[i_l]->SetLineColor(kRed);
    graph_lth8[i_l]->SetMarkerStyle(20);
    graph_lth8[i_l]->SetMarkerSize(.5);
    graph_lth8[i_l]->SetMarkerColor(kRed);
    graph_lth8[i_l]->Draw("p same");

    TLine *zero = new TLine(ptBins[0]-5, 0, ptBins[nPtBins], 0);
    zero->SetLineColor(kBlack);
    zero->SetLineStyle(kDashed);
    zero->Draw();

    TLegend *leg = new TLegend(0.82, 0.75, 0.97, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(graph_lth7[i_l], "2017", "pl");
    leg->AddEntry(graph_lth8[i_l], "2018", "pl");
    //    leg->AddEntry(g_unc, "#sigma", "f");
    leg->Draw();
  
    c->SaveAs(Form("plots/par_lth%s_years.pdf", sv_lbl[i_l].c_str()));
    
    g_unc->SetFillColorAlpha(kViolet-1, 0.5);
    g_unc->Draw("e3");

    c->SaveAs(Form("plots/par_lth%s_years_unc.pdf", sv_lbl[i_l].c_str()));
    c->Clear();
  }
  c->Destructor();
}
