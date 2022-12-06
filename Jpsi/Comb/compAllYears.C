// code to plot the fit results

void compAllYears()
{
  // get the fit results
  // get lambda values for each bin
  const int n_lbl = 2;
  string lbl[] = {"J", "NP"};
  string sv_lbl[] = {"PR", "NP"};
  string nm_lbl[] = {"prompt", "non-prompt"};
  
  TFile *fInd6 = new TFile("../2016/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth6 = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lth6[i_t] = (TGraphErrors*)fInd6->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd6->Close();
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

  for(int i_l = 0; i_l < n_lbl; i_l++) {

    // draw lambda_th(pT)
    TH1F *fl = c->DrawFrame(11, -1, 120, 1);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle("#lambda_{#theta}");
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("%s J/#psi #lambda_{#theta}", nm_lbl[i_l].c_str()));

    graph_lth6[i_l]->SetLineColor(kOrange+1);
    graph_lth6[i_l]->SetMarkerSize(.5);
    graph_lth6[i_l]->SetMarkerStyle(20);
    graph_lth6[i_l]->SetMarkerColor(kOrange+1);
    graph_lth6[i_l]->SetMarkerColor(kOrange+1);
    graph_lth6[i_l]->Draw("p same");
    
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

    TLine *zero = new TLine(11, 0, 120, 0);
    zero->SetLineColor(kBlack);
    zero->SetLineStyle(kDashed);
    zero->Draw();

    TLegend *leg = new TLegend(0.82, 0.75, 0.97, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(graph_lth6[i_l], "2016", "pl");
    leg->AddEntry(graph_lth7[i_l], "2017", "pl");
    leg->AddEntry(graph_lth8[i_l], "2018", "pl");
    leg->Draw();
  
    c->SaveAs(Form("plots/par_lth%s_allY.pdf", sv_lbl[i_l].c_str()));
    c->Clear();
  }
  c->Destructor();
}
