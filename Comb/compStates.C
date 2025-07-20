// code to plot the fit results

void compStates()
{
  // get the fit results
  // get lambda values for each bin
  const int n_lbl = 2;
  string lbl[] = {"J", "NP"}; 
  string sv_lbl[] = {"PR", "NP"};
  string nm_lbl[] = {"prompt", "non-prompt"};
  
  TFile *fIndJ = new TFile("../Jpsi/Simult/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lthJ = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lthJ[i_t] = (TGraphErrors*)fIndJ->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fIndJ->Close();
  TFile *fIndP = new TFile("../Psi2/Simult/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lthP = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lthP[i_t] = (TGraphErrors*)fIndP->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fIndP->Close();

  TFile *fIndJ_f = new TFile("../Jpsi/Simult/Systematics/mainDiffs/files/finalUnc.root");
  TGraphAsymmErrors **graph_lthJ_f = new TGraphAsymmErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lthJ_f[i_t] = (TGraphAsymmErrors*)fIndJ_f->Get(Form("lth_sys%s", sv_lbl[i_t].c_str()));
  }    
  fIndJ_f->Close();
  TFile *fIndP_f = new TFile("../Psi2/Simult/Systematics/mainDiffs/files/finalUnc.root");
  TGraphAsymmErrors **graph_lthP_f = new TGraphAsymmErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lthP_f[i_t] = (TGraphAsymmErrors*)fIndP_f->Get(Form("lth_sys%s", sv_lbl[i_t].c_str()));
  }    
  fIndP_f->Close();
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(15, -1, 125, 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);

  int col_j[2] = {kBlue, kRed};
  int col_p[2] = {kViolet+1, kPink+7};
  
  TLegend *leg = new TLegend(0.65, 0.785, 0.95, 0.985);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);

  for(int i_l = 0; i_l < 2; i_l++) {
    
    graph_lthJ[i_l]->SetMarkerSize(.75);
    graph_lthJ[i_l]->SetMarkerStyle(20);
    graph_lthJ[i_l]->SetMarkerColor(col_j[i_l]);
    graph_lthJ[i_l]->SetLineColor(col_j[i_l]);
    graph_lthJ[i_l]->Draw("p same");
    leg->AddEntry(graph_lthJ[i_l], Form("%s J/#psi", nm_lbl[i_l].c_str()), "pl");

    graph_lthP[i_l]->SetMarkerStyle(25);
    graph_lthP[i_l]->SetMarkerSize();
    graph_lthP[i_l]->SetMarkerColor(col_p[i_l]);
    graph_lthP[i_l]->SetLineColor(col_p[i_l]);
    graph_lthP[i_l]->Draw("p same");
    leg->AddEntry(graph_lthP[i_l], Form("%s #psi(2S)", nm_lbl[i_l].c_str()), "pl");
  }
  
  TLine *zero = new TLine(15, 0, 125, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  leg->Draw();
  
  c->SaveAs("par_lth.pdf");
  c->Clear();

  //define multi-error graphs, plot separately
  TGraphMultiErrors **graph_totJ = new TGraphMultiErrors*[n_lbl];
  TGraphMultiErrors **graph_totP = new TGraphMultiErrors*[n_lbl];


  // draw lambda_th(pT)
  TH1F *fl_J = c->DrawFrame(15, -1, 125, 1);
  fl_J->SetXTitle("p_{T} (GeV)");
  fl_J->SetYTitle("#lambda_{#theta}");
  fl_J->GetYaxis()->SetTitleOffset(1.3);
  fl_J->GetYaxis()->SetLabelOffset(0.01);

  TLegend *leg_J = new TLegend(0.65, 0.875, 0.95, 0.975);
  leg_J->SetTextSize(0.03);
  leg_J->SetBorderSize(0);
  leg_J->SetFillColorAlpha(kWhite,0);

  
  for(int i_l = 0; i_l < 2; i_l++) {
    graph_totJ[i_l] = new TGraphMultiErrors(Form("g_J%d", i_l), Form("g_J%d", i_l), graph_lthJ[i_l]->GetN(), graph_lthJ[i_l]->GetX(), graph_lthJ[i_l]->GetY(), graph_lthJ[i_l]->GetEX(), graph_lthJ[i_l]->GetEX(), graph_lthJ[i_l]->GetEY(), graph_lthJ[i_l]->GetEY());
    graph_totJ[i_l]->AddYError(graph_lthJ[i_l]->GetN(), graph_lthJ_f[i_l]->GetEYlow(), graph_lthJ_f[i_l]->GetEYhigh());
    
    graph_totJ[i_l]->SetMarkerSize(.5);
    graph_totJ[i_l]->SetMarkerStyle(20);
    graph_totJ[i_l]->SetMarkerColor(col_j[i_l]);
    graph_totJ[i_l]->SetLineColor(col_j[i_l]);
    graph_totJ[i_l]->GetAttLine(0)->SetLineColor(col_j[i_l]);
    graph_totJ[i_l]->GetAttFill(1)->SetFillStyle(0);
    graph_totJ[i_l]->Draw("p s same; ; 5 s=0.5");
    leg_J->AddEntry(graph_totJ[i_l], Form("%s J/#psi", nm_lbl[i_l].c_str()), "pl");

  }

  leg_J->Draw();
  zero->Draw();
  
  c->SaveAs("par_lth_bJ.pdf");
  c->Clear();  
  c->Destructor();
}
