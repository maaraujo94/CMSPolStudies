// code to plot the fit results

void compStates()
{
  // get the fit results
  // get lambda values for each bin
  const int n_lbl = 3;
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
  for(int i_l = 0; i_l < 2; i_l++) {
    
    graph_lthJ[i_l]->SetMarkerSize(.75);
    graph_lthJ[i_l]->SetMarkerStyle(20);
    graph_lthJ[i_l]->SetMarkerColor(col_j[i_l]);
    graph_lthJ[i_l]->SetLineColor(col_j[i_l]);
    graph_lthJ[i_l]->Draw("p same");
    
    graph_lthP[i_l]->SetMarkerStyle(25);
    graph_lthP[i_l]->SetMarkerSize();
    graph_lthP[i_l]->SetMarkerColor(col_p[i_l]);
    graph_lthP[i_l]->SetLineColor(col_p[i_l]);
    graph_lthP[i_l]->Draw("p same");
  }
  
  TLine *zero = new TLine(15, 0, 125, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TLegend *leg = new TLegend(0.82, 0.75, 0.97, 0.9);
  leg->SetTextSize(0.03);
  //leg->AddEntry(graph_lthJ[i_l], "J/#psi", "pl");
  //leg->AddEntry(graph_lthP[i_l], "#psi(2S)", "pl");
  //leg->Draw();
  
  c->SaveAs("par_lth.pdf");
  c->Clear();  
  c->Destructor();
}
