// code to plot the fit results for both states as lth(pT/M)
// includes both stat and sys uncertainties

void compXi()
{
  // get the fit results
  // get lambda values for each bin
  const int n_lbl = 2;
  string lbl[] = {"PR", "NP"}; 
  string nm_lbl[] = {"prompt", "non-prompt"};
  
  TFile *fIndJ = new TFile("../Jpsi/Simult/Systematics/mainDiffs/files/finalUnc.root");
  TGraphAsymmErrors **graph_lthJ = new TGraphAsymmErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lthJ[i_t] = (TGraphAsymmErrors*)fIndJ->Get(Form("lth_f%s", lbl[i_t].c_str()));
  }    
  fIndJ->Close();
  cout << "got the Jpsi results" << endl;
  TFile *fIndP = new TFile("../Psi2/Simult/Systematics/mainDiffs/files/finalUnc.root");
  TGraphAsymmErrors **graph_lthP = new TGraphAsymmErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lthP[i_t] = (TGraphAsymmErrors*)fIndP->Get(Form("lth_f%s", lbl[i_t].c_str()));
  }    
  fIndP->Close();
  cout << "got the psi(2S) results" << endl;

  // plot in pT/M
  int nJ = graph_lthJ[0]->GetN();
  for(int j = 0; j < 2; j++)
    for(int i = 0; i < nJ; i++) {
      graph_lthJ[j]->GetX()[i]/=3.097;
      graph_lthJ[j]->GetEXhigh()[i]/=3.097;
      graph_lthJ[j]->GetEXlow()[i]/=3.097;
    }
  int nP = graph_lthP[0]->GetN();
  for(int j = 0; j < 2; j++)
    for(int i = 0; i < nP; i++) {
      graph_lthP[j]->GetX()[i]/=3.686;
      graph_lthP[j]->GetEXhigh()[i]/=3.686;
      graph_lthP[j]->GetEXlow()[i]/=3.686;
    }
  
  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);

  // draw lambda_th(pT)
  TH1F *fl = c->DrawFrame(0, -1, 40, 1);
  fl->SetXTitle("p_{T}/M");
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
  
  TLine *zero = new TLine(0, 0, 50, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  
  c->SaveAs("par_lthXi.pdf");
  c->Clear();
  c->Destructor();
}
