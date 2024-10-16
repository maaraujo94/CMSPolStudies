// get lambda_tilde from 3 diff sources for NP case
// 1d HX, 2d HX, 2d CS (w/ and w/o ltp)
void lbdTil_NP()
{
  // get 1d HX results
  TFile *fcos = new TFile("../../Simult/NP_fit/files/finalFitRes.root");
  TGraphErrors *g_cos = (TGraphErrors*)fcos->Get("graph_lambda_NPc");
  fcos->Close();
  // gives beta = f(lth,lphi)
  TFile *fphi = new TFile("../../Phi_fit/NP_fit/files/finalFitRes.root");
  TGraphErrors *g_phi = (TGraphErrors*)fphi->Get("graph_lambda_NPc");
  fphi->Close();

  // now get the 2d HX
  TFile *fInd_HX = new TFile("../NP_fit/files/finalFitRes.root");
  TGraphErrors **graph_lHX = new TGraphErrors*[3];
  string lbl[3] = {"theta", "phi", "thph"};
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lHX[i_t] = (TGraphErrors*)fInd_HX->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
    graph_lHX[i_t]->SetName(Form("graph_lHX_%s", lbl[i_t].c_str()));
  }    
  fInd_HX->Close();

  // now w/o ltp
  TFile *fInd_HXz = new TFile("../NP_fit/files/finalFitRes_ltpz.root");
  TGraphErrors **graph_lHXz = new TGraphErrors*[3];
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lHXz[i_t] = (TGraphErrors*)fInd_HXz->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
    graph_lHXz[i_t]->SetName(Form("graph_lHX_%s", lbl[i_t].c_str()));
  }    
  fInd_HXz->Close();

  // now get the 2d CS
  TFile *fInd_CS = new TFile("../../Simult2d_CS/NP_fit/files/finalFitRes.root");
  TGraphErrors **graph_lCS = new TGraphErrors*[3];
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lCS[i_t] = (TGraphErrors*)fInd_CS->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
    graph_lCS[i_t]->SetName(Form("graph_lCS_%s", lbl[i_t].c_str()));
  }    
  fInd_CS->Close();

  // now w/o ltp
  TFile *fInd_CSz = new TFile("../../Simult2d_CS/NP_fit/files/finalFitRes_ltpz.root");
  TGraphErrors **graph_lCSz = new TGraphErrors*[3];
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lCSz[i_t] = (TGraphErrors*)fInd_CSz->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
    graph_lCSz[i_t]->SetName(Form("graph_lCS_%s", lbl[i_t].c_str()));
  }    
  fInd_CSz->Close();

  // calculate ltil for all cases
  const int n = g_cos->GetN();
  double v_phi[n], e_phi[n];
  double vy[5][n], ey[5][n];
  for(int i = 0; i < n; i++) {
    // 1d HX
    double v_c = g_cos->GetY()[i];
    double v_p = 0.5*g_phi->GetY()[i]*(3+v_c);
    vy[0][i] = (v_c+3.*v_p)/(1.-v_p);
    double e_c = g_cos->GetEY()[i];
    double e_p = g_phi->GetEY()[i];
    e_phi[i] = sqrt(pow((3.+v_c)/2., 2)*pow(e_p,2) + pow(g_phi->GetY()[i]/2.,2)*pow(e_c,2));    
    ey[0][i] = 1./(1.-v_p)*sqrt(pow(e_c,2)+pow((3.+v_c)*e_phi[i]/(1.-v_p),2));

    // 2d HX
    v_c = graph_lHX[0]->GetY()[i];
    v_p = graph_lHX[1]->GetY()[i];
    vy[1][i] = (v_c+3.*v_p)/(1.-v_p);
    e_c = graph_lHX[0]->GetEY()[i];
    e_p = graph_lHX[1]->GetEY()[i];
    ey[1][i] = 1./(1.-v_p)*sqrt(pow(e_c,2)+pow((3.+v_c)*e_p/(1.-v_p),2));

    // w/o ltp
    v_c = graph_lHXz[0]->GetY()[i];
    v_p = graph_lHXz[1]->GetY()[i];
    vy[2][i] = (v_c+3.*v_p)/(1.-v_p);
    e_c = graph_lHXz[0]->GetEY()[i];
    e_p = graph_lHXz[1]->GetEY()[i];
    ey[2][i] = 1./(1.-v_p)*sqrt(pow(e_c,2)+pow((3.+v_c)*e_p/(1.-v_p),2));

    // 2d CS
    v_c = graph_lCS[0]->GetY()[i];
    v_p = graph_lCS[1]->GetY()[i];
    vy[3][i] = (v_c+3.*v_p)/(1.-v_p);
    e_c = graph_lCS[0]->GetEY()[i];
    e_p = graph_lCS[1]->GetEY()[i];
    ey[3][i] = 1./(1.-v_p)*sqrt(pow(e_c,2)+pow((3.+v_c)*e_p/(1.-v_p),2));
    
    // w/o ltp
    v_c = graph_lCSz[0]->GetY()[i];
    v_p = graph_lCSz[1]->GetY()[i];
    vy[4][i] = (v_c+3.*v_p)/(1.-v_p);
    e_c = graph_lCSz[0]->GetEY()[i];
    e_p = graph_lCSz[1]->GetEY()[i];
    ey[4][i] = 1./(1.-v_p)*sqrt(pow(e_c,2)+pow((3.+v_c)*e_p/(1.-v_p),2));
  }
  
  TGraphErrors **g_t = new TGraphErrors*[5];
  for(int i = 0; i < 5; i++) {
    g_t[i] = new TGraphErrors(n, g_cos->GetX(), vy[i], g_cos->GetEX(), ey[i]);
  }

  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  TH1F *fl = c->DrawFrame(20, -1, 120, 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#tilde{#lambda}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);

  int col[] = {kBlack, kBlack, kBlack, kRed, kRed};
  int mkr[] = {25, 20, 20, 20, 20};
  int stl[] = {kDashed, kSolid, kSolid, kSolid, kSolid};
  for(int i = 0; i < 5; i++) {
    if(i%2!=0 || i==0) {
      g_t[i]->SetMarkerStyle(mkr[i]);
      g_t[i]->SetLineStyle(stl[i]);
      g_t[i]->SetMarkerColor(col[i]);
      g_t[i]->SetLineColor(col[i]);
      g_t[i]->Draw("p same");
    }
  }
  
  TLine *zero = new TLine(20, 0, 120, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  string lbl_leg[] = {"HX 1d", "HX 2d", "HX 2d (w/o #lambda_{#theta#phi})", "CS 2d", "CS 2d (w/o #lambda_{#theta#phi})"};
  
  TLegend *leg = new TLegend(0.15,0.15,0.45,0.3);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  leg->SetTextSize(0.03);
  for(int i = 0; i < 5; i++) {
    if(i%2!=0 || i == 0)
      leg->AddEntry(g_t[i], lbl_leg[i].c_str(), "pl");
  }
  leg->Draw();
  
  c->SaveAs("lNP_tilde.pdf");
  c->Clear();

  // just the w/o ltp ones
  TH1F *ftp = c->DrawFrame(20, -1, 120, 1);
  ftp->SetXTitle("p_{T} (GeV)");
  ftp->SetYTitle("#tilde{#lambda}");
  ftp->GetYaxis()->SetTitleOffset(1.3);
  ftp->GetYaxis()->SetLabelOffset(0.01);

  for(int i = 0; i < 5; i++) {
    if(i%2==0) {
      g_t[i]->SetMarkerStyle(mkr[i]);
      g_t[i]->SetLineStyle(stl[i]);
      g_t[i]->SetMarkerColor(col[i]);
      g_t[i]->SetLineColor(col[i]);
      g_t[i]->Draw("p same");
    }
  }
  
  zero->Draw();

  TLegend *legtp = new TLegend(0.15,0.15,0.45,0.3);
  legtp->SetBorderSize(0);
  legtp->SetFillColorAlpha(kWhite,0);
  legtp->SetTextSize(0.03);
  for(int i = 0; i < 5; i++) {
    if(i%2==0)
      legtp->AddEntry(g_t[i], lbl_leg[i].c_str(), "pl");
  }
  legtp->Draw();
  
  c->SaveAs("lNP_tilde_ltpz.pdf");
  c->Clear();
  c->Destructor();
}
