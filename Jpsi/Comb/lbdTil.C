void lbdTil()
{
  TFile *fcos = new TFile("../Simult/PR_fit/files/finalFitRes.root");
  TGraphErrors *g_cos = (TGraphErrors*)fcos->Get("graph_lambda_J");
  fcos->Close();
  // gives beta = f(lth,lphi)
  TFile *fphi = new TFile("../Phi_fit/PR_fit/files/finalFitRes.root");
  TGraphErrors *g_phi = (TGraphErrors*)fphi->Get("graph_lambda_J");
  fphi->Close();

  const int n = g_cos->GetN();
  double v_phi[n], e_phi[n];
  double vy[n], ey[n];
  for(int i = 0; i < n; i++) {
    double v_c = g_cos->GetY()[i];
    double v_p = 0.5*g_phi->GetY()[i]*(3+v_c);
    v_phi[i] = v_p;
    vy[i] = (v_c+3.*v_p)/(1.-v_p);
    cout << i << " " << v_c << " " << v_p << " " << vy[i] << endl;
    double e_c = g_cos->GetEY()[i];
    double e_p = g_phi->GetEY()[i];
    e_phi[i] = sqrt(pow((3.+v_c)/2., 2)*pow(e_p,2) + pow(g_phi->GetY()[i]/2.,2)*pow(e_c,2));    
    ey[i] = 1./(1.-v_p)*sqrt(pow(e_c,2)+pow((3.+v_c)*e_phi[i]/(1.-v_p),2));
  }
  TGraphErrors *g_t = new TGraphErrors(n, g_cos->GetX(), vy, g_cos->GetEX(), ey);
  TGraphErrors *g_p = new TGraphErrors(n, g_cos->GetX(), v_phi, g_cos->GetEX(), e_phi);

  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  TH1F *fl = c->DrawFrame(20, -1, 120, 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  
  g_cos->SetMarkerColor(kBlue);  
  g_cos->SetLineColor(kBlue);
  g_cos->SetLineStyle(kDashed);
  g_cos->Draw("p same");

  g_p->SetMarkerColor(kRed);  
  g_p->SetLineColor(kRed);
  g_p->SetLineStyle(kDashed);
  g_p->Draw("p same");

  g_t->SetMarkerColor(kBlack);  
  g_t->SetLineColor(kBlack);
  g_t->Draw("p same");

  TLine *zero = new TLine(20, 0, 120, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TF1 *fth = new TF1("fth", "[0]+[1]*x", 25, 120);
  fth->SetLineColor(kBlue);
  fth->SetLineStyle(kDashed);
  fth->SetParameters(0.1, 0.1);
  g_cos->Fit(fth);

  TF1 *ftil = new TF1("ftil", "[0]+[1]*x", 25, 120);
  ftil->SetLineColor(kBlack);
  ftil->SetLineStyle(kDashed);
  ftil->SetParameters(0.1,0.1);
  g_t->Fit(ftil);

  TLegend *leg = new TLegend(0.8,0.75,1.1,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  leg->SetTextSize(0.03);
  leg->AddEntry(g_cos, "#lambda_{#theta}", "pl");
  leg->AddEntry(g_p, "#lambda_{#phi}", "pl");
  leg->AddEntry(g_t, "#tilde{#lambda}", "pl");
  leg->Draw();
  
  c->SaveAs("plots/lbd_tilde.pdf");
  c->Destructor();
}
