void lbdTil()
{
  TFile *fcos = new TFile("../Simult/PR_fit/files/finalFitRes.root");
  TGraphErrors *g_cos = (TGraphErrors*)fcos->Get("graph_lambda_J");
  fcos->Close();
  TFile *fphi = new TFile("../Phi_fit/PR_fit/files/finalFitRes.root");
  TGraphErrors *g_phi = (TGraphErrors*)fphi->Get("graph_lambda_J");
  fphi->Close();

  const int n = g_cos->GetN();
  double vy[n], ey[n];
  for(int i = 0; i < n; i++) {
    double v_c = g_cos->GetY()[i];
    double v_p = g_phi->GetY()[i];
    vy[i] = (v_c+3.*v_p)/(1.-v_p);
    cout << i << " " << v_c << " " << v_p << " " << vy[i] << endl;
    double e_c = g_cos->GetEY()[i];
    double e_p = g_phi->GetEY()[i];
    ey[i] = 1./(1.-v_p)*sqrt(pow(e_c,2)+pow((3.+v_c)*e_p/(1.-v_p),2));
  }
  TGraphErrors *g_t = new TGraphErrors(n, g_cos->GetX(), vy, g_cos->GetEX(), ey);

  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  TH1F *fl = c->DrawFrame(15, -0.5, 100, 0.5);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  
  g_cos->SetMarkerColor(kBlue);  
  g_cos->SetLineColor(kBlue);
  g_cos->Draw("p same");

  g_phi->SetMarkerColor(kRed);  
  g_phi->SetLineColor(kRed);
  g_phi->Draw("p same");

  g_t->SetMarkerColor(kBlack);  
  g_t->SetLineColor(kBlack);
  g_t->Draw("p same");

  TLine *zero = new TLine(15, 0, 100, 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();

  TF1 *flin = new TF1("flin", "[0]+[1]*x", 15, 100);
  flin->SetLineColor(kViolet);
  flin->SetLineStyle(kDashed);
  flin->SetParameters(0.1, 0.1);
  g_cos->Fit(flin, "0");
  g_phi->Fit(flin, "0");
  g_t->Fit(flin, "0");
  
  c->SaveAs("plots/lbd_tilde.pdf");
  c->Destructor();
}
