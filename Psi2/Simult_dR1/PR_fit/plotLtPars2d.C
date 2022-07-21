double sig_m(double f, double sig1, double sig2)
{
  double sig1_s = pow(sig1, 2);
  double sig2_s = pow(sig2, 2);
  f /= 100.;

  return sqrt(f*sig1_s + (1.-f)*sig2_s);
}

double esig_m(double f, double sig1, double sig2, double ef, double esig1, double esig2)
{
  double sig1_s = pow(sig1, 2);
  double sig2_s = pow(sig2, 2);
  f/= 100.;
  ef /= 100.;

  return sqrt(pow((sig1_s-sig2_s)*ef/2, 2) + pow(f * sig1 * esig1, 2) + pow((1.-f) * sig2 * esig2, 2))/sqrt(f*sig1_s + (1.-f)*sig2_s);
}

void plotLtPars2d()
{
  // aux arrays
  const int n_p = 7;

  string parlab[] = {"N_PR", "N_NP", "f", "mu", "sig1", "sig2", "lambda"};
  string partit[] = {"N_{PR}", "N_{NP}", "f", "#mu", "#sigma_{1}", "#sigma_{2}", "t_{NP}"};
  string par_unit[] = {" per 1 GeV", " per 1 GeV", " (%)", " (#mum)", " (#mum)", " (#mum)", " (#mum)"};

  double parmin[] = {2e1, 7e1, 0,  -5, 0,  0,  300};
  double parmax[] = {2e4, 5e4, 100, 5, 20, 40, 400};

  // initialize tgraphs for parameters
  TGraphErrors **g_par = new TGraphErrors*[n_p];
  TFile *fin = new TFile(Form("files/ltfitres2d.root"));
  for(int i_p = 0; i_p < n_p; i_p++) {
    fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_p]);
  }
  fin->Close();

  double *xv = g_par[0]->GetX();
  double *xe = g_par[0]->GetEX();
  int n = g_par[0]->GetN();
  // scale all graphs for plotting
  TGraphErrors **g_par_s = new TGraphErrors*[n_p];
  double mult[] = {1., 1., 100., 1e3, 1e3, 1e3, 1e3};
  for(int i = 0; i < n_p; i++) {
    double *yv = g_par[i]->GetY();
    double *ye = g_par[i]->GetEY();

    for(int j = 0; j < n; j++) {
      if(i < 2) {
	yv[j] /= (2.*xe[j]);
	ye[j] /= (2.*xe[j]);
      }
      else {
	yv[j] *= mult[i];
	ye[j] *= mult[i];
      }
    }
    g_par_s[i] = new TGraphErrors(n, xv, yv, xe, ye);
  }

  double pt_min = xv[0]-xe[0]-5;
  double pt_max = xv[n-1]+xe[n-1];
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);

  for(int i = 0; i < n_p; i++) {

    if(i < 2) c->SetLogy();
    else c->SetLogy(0);

    TH1F *fp = c->DrawFrame(pt_min, parmin[i], pt_max, parmax[i]);
    fp->SetXTitle("p_{T} (GeV)");
    fp->SetYTitle(Form("%s%s", partit[i].c_str(), par_unit[i].c_str()));
    fp->GetYaxis()->SetTitleOffset(1.5);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("Run 2 %s", partit[i].c_str()));

    // constant pars
    if( i == 2 || i == 3) {
      g_par_s[i]->SetLineColor(kBlack);
      g_par_s[i]->SetFillColorAlpha(kBlack, 0.5);
      g_par_s[i]->Draw("ce3");
    }

    else {
      g_par_s[i]->SetMarkerStyle(20);
      g_par_s[i]->SetMarkerSize(.75);
      g_par_s[i]->SetMarkerColor(kBlack);
      g_par_s[i]->SetLineColor(kBlack);
      g_par_s[i]->Draw("psame");
    }
  
    c->SaveAs(Form("plots/lifetime2d/par_%s.pdf", parlab[i].c_str()));
    c->Clear();
  }
    
  // do the sigma plots to compare with rms from previous ANs
  double sig_avg[n], esig_avg[n], rms[n], erms[n], sig_rat[n], esig_rat[n];

  double *fv = g_par[2]->GetY();
  double *s1v = g_par[4]->GetY();
  double *s2v = g_par[5]->GetY();
  double *fe = g_par[2]->GetEY();
  double *s1e = g_par[4]->GetEY();
  double *s2e = g_par[5]->GetEY();

  for(int i_p = 0; i_p < n; i_p++){
    sig_avg[i_p] = sig_m(fv[i_p], s1v[i_p], s2v[i_p]);
    esig_avg[i_p] = esig_m(fv[i_p], s1v[i_p], s2v[i_p], fe[i_p], s1e[i_p], s2e[i_p]);

    rms[i_p] = sig_avg[i_p] * xv[i_p]/3.097;
    erms[i_p] = esig_avg[i_p] * xv[i_p]/3.097;

    sig_rat[i_p] = s2v[i_p]/s1v[i_p];
    esig_rat[i_p] = sig_rat[i_p] * sqrt(pow(s1e[i_p]/s1v[i_p], 2) + pow(s2e[i_p]/s2v[i_p], 2));
  }

  TGraphErrors *g_sig = new TGraphErrors(n, xv, sig_avg, xe, esig_avg);

  TH1F *fsig = c->DrawFrame(pt_min, 0., pt_max, 20);
  fsig->SetXTitle("p_{T} (GeV)");
  fsig->SetYTitle("#sigma (#mum)");
  fsig->GetYaxis()->SetTitleOffset(1.3);
  fsig->GetYaxis()->SetLabelOffset(0.01);
  fsig->SetTitle("Run 2 #sigma_{m}");
  
  g_sig->SetMarkerStyle(20);
  g_sig->SetMarkerSize(.75);
  g_sig->SetMarkerColor(kBlack);
  g_sig->SetLineColor(kBlack);
  g_sig->Draw("psame");

  c->SaveAs(Form("plots/lifetime2d/par_sigma.pdf"));
  c->Clear();

  TGraphErrors *g_rms = new TGraphErrors(n, xv, rms, xe, erms);

  TH1F *frms = c->DrawFrame(pt_min, 0., pt_max, 500);
  frms->SetXTitle("p_{T} (GeV)");
  frms->SetYTitle("rms (Lxy) (#mum)");
  frms->GetYaxis()->SetTitleOffset(1.3);
  frms->GetYaxis()->SetLabelOffset(0.01);
  frms->SetTitle("Run 2 rms (Lxy)");
  
  g_rms->SetMarkerStyle(20);
  g_rms->SetMarkerSize(.75);
  g_rms->SetMarkerColor(kBlack);
  g_rms->SetLineColor(kBlack);
  g_rms->Draw("psame");
  
  c->SaveAs(Form("plots/lifetime2d/par_rms.pdf"));
  c->Clear();

  TGraphErrors *g_rat = new TGraphErrors(n, xv, sig_rat, xe, esig_rat);

  TH1F *frat = c->DrawFrame(pt_min, 1., pt_max, 3.);
  frat->SetXTitle("p_{T} (GeV)");
  frat->SetYTitle("#sigma_{2}/#sigma_{1}");
  frat->GetYaxis()->SetTitleOffset(1.3);
  frat->GetYaxis()->SetLabelOffset(0.01);
  frat->SetTitle("Run 2 #sigma_{2}/#sigma_{1}");
  
  g_rat->SetMarkerStyle(20);
  g_rat->SetMarkerSize(.75);
  g_rat->SetMarkerColor(kBlack);
  g_rat->SetLineColor(kBlack);
  g_rat->Draw("psame");

  c->SaveAs(Form("plots/lifetime2d/par_sigRat.pdf"));
  c->Clear();

  c->Destructor();

}
