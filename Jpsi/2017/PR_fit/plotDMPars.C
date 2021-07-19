// macro to draw all the fit parameters
void plotDMPars()
{
  // aux arrays
  const int n_p = 10;

  string parlab[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda", "fBG"};
  string partit[] = {"f", "N_{SR}", "#mu", "#sigma_{1}", "#sigma_{2}", "n", "#alpha", "N_{BG}", "#lambda", "f_{bkg}"};
  string parax[] = {"f (%)", "N_{SR} per 1 GeV", "#mu (MeV)", "#sigma_{1} (MeV)", "#sigma_{2} (MeV)", "n", "#alpha", "N_{BG} per 1 GeV", "#lambda (MeV)", "f_{bkg} (%)"};
  
  double parmin[] = {0,    5e0, 3093.4, 15, 30, 1.0, 2.0, 3e0, 0, 0.};
  double parmax[] = {100., 1e4, 3094.4, 30, 45, 1.4, 2.3, 2e4, 5000, 15.};

  // initialize tgraphs for parameters
  TGraphErrors **g_par = new TGraphErrors*[n_p];
  TFile *fin = new TFile(Form("files/mfit.root"));
  for(int i_p = 0; i_p < n_p; i_p++) {
    fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_p]);
  }
  fin->Close();

  double *xv = g_par[0]->GetX();
  double *xe = g_par[0]->GetEX();
  int n = g_par[0]->GetN();
  // scale all graphs for plotting
  TGraphErrors **g_par_s = new TGraphErrors*[n_p];
  double mult[] = {100., 1., 1e3, 1e3, 1e3, 1., 1., 1., 1e3, 100.};
  for(int i = 0; i < n_p; i++) {
    double *yv = g_par[i]->GetY();
    double *ye = g_par[i]->GetEY();

    for(int j = 0; j < n; j++) {
      if(i == 1 || i == 7) {
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
  double pt_max = xv[n-1]+xe[n-1]+5;
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);

  for(int i_p = 0; i_p < n_p; i_p++) {

    if(i_p == 1 || i_p == 7) c->SetLogy();
    else c->SetLogy(0);
    
    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.5);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("2017 %s", partit[i_p].c_str()));

    // free pars
    if(i_p == 1 || i_p >= 5) {
      g_par_s[i_p]->SetMarkerStyle(20);
      g_par_s[i_p]->SetMarkerSize(.75);
      g_par_s[i_p]->SetLineColor(kBlue);
      g_par_s[i_p]->SetMarkerColor(kBlue);
      g_par_s[i_p]->Draw("p");

      }
    // linear or constant pars 
    else {
      g_par_s[i_p]->SetLineColor(kBlue);
      g_par_s[i_p]->SetFillColorAlpha(kBlue, 0.5);
      g_par_s[i_p]->Draw("ce3");
    }

    c->SaveAs(Form("plots/mass/par_%s.pdf", parlab[i_p].c_str()));
    c->Clear();
  }
  
  c->Destructor();
}
