// macro to compare fit parameters w/ and w/o Gaussian
// get relative position on an axis (pi, pf)
void plotDMPars_NP()
{
  // aux arrays
  int pc[] = {kBlack, kBlue, kViolet};
  const int n_p = 12, n_m = 1;
  string modn[] = {""};
  string legn[] = {"no G", "with G"};

  string parlab[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda", "fBG", "fG", "sigG"};
  string parsave[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "tbkg", "fBG", "fG", "sigG"};
  string partit[] = {"f", "N_{SR}", "#mu", "#sigma", "#sigma_{2}", "n", "#alpha", "N_{BG}", "t_{Bg}", "f_{Bg}", "f_{G}", "#sigma_{G}"};
  string parax[] = {"f (%)", "N_{SR} per 1 GeV", "#mu (MeV)", "#sigma (MeV)", "#sigma_{2} (MeV)", "n", "#alpha", "N_{BG} per 1 GeV", "t_{Bg} (GeV)", "f_{Bg} (%)", "f_{G} (%)", "#sigma_{G} (MeV)"};
  
  double parmin[] = {0,    6e0, 3090, 0,   32, 2.0, 1.0, 1e4, 0, 0.,  0,   0};
  double parmax[] = {100., 2e4, 3100, 200, 46, 3.0, 2.3, 8e6, 1, 15., 100, 100};
 
  // initialize tgraphs for parameters
  TGraphErrors ***g_par = new TGraphErrors**[n_m];
  for(int i_m = 0; i_m < n_m; i_m++) {
    g_par[i_m] = new TGraphErrors*[n_p];
    TFile *fin = new TFile(Form("files/mfit_NP%s.root", modn[i_m].c_str()));
    int n_pv = n_p;
    for(int i_p = 0; i_p < n_pv; i_p++) {
      fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_m][i_p]);
    }
    fin->Close();
  }
  
  // scale all graphs for plotting
  TGraphErrors ***g_par_s = new TGraphErrors**[n_m];
  double mult[] = {100., 1., 1e3, 1e3, 1e3, 1., 1., 1., 1., 100., 100., 1e3};
  double pt_min, pt_max;
  for(int i_m = 0; i_m < n_m; i_m++) {
    double *xv = g_par[i_m][0]->GetX();
    double *xe = g_par[i_m][0]->GetEX();
    int n = g_par[i_m][0]->GetN();
    pt_min = xv[0]-xe[0]-5;
    pt_max = xv[n-1]+xe[n-1]+5;
    g_par_s[i_m] = new TGraphErrors*[n_p];
    //int n_pv = i_m == 0 ? 10 : 12;
    int n_pv = n_p;
    for(int i = 0; i < n_pv; i++) {
      double *yv = g_par[i_m][i]->GetY();
      double *ye = g_par[i_m][i]->GetEY();

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
      g_par_s[i_m][i] = new TGraphErrors(n, xv, yv, xe, ye);
    }
  }

  // for the NB estimate, get histogram maximum and integral as well
  TH2D *h_d2d = new TH2D();
  TFile *finh = new TFile("files/mStore.root");
  finh->GetObject("mH_NP", h_d2d);
  h_d2d->SetDirectory(0);
  finh->Close();

  const int nPtBins = h_d2d->GetNbinsY();
  const double *ptBins = h_d2d->GetYaxis()->GetXbins()->GetArray();
  
  // Make 1d histos
  TH1D **h_d1d = new TH1D*[nPtBins];
  TH1D *h_di = new TH1D("h_di", "name", nPtBins, ptBins);
  for(int i = 0; i < nPtBins; i++) {
    h_d1d[i] = h_d2d->ProjectionX(Form("mH%.0f", ptBins[i]), i+1, i+1);
    double m_c = h_d1d[i]->Integral()/(ptBins[i+1]-ptBins[i]);
    //m_c /= (1.5*i+1);
    m_c*=4;
    h_di->SetBinContent(i+1, m_c);
  }

  
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.12);

  for(int i_p = 0; i_p < 10; i_p++) {

    if(i_p == 4) continue;

    if(i_p == 1 || i_p == 7) c->SetLogy();
    else c->SetLogy(0);
    
    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.8);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("Run 2 %s", partit[i_p].c_str()));

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);


    // drawing base+G model
    // free pars - NS, NB, lambda
    if(i_p == 1 || i_p > 6) {
      for(int i_n = 0; i_n < n_m; i_n++) {
	g_par_s[i_n][i_p]->SetMarkerStyle(20);
	g_par_s[i_n][i_p]->SetMarkerSize(.75);
	g_par_s[i_n][i_p]->SetLineColor(pc[i_n]);
	g_par_s[i_n][i_p]->SetMarkerColor(pc[i_n]);
	g_par_s[i_n][i_p]->Draw("p");
      }

      //if(i_p == 7) h_di->Draw("histo same");

    }
    // linear or constant pars 
    else if(i_p != 3) {
      for(int i_n = 0; i_n < n_m; i_n++) {
	g_par_s[i_n][i_p]->SetMarkerStyle(20);
	g_par_s[i_n][i_p]->SetMarkerSize(.75);
	g_par_s[i_n][i_p]->SetLineColor(pc[i_n]);
	g_par_s[i_n][i_p]->SetMarkerColor(pc[i_n]);
	g_par_s[i_n][i_p]->SetFillColorAlpha(pc[i_n], 0.5);
	g_par_s[i_n][i_p]->Draw("pce3");
      }
    }

    // SPECIAL CASES - COMBINED PLOTS
    // if we're plotting f, add fG
    if( i_p == 0) {      
      for(int i_n = 0; i_n < n_m; i_n++) {
	g_par_s[i_n][10]->SetMarkerStyle(22);
	g_par_s[i_n][10]->SetLineColor(pc[i_n]);
	g_par_s[i_n][10]->SetMarkerColor(pc[i_n]);
	g_par_s[i_n][10]->SetFillColorAlpha(pc[i_n], 0.5);
	g_par_s[i_n][10]->Draw("pce3");
      }
      
      leg->AddEntry(g_par_s[0][0], "f_{CB1}", "pl");
      leg->AddEntry(g_par_s[0][10], "f_{G}", "pl");
      leg->Draw();
    }
 
    // if we're plotting par sig1, add sig2 and sigG
    else if( i_p == 3) {
      for(int i_n = 0; i_n < n_m; i_n++) {
	g_par_s[i_n][3]->SetMarkerStyle(20);
	g_par_s[i_n][3]->SetMarkerSize(.75);
	g_par_s[i_n][3]->SetMarkerColor(pc[i_n]);
	g_par_s[i_n][3]->SetLineColor(pc[i_n]);
	g_par_s[i_n][3]->SetFillColorAlpha(pc[i_n], 0.5);
	g_par_s[i_n][3]->Draw("pce3");

	g_par_s[i_n][4]->SetMarkerStyle(22);
	g_par_s[i_n][4]->SetMarkerColor(pc[i_n]);
	g_par_s[i_n][4]->SetLineColor(pc[i_n]);
	g_par_s[i_n][4]->SetFillColorAlpha(pc[i_n], 0.5);
	g_par_s[i_n][4]->Draw("pce3");

	g_par_s[i_n][11]->SetMarkerStyle(29);
	g_par_s[i_n][11]->SetMarkerColor(pc[i_n]);
	g_par_s[i_n][11]->SetLineColor(pc[i_n]);
	g_par_s[i_n][11]->SetFillColorAlpha(pc[i_n], 0.5);
	g_par_s[i_n][11]->Draw("pce3");
      }
      
      leg->AddEntry(g_par_s[0][3], "#sigma_{1}", "pl");
      leg->AddEntry(g_par_s[0][4], "#sigma_{2}", "pl");
      leg->AddEntry(g_par_s[0][11], "#sigma_{G}", "pl");
      leg->Draw();
    }

    
    int isLog = 0;
    if(i_p == 1 || i_p == 7 ) isLog = 1;

    c->SaveAs(Form("plots/massNP/par_%s.pdf", parsave[i_p].c_str()));
    c->Clear();
  }
  
  c->Destructor();

}
