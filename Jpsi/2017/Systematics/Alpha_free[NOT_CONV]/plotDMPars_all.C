// macro to compare fit parameters w/ alpha free or fixed
void plotDMPars_all()
{
  // aux arrays
  int pc[] = {kBlue, kGreen};
  int sigc[] = {kBlack, kViolet};
  const int n_p = 12, n_m = 2;
  string modn[] = {"../../PR_fit/", ""};
  
  string parlab[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda", "fBG", "fG", "sigG"};
  string partit[] = {"f", "N_{SR}", "#mu", "#sigma", "#sigma_{2}", "n", "#alpha", "N_{BG}", "#lambda", "f_{bkg}", "f_{G}", "#sigma_{G}"};
  string parax[] = {"f (%)", "N_{SR} per 1 GeV", "#mu (MeV)", "#sigma_{1} (MeV)", "#sigma_{2} (MeV)", "n", "#alpha", "N_{BG} per 1 GeV", "#lambda (GeV)", "f_{bkg} (%)", "f_{G} (%)", "#sigma_{G} (MeV)"};
  
  double parmin[] = {0,    1e1, 3090, 0,   32, 1.0, 2.0,  3e0, 0,  0.,  0,   0};
  double parmax[] = {100., 1e4, 3100, 100, 46, 1.4, 2.25, 3e4, 60, 15., 100, 100};
 
  // initialize tgraphs for parameters
  TGraphErrors ***g_par = new TGraphErrors**[n_m];
  for(int i_m = 0; i_m < n_m; i_m++) {
    g_par[i_m] = new TGraphErrors*[n_p];
    TFile *fin = new TFile(Form("%sfiles/mfit.root", modn[i_m].c_str()));
    for(int i_p = 0; i_p < n_p; i_p++) {
      fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_m][i_p]);
    }
    fin->Close();
  }
  
  double *xv = g_par[0][0]->GetX();
  double *xe = g_par[0][0]->GetEX();
  int n = g_par[0][0]->GetN();
  // scale all graphs for plotting
  TGraphErrors ***g_par_s = new TGraphErrors**[n_m];
  double mult[] = {100., 1., 1e3, 1e3, 1e3, 1., 1., 1., 1., 100., 100., 1e3};
  for(int i_m = 0; i_m < n_m; i_m++) {
    g_par_s[i_m] = new TGraphErrors*[n_p];
    for(int i = 0; i < n_p; i++) {
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
  
  double pt_min = xv[0]-xe[0]-5;
  double pt_max = xv[n-1]+xe[n-1]+5;
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);

  for(int i_p = 0; i_p < 10; i_p++) {

    if(i_p == 4) continue;

    if(i_p == 1 || i_p == 7) c->SetLogy();
    else c->SetLogy(0);
    
    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.5);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("2018 %s", partit[i_p].c_str()));

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);

    // drawing both models
    for(int i_m = 0; i_m < n_m; i_m++) {
      // free pars
      if(i_p == 1 || i_p == 7 || i_p == 8 || i_p == 9) {
	g_par_s[i_m][i_p]->SetMarkerStyle(20);
	g_par_s[i_m][i_p]->SetMarkerSize(.75);
	g_par_s[i_m][i_p]->SetLineColor(pc[i_m]);
	g_par_s[i_m][i_p]->SetMarkerColor(pc[i_m]);
	g_par_s[i_m][i_p]->Draw("p");	
      }
      // linear or constant pars 
      else if(i_p != 3) {
	g_par_s[i_m][i_p]->SetMarkerStyle(20);
	g_par_s[i_m][i_p]->SetMarkerSize(.75);
	g_par_s[i_m][i_p]->SetLineColor(pc[i_m]);
	g_par_s[i_m][i_p]->SetMarkerColor(pc[i_m]);
	g_par_s[i_m][i_p]->SetFillColorAlpha(pc[i_m], 0.5);
	g_par_s[i_m][i_p]->Draw("ce3");
      }

      // SPECIAL CASES - COMBINED PLOTS
      // if we're plotting f, add fG
      if( i_p == 0) {      
	g_par_s[i_m][10]->SetMarkerStyle(22);
	g_par_s[i_m][10]->SetLineColor(kBlack);
	g_par_s[i_m][10]->SetMarkerColor(kBlack);
	g_par_s[i_m][10]->SetFillColorAlpha(kBlack, 0.5);
	g_par_s[i_m][10]->Draw("ce3");
      }
 
      // if we're plotting par sig1, add sig2 and sigG
      else if( i_p == 3) {
	g_par_s[i_m][i_p]->SetMarkerStyle(20);
	g_par_s[i_m][i_p]->SetMarkerSize(.75);
	g_par_s[i_m][i_p]->SetMarkerColor(pc[i_m]);
	g_par_s[i_m][i_p]->SetLineColor(pc[i_m]);
	g_par_s[i_m][i_p]->SetFillColorAlpha(pc[i_m], 0.5);
	g_par_s[i_m][i_p]->Draw("ce3");

	g_par_s[i_m][4]->SetMarkerStyle(22);
	g_par_s[i_m][4]->SetMarkerColor(sigc[i_m]);
	g_par_s[i_m][4]->SetLineColor(sigc[i_m]);
	g_par_s[i_m][4]->SetFillColorAlpha(sigc[i_m], 0.5);
	g_par_s[i_m][4]->Draw("ce3");

	g_par_s[i_m][11]->SetMarkerStyle(29);
	g_par_s[i_m][11]->SetMarkerColor(kRed);
	g_par_s[i_m][11]->SetLineColor(kRed);
	g_par_s[i_m][11]->SetFillColorAlpha(kRed, 0.5);
	g_par_s[i_m][11]->Draw("ce3");
      }
    }

    if(i_p == 0) {
      leg->AddEntry(g_par_s[0][0], "f_{CB1}", "l");
      leg->AddEntry(g_par_s[0][10], "f_{G}", "l");
      leg->Draw();
    }
    else if(i_p == 3) {
      leg->AddEntry(g_par_s[0][3], "#sigma_{1}", "l");
      leg->AddEntry(g_par_s[0][4], "#sigma_{2}", "l");
      leg->AddEntry(g_par_s[0][11], "#sigma_{G}", "l");
      leg->Draw();
    }
    
    if(i_p == 2) {
      TLine *pdg = new TLine(pt_min, 3096.9, pt_max, 3096.9);
      pdg->SetLineStyle(kDashed);
      pdg->Draw();
    }

    int isLog = 0;
    if(i_p == 1 || i_p == 7 ) isLog = 1;

    c->SaveAs(Form("plots/mass/parA_%s.pdf", parlab[i_p].c_str()));
    c->Clear();
  }
  
  c->Destructor();

}
