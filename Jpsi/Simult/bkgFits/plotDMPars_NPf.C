// macro to compare fit parameters w/ and w/o Gaussian
// get relative position on an axis (pi, pf)
void plotDMPars_NPf()
{
  // aux arrays
  int pc[] = {kBlack, kBlue, kViolet};
  const int n_p = 12, n_m = 1;
  string modn[] = {"2con"};
  string legn[] = {"no G", "with G"};

  string parlab[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda", "fBG", "fG", "sigG"};
  string partit[] = {"f", "N_{SR}", "#mu", "#sigma", "#sigma_{2}", "n", "#alpha", "N_{BG}", "t", "f_{bkg}", "f_{G}", "#sigma_{G}"};
  string parax[] = {"f (%)", "N_{SR} per 1 GeV", "#mu (MeV)", "#sigma (MeV)", "#sigma_{2} (MeV)", "n", "#alpha", "N_{BG} per 1 GeV", "t (GeV)", "f_{bkg} (%)", "f_{G} (%)", "#sigma_{G} (MeV)"};
  
  double parmin[] = {0,    6e0, 3090, 0,   32, 1.0, 1.0, 8e3, 0, 0.,  0,   0};
  double parmax[] = {100., 2e4, 3100, 100, 46, 1.4, 2.3, 4e6, 1, 15., 100, 100};
 
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
    // free pars - NS, NB, alpha, lambda, fBG
    if(i_p == 1 || i_p > 6) {
      for(int i_n = 0; i_n < n_m; i_n++) {
	g_par_s[i_n][i_p]->SetMarkerStyle(20);
	g_par_s[i_n][i_p]->SetMarkerSize(.75);
	g_par_s[i_n][i_p]->SetLineColor(pc[i_n]);
	g_par_s[i_n][i_p]->SetMarkerColor(pc[i_n]);
	g_par_s[i_n][i_p]->Draw("p");
      }
    }
    else if(i_p == 6) {
      g_par_s[0][i_p]->SetMarkerStyle(20);
      g_par_s[0][i_p]->SetMarkerSize(.75);
      g_par_s[0][i_p]->SetLineColor(pc[0]);
      g_par_s[0][i_p]->SetMarkerColor(pc[0]);
      g_par_s[0][i_p]->Draw("p");

      for(int i_n = 1; i_n < n_m; i_n++) {
	g_par_s[i_n][i_p]->SetMarkerStyle(20);
	g_par_s[i_n][i_p]->SetMarkerSize(.75);
	g_par_s[i_n][i_p]->SetLineColor(pc[i_n]);
	g_par_s[i_n][i_p]->SetMarkerColor(pc[i_n]);
	g_par_s[i_n][i_p]->SetFillColorAlpha(pc[i_n], 0.5);
	g_par_s[i_n][i_p]->Draw("pce3");
      }
      
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
	g_par_s[i_n][10]->Draw("p");
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
	g_par_s[i_n][11]->Draw("p");
      }
      
      leg->AddEntry(g_par_s[0][3], "#sigma_{1}", "pl");
      leg->AddEntry(g_par_s[0][4], "#sigma_{2}", "pl");
      leg->AddEntry(g_par_s[0][11], "#sigma_{G}", "pl");
      leg->Draw();
    }

    
    int isLog = 0;
    if(i_p == 1 || i_p == 7 ) isLog = 1;

    c->SaveAs(Form("plots/massNP/parf_%s.pdf", parlab[i_p].c_str()));
    c->Clear();
  }
  
  c->Destructor();

}
