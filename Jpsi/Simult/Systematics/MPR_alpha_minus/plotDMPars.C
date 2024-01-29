// macro to compare fit parameters w/ and w/o Gaussian
// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

void plotDMPars()
{
  // aux arrays
  int pc[] = {kBlack, kBlue, kViolet};
  int mc[] = {20, 25, 34};
  const int n_p = 12, n_m = 2;
  string modn[] = {"bkgFits", "Systematics/MPR_alpha_minus"};
  string legn[] = {"with G"};

  string parlab[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda", "fG", "sigG", "fBG"};
  string parsave[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "tbkg", "fG", "sigG", "fBG"};
  string partit[] = {"f_{CB_{1}}", "N_{SR}", "#mu", "#sigma", "#sigma_{2}", "n", "#alpha", "A_{BG}", "t_{Bg}", "f_{G}", "#sigma_{G}", "f_{Bg}"};
  string parax[] = {"f_{CB_{1}} (%)", "N_{SR} per 1 GeV", "#mu (MeV)", "#sigma (MeV)", "#sigma_{2} (MeV)", "n", "#alpha", "A_{BG} per 1 GeV", "t_{Bg} (GeV)", "f_{G} (%)", "#sigma_{G} (MeV)", "f_{Bg} (%)"};
  
  double parmin[] = {0,    3e0, 3090, 0,  32, 2.0, 1.3, 2e0, 0, 0,   0,  0.};
  double parmax[] = {100., 1e4, 3100, 80, 46, 3.0, 2.3, 2e6, 4, 100, 80, 15.};
 
  // initialize tgraphs for parameters
  TGraphErrors ***g_par = new TGraphErrors**[n_m];
  for(int i_m = 0; i_m < n_m; i_m++) {
    g_par[i_m] = new TGraphErrors*[n_p];
    TFile *fin = new TFile(Form("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/%s/files/mfit_2.root", modn[i_m].c_str()));
    int n_pv = n_p;
    for(int i_p = 0; i_p < n_pv; i_p++) {
      fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_m][i_p]);
    }
    fin->Close();
  }
  
  // scale all graphs for plotting
  TGraphErrors ***g_par_s = new TGraphErrors**[n_m];
  double mult[] = {100., 1., 1e3, 1e3, 1e3, 1., 1., 1., 1, 100, 1e3, 100.};
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

  for(int i_p = 0; i_p < 12; i_p++) {

    if(i_p == 4 || i_p == 9 || i_p == 10) continue;

    if(i_p == 1 || i_p == 7) c->SetLogy();
    else c->SetLogy(0);
    
    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.8);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("Prompt %s vs p_{T}", partit[i_p].c_str()));

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);

    // drawing base+G model
    // free pars - NS, NB, lambda, f_Bg
    if(i_p == 1 || i_p > 6) {
      for(int i_m = 0; i_m < n_m; i_m++) {
	g_par_s[i_m][i_p]->SetMarkerStyle(mc[i_m]);
	g_par_s[i_m][i_p]->SetMarkerSize(.75);
	g_par_s[i_m][i_p]->SetLineColor(pc[i_m]);
	g_par_s[i_m][i_p]->SetMarkerColor(pc[i_m]);
	g_par_s[i_m][i_p]->Draw("p");
      }
    }
    // linear or constant pars 
    else {
      for(int i_m = 0; i_m < n_m; i_m++) {
	g_par_s[i_m][i_p]->SetMarkerStyle(mc[i_m]);
	g_par_s[i_m][i_p]->SetMarkerSize(.75);
	g_par_s[i_m][i_p]->SetLineColor(pc[i_m]);
	g_par_s[i_m][i_p]->SetMarkerColor(pc[i_m]);
	g_par_s[i_m][i_p]->SetFillColorAlpha(pc[i_m], 0.5);
	g_par_s[i_m][i_p]->Draw("pce3");
      }
    }

    // SPECIAL CASES - COMBINED PLOTS 
    // if we're plotting par f_CB, add f_G
    if( i_p == 0) {
      for(int i_m = 0; i_m < n_m; i_m++) {
	g_par_s[i_m][9]->SetMarkerStyle(22);
	g_par_s[i_m][9]->SetMarkerColor(pc[i_m]);
	g_par_s[i_m][9]->SetLineColor(pc[i_m]);
	g_par_s[i_m][9]->SetFillColorAlpha(pc[i_m], 0.5);
	g_par_s[i_m][9]->Draw("pce3");
      }
      
      TLatex lc;
      lc.SetTextSize(0.03);
      lc.DrawLatex(30, 60, "f_{CB_{1}}");
      lc.DrawLatex(30, 10, "f_{G}");
    }
    
    // if we're plotting par sig1, add sig2 and sigG
    else if( i_p == 3) {
      for(int i_m = 0; i_m < n_m; i_m++) {
	g_par_s[i_m][4]->SetMarkerStyle(22);
	g_par_s[i_m][4]->SetMarkerColor(pc[i_m]);
	g_par_s[i_m][4]->SetLineColor(pc[i_m]);
	g_par_s[i_m][4]->SetFillColorAlpha(pc[i_m], 0.5);
	g_par_s[i_m][4]->Draw("pce3");

      	g_par_s[i_m][10]->SetMarkerStyle(25);
	g_par_s[i_m][10]->SetMarkerColor(pc[i_m]);
	g_par_s[i_m][10]->SetLineColor(pc[i_m]);
	g_par_s[i_m][10]->SetFillColorAlpha(pc[i_m], 0.5);
	g_par_s[i_m][10]->Draw("pce3");
      }
      
      TLatex lc;
      lc.SetTextSize(0.03);
      lc.DrawLatex(30, 23, "#sigma_{1}");
      lc.DrawLatex(30, 35, "#sigma_{2}");
      lc.DrawLatex(30, 65, "#sigma_{G}");
    }

    
    int isLog = 0;
    if(i_p == 1 || i_p == 7) isLog = 1;
    
    TLatex lc;
    lc.SetTextSize(0.03);
    for(int i_m = 0; i_m < n_m; i_m++) {
      lc.SetTextColor(pc[i_m]);
      double xp = 100;
      lc.DrawLatex(xp, getPos(parmin[i_p], parmax[i_p], 0.9-0.07*i_m, isLog), legn[i_m].c_str());
    }

    c->SaveAs(Form("plots/mass/par_%s.pdf", parsave[i_p].c_str()));
    c->Clear();
  }

  c->Destructor();

}
