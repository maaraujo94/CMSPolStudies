// macro to compare fit parameters w/ and w/o Gaussian

#import "../ptbins.C"

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

double doLin(double x, double *par) {
  return par[0]*x + par[1];
}

//aux func calculates propagated uncertainty
double parErr(const int npar, double (*fpar)(double, double *), double x, double *par, double *epar, double cov[])
{
  double ln = 1e4;
  double fval = fpar(x, par);
  double ferr = 0;

  double devp[npar];
  double par_var[npar];
  // get the aux array of pars without affecting main one
  for(int i_p = 0; i_p < npar; i_p++)
    par_var[i_p] = par[i_p];
  // get the array of deviations
  for(int i_p = 0; i_p < npar; i_p++) {
    par_var[i_p] += epar[i_p]/ln;
    devp[i_p] = (fpar(x, par_var)-fval)/(epar[i_p]/ln);
    par_var[i_p] -= epar[i_p]/ln;
  }

  // get the uncertainty as the sum of contributions
  for(int i_p = 0; i_p < npar; i_p++) {
    for(int j_p = 0; j_p < npar; j_p++) {
      ferr += devp[i_p]*devp[j_p]*cov[i_p*npar+j_p];
    }
  }
  ferr = sqrt(ferr);
  return ferr;
}

void plotDMPars()
{
  // aux arrays
  int pc[] = {kBlack, kRed};
  const int n_p = 9, n_m = 2;
  string modn[] = {"_0", "_1"};
  string legn[] = {"#alpha free", "#alpha constant"};

  string parlab[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda"};
  string parsave[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "tbkg"};
  string partit[] = {"f", "N_{SR}", "#mu", "#sigma", "#sigma_{2}", "n", "#alpha", "N_{BG}", "t_{Bg}"};
  string parax[] = {"f (%)", "N_{SR} per 1 GeV", "#mu (MeV)", "#sigma (MeV)", "#sigma_{2} (MeV)", "n", "#alpha", "N_{BG} per 1 GeV", "t_{Bg} (GeV)"};
  
  double parmin[] = {0,    2e0, 3620, 0,   32, 2.0, 1.3, 6e1, 0};
  double parmax[] = {100., 2e3, 3720, 60, 46, 3.0, 5., 3e4, 4};
 
  // initialize tgraphs for parameters
  TGraphErrors ***g_par = new TGraphErrors**[n_m];
  TFitResult **fitres = new TFitResult*[n_m];
  for(int i_m = 0; i_m < n_m; i_m++) {
    g_par[i_m] = new TGraphErrors*[n_p];
    TFile *fin = new TFile(Form("files/mfit%s.root", modn[i_m].c_str()));
    fitres[i_m] = (TFitResult*)fin->Get("fitres");
    for(int i_p = 0; i_p < n_p; i_p++) {
      fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_m][i_p]);
    }
    fin->Close();
  }
  
  // scale all graphs for plotting
  TGraphErrors ***g_par_s = new TGraphErrors**[n_m];
  double pt_min, pt_max;
  double sigL[2][nPtBins], esigL[2][nPtBins], cov_d[2][nPtBins];
  // f_CB1 is percentage
  // mu_m, sigma_CB1, sigma_CB2 GeV->MeV
  // N_SR and N_BG scaled by bin width
  // n, alpha, t_Bg (GeV) not scaled from original
  double mult[] = {100., 1., 1e3, 1e3, 1e3, 1., 1., 1., 1.};
  
  for(int i_m = 0; i_m < n_m; i_m++) {
    double *xv = g_par[i_m][0]->GetX();
    double *xe = g_par[i_m][0]->GetEX();
    int n = g_par[i_m][0]->GetN();
    pt_min = xv[0]-xe[0]-5;
    pt_max = xv[n-1]+xe[n-1]+5;
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
      if(i != 3 && i != 4)      
	g_par_s[i_m][i] = new TGraphErrors(n, xv, yv, xe, ye);
    }
    // also get sigmas from the linear fits
    
    // get diagonal of cov matrix
    for(int i = 0; i < nPtBins; i++) {
      cov_d[0][i] = fitres[i_m]->GetCovarianceMatrix()[3*nPtBins][3*nPtBins+1];
      cov_d[1][i] = fitres[i_m]->GetCovarianceMatrix()[3*nPtBins][4*nPtBins+1];
    }
    
    for(int i = 0; i < nPtBins; i++){
      double pt_val = 0.5*(ptBins[i+1]+ptBins[i]);
      for(int j = 0; j < 2; j++) {
	// define input for the dynamic function
	double par_vec[] = {fitres[i_m]->Parameter(3*nPtBins), fitres[i_m]->Parameter((3+j)*nPtBins+1)};
	double epar_vec[] = {fitres[i_m]->ParError(3*nPtBins), fitres[i_m]->ParError((3+j)*nPtBins+1)};
	double cov_mat[2][2] = {{epar_vec[0]*epar_vec[0], cov_d[j][i]}, {cov_d[j][i], epar_vec[1]*epar_vec[1]}};
	double *cov_ptr = cov_mat[0];
	
	sigL[j][i] = doLin(pt_val, par_vec)*mult[j+3];
	esigL[j][i] = parErr(2, doLin, pt_val, par_vec, epar_vec, cov_ptr)*mult[j+3];
      }
    }
    g_par_s[i_m][3] = new TGraphErrors(n, xv, sigL[0], xe, esigL[0]);
    g_par_s[i_m][4] = new TGraphErrors(n, xv, sigL[1], xe, esigL[1]);
  }

  // get extra graph for alpha to make it prettier
  int n = g_par_s[1][6]->GetN();
  double *xv = g_par_s[1][6]->GetX();
  double *xe = g_par_s[1][6]->GetEX();
  double *yv = g_par_s[1][6]->GetY();
  double *ye = g_par_s[1][6]->GetEY();

  double xvf[n+2], xef[n+2], yvf[n+2], yef[n+2];
  for(int i = 0; i < n+2; i++) {
    if (i == 0) {
      xvf[i] = xv[0]-xe[0];
      xef[i] = 0;
      yvf[i] = yv[0];
      yef[i] = ye[0];
    }
    else if (i == n+1) {
      xvf[i] = xv[n-1]+xe[n-1];
      xef[i] = 0;
      yvf[i] = yv[n-1];
      yef[i] = ye[n-1];
    }
    else {
      xvf[i] = xv[i-1];
      xef[i] = xe[i-1];
      yvf[i] = yv[i-1];
      yef[i] = ye[i-1];
    }
  }
  TGraphErrors *g_alpha = new TGraphErrors(n+2, xvf, yvf, xef, yef);
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.03);
  
  for(int i_p = 0; i_p < n_p; i_p++) {

    if(i_p == 4) continue;

    if(i_p == 1 || i_p == 7) c->SetLogy();
    else c->SetLogy(0);
    
    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.8);
    if(i_p == 8)    fl->GetYaxis()->SetTitleOffset(1.6);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("Prompt %s vs p_{T}", partit[i_p].c_str()));

    TLegend *leg = new TLegend(0.775, 0.775, 1.075, 0.875);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(kWhite,0);

    // drawing free vs const alpha
    // free par NS
    if(i_p == 1) {
      for(int i_n = 0; i_n < n_m; i_n++) {
	g_par_s[i_n][i_p]->SetMarkerStyle(20);
	g_par_s[i_n][i_p]->SetMarkerSize(.75);
	g_par_s[i_n][i_p]->SetLineColor(pc[i_n]);
	g_par_s[i_n][i_p]->SetMarkerColor(pc[i_n]);
	g_par_s[i_n][i_p]->Draw("p");
      }
    }

    // free pars NB, lambda (just last fit)
    else if(i_p == 7 || i_p == 8) {
      for(int i_n = n_m-1; i_n < n_m; i_n++) {
	g_par_s[i_n][i_p]->SetMarkerStyle(20);
	g_par_s[i_n][i_p]->SetMarkerSize(.75);
	g_par_s[i_n][i_p]->SetLineColor(pc[0]);
	g_par_s[i_n][i_p]->SetMarkerColor(pc[0]);
	g_par_s[i_n][i_p]->Draw("p");
      }

    }

    // mixed free-constant - alpha
    else if(i_p == 6) {
      g_par_s[0][i_p]->SetMarkerStyle(20);
      g_par_s[0][i_p]->SetMarkerSize(.75);
      g_par_s[0][i_p]->SetLineColor(pc[0]);
      g_par_s[0][i_p]->SetMarkerColor(pc[0]);
      g_par_s[0][i_p]->Draw("p");
      
      g_alpha->SetMarkerStyle(20);
      g_alpha->SetMarkerSize(.75);
      g_alpha->SetLineColor(kRed);
      g_alpha->SetMarkerColor(kRed);
      g_alpha->SetFillColorAlpha(kRed, 0.5);
      g_alpha->Draw("ce3");
    }
    
    // linear or constant pars - f, mu, n
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
    // if we're plotting par sig1, add sig2
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
      }
      
      leg->AddEntry(g_par_s[0][3], "#sigma_{1}", "pl");
      leg->AddEntry(g_par_s[0][4], "#sigma_{2}", "pl");
      leg->Draw();
    }

    int isLog = 0;
    if(i_p == 1 || i_p == 7 ) isLog = 1;

    TLatex lc;
    lc.SetTextSize(0.03);
    double xp = getPos(pt_min, pt_max, 0.1, 0);
    if(i_p == 7 || i_p == 1) xp = getPos(pt_min, pt_max, 0.7, 0);
    for(int i_m = 0; i_m < n_m; i_m++) {
      lc.SetTextColor(pc[i_m]);
      double yp = getPos(parmin[i_p], parmax[i_p], 0.9-0.07*i_m, isLog);
      if(i_p != 7 && i_p != 8) 
	lc.DrawLatex(xp, yp, legn[i_m].c_str());
    }

    c->SaveAs(Form("plots/mass/par_%s.pdf", parsave[i_p].c_str()));
    c->Clear();
  }
  
  c->Destructor();

}
