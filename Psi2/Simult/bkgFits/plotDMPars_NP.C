// macro to compare fit parameters w/ diff ctau lims
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


void plotDMPars_NP()
{
  // aux arrays
  int pc[] = {kBlack, kBlue, kViolet};
  const int n_p = 9, n_m = 1;
  string modn[] = {""};
  string legn[] = {""};

  string parlab[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "lambda"};
  string parsave[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "NB", "tbkg"};
  string partit[] = {"f", "N_{SR}", "#mu", "#sigma", "#sigma_{2}", "n", "#alpha", "N_{BG}", "t_{Bg}"};
  string parax[] = {"f (%)", "N_{SR} per 1 GeV", "#mu (MeV)", "#sigma (MeV)", "#sigma_{2} (MeV)", "n", "#alpha", "N_{BG} per 1 GeV", "t_{Bg} (MeV)"};
  
  double parmin[] = {0,    1e0, 3620, 0,   32, 2.0, 1.0, 2e5, 0};
  double parmax[] = {100., 2e3, 3720, 80, 46, 3.0, 2.3, 3e8, 600};
 
  // initialize tgraphs for parameters
  TGraphErrors ***g_par = new TGraphErrors**[n_m];
  TFitResult **fitres = new TFitResult*[n_m];
  for(int i_m = 0; i_m < n_m; i_m++) {
    g_par[i_m] = new TGraphErrors*[n_p];
    TFile *fin = new TFile(Form("%sfiles/mfit_NP.root", modn[i_m].c_str()));
    fitres[i_m] = (TFitResult*)fin->Get("fitres");
    for(int i_p = 0; i_p < n_p; i_p++) {
      fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_m][i_p]);
    }
    fin->Close();
  }
  
  // scale all graphs for plotting
  TGraphErrors ***g_par_s = new TGraphErrors**[n_m];
  double mult[] = {100., 1., 1e3, 1e3, 1e3, 1., 1., 1., 1e3};
  double pt_min, pt_max;
  double sigL[2][nPtBins], esigL[2][nPtBins], cov_d[2][nPtBins];
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
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.12);
  
  for(int i_p = 0; i_p < n_p; i_p++) {

    if(i_p == 4) continue;

    if(i_p == 1 || i_p == 7) c->SetLogy();
    else c->SetLogy(0);
    
    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.8);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("Non-prompt %s vs p_{T}", partit[i_p].c_str()));

    TLegend *leg = new TLegend(0.775, 0.775, 1.075, 0.875);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(kWhite,0);
    
    // drawing base + extended ctau model
    // free pars - NS, NB, lambda
    if(i_p == 1 || i_p > 6) {
      for(int i_n = 0; i_n < n_m; i_n++) {
	g_par_s[i_n][i_p]->SetMarkerStyle(20);
	g_par_s[i_n][i_p]->SetMarkerSize(.75);
	g_par_s[i_n][i_p]->SetLineColor(pc[i_n]);
	g_par_s[i_n][i_p]->SetMarkerColor(pc[i_n]);
	g_par_s[i_n][i_p]->Draw("p");
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
      //lc.DrawLatex(xp, yp, legn[i_m].c_str());
     }
     
    c->SaveAs(Form("plots/massNP/par_%s.pdf", parsave[i_p].c_str()));
    c->Clear();
  }
  
  c->Destructor();

}
