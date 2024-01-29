#import "../../ptbins.C"

double doProd(double x, double *par) {
  return par[0]*par[1];
}

double dof3(double x, double *par) {
  return 100.-par[0]-par[1];
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

void plotLtPars2d()
{
  const int n_p = 9;
  double cov_df[nPtBins];

  string parlab[] = {"N_PR", "N_NP", "f", "f_2", "mu", "sig", "sig2", "sig3", "t_NP"};
  string par_read[] = {"N_PR", "N_NP", "f1", "f2", "mu", "sig1", "sigR21", "sigR31", "t_NPSR"};
  string partit[] = {"N_{PR}", "N_{NP}", "f_{G}", "f_{G_{2}}", "#mu_{c#tau}", "#sigma", "#sigma_{G_{2}}", "#sigma_{G_{3}}", "t_{NP_{#psi}}"};
  string par_unit[] = {" per 1 GeV", " per 1 GeV", " (%)", " (%)", " (#mum)", " (#mum)", " (#mum)", " (#mum)", " (#mum)"};

  double parmin[] = {1e2, 2e2, 0,   0,  -1,  0,  0,  0, 0,  };
  double parmax[] = {6e4, 2e5, 100, 100, 1,  140, 40, 40, 600.};
  double mult[] = {1, 1, 100., 100, 1e3, 1e3, 1, 1, 1e3};

  TGraphErrors **g_parR = new TGraphErrors*[n_p];
  TFile *fin = new TFile("files/ltfitres2d.root");
  TFitResult *fitres = (TFitResult*)fin->Get("fitres");
  for(int j = 0; j < n_p; j++) {
    g_parR[j] = (TGraphErrors*)fin->Get(Form("fit_%s", par_read[j].c_str()));
  }
  fin->Close();
  for(int i = 0; i < nPtBins; i++) // f1, f2 constant
    cov_df[i] = fitres->GetCovarianceMatrix()[2*nPtBins][3*nPtBins]*mult[2]*mult[3];

  double *xv = g_parR[0]->GetX();
  double *xe = g_parR[0]->GetEX();
  int n = g_parR[0]->GetN();
  // scale all graphs for plotting
  TGraphErrors **g_par_s = new TGraphErrors*[n_p];

  for(int i = 0; i < n_p; i++) {
    if(i != 6 && i != 7) { // not doing sig2, sig3
      double *yv = g_parR[i]->GetY();
      double *ye = g_parR[i]->GetEY();

      for(int j = 0; j < n; j++) {
	if(i < 2) { // scaling N_PR and N_NP
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
  }

  // also get sigmas from the fitted ratios
  double sigP[2][nPtBins], esigP[2][nPtBins], cov_d[2][nPtBins];
  // get diagonal of cov matrix
  for(int i = 0; i < nPtBins; i++) {
    cov_d[0][i] = fitres->GetCovarianceMatrix()[5*nPtBins+i][6*nPtBins]*mult[5]*mult[6];
    cov_d[1][i] = fitres->GetCovarianceMatrix()[5*nPtBins+i][7*nPtBins]*mult[5]*mult[7];
  }
  
  for(int i = 0; i < nPtBins; i++){
    for(int j = 0; j < 2; j++) {
      // define input for the dynamic function
      double par_vec[] = {g_parR[6+j]->GetY()[i], g_par_s[5]->GetY()[i]};
      double epar_vec[] = {g_parR[6+j]->GetEY()[i], g_par_s[5]->GetEY()[i]};
      double cov_mat[2][2] = {{epar_vec[0]*epar_vec[0], cov_d[j][i]}, {cov_d[j][i], epar_vec[1]*epar_vec[1]}};
      double *cov_ptr = cov_mat[0];

      sigP[j][i] = doProd(0, par_vec);
      esigP[j][i] = parErr(2, doProd, 0, par_vec, epar_vec, cov_ptr);
    }
  }
  g_par_s[6] = new TGraphErrors(n, xv, sigP[0], xe, esigP[0]);
  g_par_s[7] = new TGraphErrors(n, xv, sigP[1], xe, esigP[1]);

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.11);

  // get also f3
  double v_f3[nPtBins], e_f3[nPtBins];
  for(int i = 0 ; i < nPtBins ; i++) {
    // define input for the dynamic function
    double par_vec[] = {g_par_s[2]->GetY()[i], g_par_s[3]->GetY()[i]};
    double epar_vec[] = {g_par_s[2]->GetEY()[i], g_par_s[3]->GetEY()[i]};
    double cov_mat[2][2] = {{epar_vec[0]*epar_vec[0], cov_df[i]}, {cov_df[i], epar_vec[1]*epar_vec[1]}};
    double *cov_ptr = cov_mat[0];

    v_f3[i] = dof3(0, par_vec);
    e_f3[i] = parErr(2, dof3, 0, par_vec, epar_vec, cov_ptr);
  }
  TGraphErrors *g_f3_s = new TGraphErrors(n, xv, v_f3, xe, e_f3);

  for(int j = 0; j < n_p; j++) {

    if(j < 2 ) c->SetLogy();
    else c->SetLogy(0);
    if(j == 3 || j == 6 || j == 7) continue;

    TH1F *fp = c->DrawFrame(ptBins[0]-5, parmin[j], ptBins[nPtBins]+5, parmax[j]);
    fp->SetXTitle("p_{T} (GeV)");
    fp->SetYTitle(Form("%s%s", partit[j].c_str(), par_unit[j].c_str()));
    fp->GetYaxis()->SetTitleOffset(1.5);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("%s vs p_{T}", partit[j].c_str()));

    g_par_s[j]->SetMarkerStyle(20);
    g_par_s[j]->SetMarkerSize(.75);
    g_par_s[j]->SetMarkerColor(kBlack);
    g_par_s[j]->SetLineColor(kBlack);
    g_par_s[j]->Draw("psame");

    if(j == 2) {
      int k = j+1;
      g_par_s[k]->SetMarkerStyle(25);
      g_par_s[k]->SetMarkerSize(.75);
      g_par_s[k]->SetMarkerColor(kBlue);
      g_par_s[k]->SetLineColor(kBlue);
      g_par_s[k]->Draw("psame");

      // draw f3 also
      g_f3_s->SetMarkerStyle(29);
      g_f3_s->SetMarkerColor(kViolet+2);
      g_f3_s->SetLineColor(kViolet+2);
      g_f3_s->Draw("psame");

      TLegend* leg = new TLegend(0.77, 0.725, 1.07, 0.875);
      leg->SetTextSize(0.03);
      leg->SetBorderSize(0);
      leg->SetFillColorAlpha(kWhite,0);
      leg->AddEntry(g_par_s[j], "f_{G_{1}}", "pl");
      leg->AddEntry(g_par_s[k], partit[k].c_str(), "pl");
      leg->AddEntry(g_f3_s, "f_{G_{3}}", "pl");
      leg->Draw();
    }

    if(j == 5) {
      int k = j+1;
      g_par_s[k]->SetMarkerStyle(25);
      g_par_s[k]->SetMarkerSize(.75);
      g_par_s[k]->SetMarkerColor(kBlue);
      g_par_s[k]->SetLineColor(kBlue);
      g_par_s[k]->Draw("psame");

      int l = j+2;
      g_par_s[l]->SetMarkerStyle(29);
      //g_par_s[l]->SetMarkerSize(.75);
      g_par_s[l]->SetMarkerColor(kViolet+2);
      g_par_s[l]->SetLineColor(kViolet+2);
      g_par_s[l]->Draw("psame");

      TLegend* leg = new TLegend(0.77, 0.725, 1.07, 0.875);
      leg->SetTextSize(0.03);
      leg->SetBorderSize(0);
      leg->SetFillColorAlpha(kWhite,0);
      leg->AddEntry(g_par_s[j], "#sigma_{G_{1}}", "pl");
      leg->AddEntry(g_par_s[k], partit[k].c_str(), "pl");
      leg->AddEntry(g_par_s[l], partit[l].c_str(), "pl");
      leg->Draw();
    }

    c->SaveAs(Form("plots/lifetime2d/par_%s.pdf", parlab[j].c_str()));
    c->Clear();
  }

  // draw zoomed in sigma plot with sigma_1,2
  int j = 5;
  c->SetLogy(0);
  
  TH1F *fp = c->DrawFrame(ptBins[0]-5, parmin[j], ptBins[nPtBins]+5, parmax[j+1]);
  fp->SetXTitle("p_{T} (GeV)");
  fp->SetYTitle(Form("%s%s", partit[j].c_str(), par_unit[j].c_str()));
  fp->GetYaxis()->SetTitleOffset(1.5);
  fp->GetYaxis()->SetLabelOffset(0.01);
  fp->SetTitle(Form("%s vs p_{T}", partit[j].c_str()));
  
  g_par_s[j]->SetMarkerStyle(20);
  g_par_s[j]->SetMarkerSize(.75);
  g_par_s[j]->SetMarkerColor(kBlack);
  g_par_s[j]->SetLineColor(kBlack);
  g_par_s[j]->Draw("psame");

  int k = j+1;
  g_par_s[k]->SetMarkerStyle(25);
  g_par_s[k]->SetMarkerSize(.75);
  g_par_s[k]->SetMarkerColor(kBlue);
  g_par_s[k]->SetLineColor(kBlue);
  g_par_s[k]->Draw("psame");

  
  TLegend* leg = new TLegend(0.77, 0.725, 1.07, 0.875);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  leg->AddEntry(g_par_s[j], "#sigma_{G_{1}}", "pl");
  leg->AddEntry(g_par_s[k], partit[k].c_str(), "pl");
  leg->Draw();
  
  c->SaveAs(Form("plots/lifetime2d/par_%s_zoom.pdf", parlab[j].c_str()));
  c->Clear();
  c->Destructor();

}
