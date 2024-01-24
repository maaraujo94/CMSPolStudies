// macro to propagate the uncertainties to f_bkg analytically

#import "../ptbins.C"

TF1 *fbkg;
// mass region ranges
double m_min[] = {2.94, 3.0, 3.21};
double m_max[] = {2.95, 3.2, 3.26};

// background function
double bkg_exp(double m, double p1, double p2)
{
  return p1 * exp( - m / p2 );
}

double intfbkg(double x, double *par)
{
  double N_Bg = par[0];
  double lambda = par[1];

  fbkg->SetParameters(N_Bg, lambda);
  return fbkg->Integral(m_min[1], m_max[1]);
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


void fbkgProp()
{
  // PART 1: get the uncertainty of the f_bkg

  // prepare mass histograms
  TH1D **h_d1d = new TH1D*[nPtBins];
  TH2D *h_d2d = new TH2D();
  TFile *fin = new TFile("files/mStore.root");
  fin->GetObject("mH", h_d2d);
  h_d2d->SetDirectory(0);
  fin->Close();
  for(int ip = 0; ip < nPtBins; ip++) {
    h_d1d[ip] = h_d2d->ProjectionX(Form("mH%.0f", ptBins[ip]), ip+1, ip+1);
  }

  // define the bkg function
  fbkg = new TF1("fbkg", "bkg_exp(x,[0],[1])", m_min[0], m_max[2]);
  
  // get fit parameters - need to know which params are being used
  TFile *inBG = new TFile("files/mfit_2.root");
  int n_par = 2;
  // get fitres to get cov in each pT bin
  TFitResult* fitres = (TFitResult*)inBG->Get("fitres");
  inBG->Close();
    
  // start cycle of calculations
  double f_Bg[nPtBins], ef_Bg[nPtBins];
  double evt_all[nPtBins], int_Bg[nPtBins], eint_Bg[nPtBins];

  TH1D *h_fbkg = new TH1D("h_fbkg", "Run 2 f_{Bg}", nPtBins, ptBins);
  TH1D *h_intBg = new TH1D("h_intBg", "integ Bg in SR", nPtBins, ptBins);
  TH1D *h_nevt = new TH1D("h_nevt", "Evt count in SR", nPtBins, ptBins);

  double min_bin = h_d1d[0]->GetXaxis()->FindBin(m_min[1]+1e-6);
  double max_bin = h_d1d[0]->GetXaxis()->FindBin(m_max[1]-1e-6);

  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    // N_histo has no uncertainty: just histogram integral
    evt_all[i_pt] = h_d1d[i_pt]->Integral(min_bin, max_bin, "width");

    // get vectors of fit parameters
    double par_vec[n_par], epar_vec[n_par], cov_mat[n_par][n_par];
    for(int i_p = 0; i_p < n_par; i_p++) {
      par_vec[i_p] = fitres->Parameter((i_p+7)*nPtBins + i_pt);
      epar_vec[i_p] = fitres->ParError((i_p+7)*nPtBins + i_pt);
      for(int j_p = 0; j_p < n_par; j_p++) {
	cov_mat[i_p][j_p] = fitres->GetCovarianceMatrix()[(i_p+7)*nPtBins + i_pt][(j_p+7)*nPtBins + i_pt];
      }
    }
    double *cov_ptr = cov_mat[0];

    int_Bg[i_pt] = intfbkg(0, par_vec);
    f_Bg[i_pt] = int_Bg[i_pt]/evt_all[i_pt] * 100.;
    eint_Bg[i_pt] = parErr(n_par, intfbkg, 0, par_vec, epar_vec, cov_ptr);
    ef_Bg[i_pt] = eint_Bg[i_pt]/evt_all[i_pt] * 100.;
    
    h_fbkg->SetBinContent(i_pt+1, f_Bg[i_pt]);
    h_fbkg->SetBinError(i_pt+1, ef_Bg[i_pt]);
    h_intBg->SetBinContent(i_pt+1, int_Bg[i_pt]);
    h_intBg->SetBinError(i_pt+1, eint_Bg[i_pt]);
    h_nevt->SetBinContent(i_pt+1, evt_all[i_pt]);
    h_nevt->SetBinError(i_pt+1, 0);
   }

  // plotting in pT
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);

  TH1F *fr1 = c->DrawFrame(ptBins[0]-5, 0, ptBins[nPtBins]+5, 15);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f_{Bg} (%)");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("f_{Bg} vs p_{T}");

  h_fbkg->SetStats(0);
  h_fbkg->SetMarkerStyle(20);
  h_fbkg->SetMarkerSize(.5);
  h_fbkg->SetMarkerColor(kBlack);
  h_fbkg->SetLineColor(kBlack);
  h_fbkg->Draw("e1 same");
  
  c->SaveAs("plots/fBG_unc.pdf");

  // PART 2: generate f_bkg histo

  // get costh binning from the stored data histos
  TFile *infile = new TFile("../PR_fit/files/histoStore.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("PRH"), hist);

  // get the binning
  int nBinsX = hist->GetNbinsX();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/nBinsX;

  infile->Close();

  // f_bkg(pT) but generating 2d map so it's easier to apply uncertainties
  TH2D *h_fbkg2d = new TH2D("h_fbkg2d", "Run 2 f_{Bg}", nBinsX, minX, maxX, nPtBins, ptBins);
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    // same result for all costh bins
    for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
      h_fbkg2d->SetBinContent(i_cos+1, i_pt+1, h_fbkg->GetBinContent(i_pt+1));
      h_fbkg2d->SetBinError(i_cos+1, i_pt+1, h_fbkg->GetBinError(i_pt+1));
    }

  }
  c->Destructor();

  // scale fractions down from percentage
  h_fbkg->Scale(1./100.);
  h_fbkg2d->Scale(1./100.);

  TFile *fout = new TFile("files/bkgFrac.root", "recreate");
  h_fbkg->SetName("fbkg_unc");
  h_fbkg->Write();
  h_fbkg2d->SetName("h_fbkg");
  h_fbkg2d->Write();
  h_intBg->Write();
  h_nevt->Write();
  fout->Close();

}
