// macro to propagate the uncertainties to f_NP analytically
// then save it as a 2d histo unc band

#import "../ptbins.C"

// global definitions
TF1 *fNP_SR;
double pr_lim = 0.05;

// define negative exponential only for positive x
double pos_exp(double x, double ld)
{
  if(x > 0) return exp(-x/ld);
  else return 0;
}

double intfNP(double x, double *par)
{
  double N_NP = par[0];
  double f1 = par[1];
  double f2 = par[2];
  double mu = par[3];
  double sig1 = par[4];
  double sigR21 = par[5]; 
  double sigR31 = par[6]; 
  double tnp_SR = par[7];

  double sig2 = sigR21*sig1;
  double sig3 = sigR31*sig1;

  fNP_SR->SetParameters(N_NP, f1, f2, mu, sig1, sig2, sig3, tnp_SR);
  return fNP_SR->Integral(-pr_lim, pr_lim);
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

void fnpProp()
{
  // PART 1: get the uncertainty of the f_NP

  // prepare lt histograms
  TH2D *h_d2d = new TH2D();
  TH1D **h_d1d = new TH1D*[nPtBins];
  TFile *fin = new TFile("files/ltStore.root");
  fin->GetObject("ltH_SR", h_d2d);
  h_d2d->SetDirectory(0);
  fin->Close();
  for(int ip = 0; ip < nPtBins; ip++) {
    h_d1d[ip] = h_d2d->ProjectionX(Form("ltH%.0f", ptBins[ip]), ip+1, ip+1);
  }

  // get all needed input
  int tbins = h_d1d[0]->GetNbinsX();
  double lowt = h_d1d[0]->GetBinLowEdge(1);
  double hit = h_d1d[0]->GetBinLowEdge(tbins+1);

  // define the resolution (=PR) function
  TF1 *fres = new TF1("fres", "[0]*([1]*TMath::Gaus(x, [3],[4]) + [2]*TMath::Gaus(x, [3], [5])+ (1.-[1]-[2])*TMath::Gaus(x, [3], [6]))", 5*lowt, 5*hit);

  // define the NP function by convolution
  TF1 *fexp = new TF1("fexp", "pos_exp(x,[0])", 5*lowt, 5*hit);
  TF1Convolution *fcNP = new TF1Convolution(fres, fexp, 5*lowt, 5*hit);
  fcNP->SetRange(5*lowt, 5*hit);
  fcNP->SetNofPointsFFT(1000);
  fNP_SR = new TF1("fNP_SR", *fcNP, lowt, hit, fcNP->GetNpar());

  // get fit parameters from fitres
  TFile *inNP = new TFile("files/ltfitres2d.root");
  int n_par = 13; 
  // get fitres to get cov in each pT bin
  TFitResult *fitres = new TFitResult();
  inNP->GetObject("fitres", fitres);
  inNP->Close();

  // start cycle of calculations
  TH1D *h_fnp_psi = new TH1D("h_fnp_psi", "Run 2 f_{NP}", nPtBins, ptBins);
  double fnp[nPtBins], efnp[nPtBins], evt_all[nPtBins];
  int di[] = {1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1}; // determines which pars are constant  
  double min_bin = h_d1d[0]->GetXaxis()->FindBin(-pr_lim+1e-6);
  double max_bin = h_d1d[0]->GetXaxis()->FindBin(pr_lim-1e-6);
  
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    // N_histo has no uncertainty: just histogram integral
    evt_all[i_pt] = h_d1d[i_pt]->Integral(min_bin, max_bin, "width");

    // get vectors of fit parameters
    double par_vec[8], epar_vec[8], cov_mat[8][8];
    for(int i_p = 0; i_p < 8; i_p++) {
      par_vec[i_p] = fitres->Parameter((i_p+1)*nPtBins + di[i_p+1] * i_pt);
      epar_vec[i_p] = fitres->ParError((i_p+1)*nPtBins + di[i_p+1] * i_pt);
      for(int j_p = 0; j_p < 8; j_p++) {
	cov_mat[i_p][j_p] = fitres->GetCovarianceMatrix()[(i_p+1)*nPtBins + di[i_p+1] * i_pt][(j_p+1)*nPtBins + di[j_p+1] * i_pt];
      }
    }
    double *cov_ptr = cov_mat[0];

    fnp[i_pt] = intfNP(0, par_vec)/evt_all[i_pt] * 100.;
    efnp[i_pt] = parErr(8, intfNP, 0, par_vec, epar_vec, cov_ptr)/evt_all[i_pt] * 100.;
    
    h_fnp_psi->SetBinContent(i_pt+1, fnp[i_pt]);
    h_fnp_psi->SetBinError(i_pt+1, efnp[i_pt]);
  }

  // now get the fNP_bkg (no prop uncert bc it's a 2-step fit)
  TF1 *fexp_bkg = new TF1("fexp_bkg", "[0]*pos_exp(x,[1]) + (1.-[0])*pos_exp(x,[2])", 5*lowt, 5*hit);
  TF1Convolution *fcNP_bkg = new TF1Convolution(fres, fexp_bkg, 5*lowt, 5*hit);
  fcNP_bkg->SetRange(5*lowt, 5*hit);
  fcNP_bkg->SetNofPointsFFT(1000);
  TF1 *fNP_bkg = new TF1("fNP_bkg", *fcNP_bkg, lowt, hit, fcNP_bkg->GetNpar());
  // define the bkg exp for N_bkg adjustment
  TF1 *fexp_nbkg = new TF1("fexp_nbkg", "[0]*[1]*pos_exp(x,[2]) + [0]*(1.-[1])*pos_exp(x,[3])", 5*lowt, 5*hit);

  // get corrected N_bkg, then get fNP_bkg
  double b_pars[10];
  TH1D *h_fnp_bkg = new TH1D("h_fnp_bkg", "Run 2 f_{NP_{bkg}}", nPtBins, ptBins);
  double fnp_bkg[nPtBins], efnp_bkg[nPtBins], xv[nPtBins], xe[nPtBins];

  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    // setting parameters
    b_pars[0] = 1; // N_bkg (later)
    for(int i_p = 1; i_p < 5; i_p++) // up to sig1
      b_pars[i_p] = fitres->Parameter((i_p+1)*nPtBins + di[i_p+1] * i_pt);
    for(int i_p = 5; i_p < 7; i_p++) // sig2,3
      b_pars[i_p] = fitres->Parameter((i_p+1)*nPtBins + di[i_p+1] * i_pt) * fitres->Parameter(5*nPtBins + di[5] * i_pt);
    for(int i_p = 7; i_p < 10; i_p++) // bkg pars
      b_pars[i_p] = fitres->Parameter((i_p+3)*nPtBins + di[i_p+3] * i_pt);

    // scale N_bkg par
    for(int i_p = 0; i_p < 10; i_p++) 
      fNP_bkg->SetParameter(i_p, b_pars[i_p]);
    fexp_nbkg->SetParameters(fitres->Parameter(9*nPtBins + i_pt), b_pars[7], b_pars[8], b_pars[9]);
    double N_bkg_conv = fexp_nbkg->Eval(0.4)/fNP_bkg->Eval(0.4);
    fNP_bkg->SetParameter(0, N_bkg_conv);

    // integrate, get fraction
    fnp_bkg[i_pt] = fNP_bkg->Integral(-pr_lim, pr_lim)/evt_all[i_pt] * 100;
    efnp_bkg[i_pt] = 0;
    
    h_fnp_bkg->SetBinContent(i_pt+1, fnp_bkg[i_pt]);
    h_fnp_bkg->SetBinError(i_pt+1, efnp_bkg[i_pt]);

    // also getting pt values for tgraph
    xv[i_pt] = 0.5*(ptBins[i_pt+1]+ptBins[i_pt]);
    xe[i_pt] = 0.5*(ptBins[i_pt+1]-ptBins[i_pt]);
  }
  
  // plotting in pT
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.11);
  c->SetTopMargin(0.015);
  
  TH1F *fr1 = c->DrawFrame(ptBins[0]-5, 0, ptBins[nPtBins]+5, 25);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f_{NP#psi} (%)");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("");

  h_fnp_psi->SetStats(0);
  h_fnp_psi->SetMarkerStyle(20);
  h_fnp_psi->SetMarkerSize(.5);
  h_fnp_psi->SetMarkerColor(kBlack);
  h_fnp_psi->SetLineColor(kBlack);
  h_fnp_psi->Draw("e1 same");

  TGraphErrors *g_fnp_bkg = new TGraphErrors(nPtBins, xv, fnp_bkg, xe, efnp_bkg);
  //g_fnp_bkg->SetStats(0);
  g_fnp_bkg->SetMarkerStyle(20);
  g_fnp_bkg->SetMarkerSize(.5);
  g_fnp_bkg->SetMarkerColor(kGreen+2);
  g_fnp_bkg->SetLineColor(kGreen+2);
  //g_fnp_bkg->Draw("psame");

  TH1D *h_fnp_sum = new TH1D("h_fnp_sum", "", nPtBins, ptBins);
  h_fnp_sum->Sumw2();
  h_fnp_sum->Add(h_fnp_psi, h_fnp_bkg, 1, 1);
  h_fnp_sum->SetStats(0);
  h_fnp_sum->SetMarkerStyle(20);
  h_fnp_sum->SetMarkerSize(.5);
  h_fnp_sum->SetMarkerColor(kRed);
  h_fnp_sum->SetLineColor(kRed);
  //h_fnp_sum->Draw("e1 same");

  TLegend* leg = new TLegend(0.77, 0.825, 1.07, 0.975);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  leg->AddEntry(h_fnp_psi, "f_{NP#psi}", "pl");
  leg->AddEntry(g_fnp_bkg, "f_{NPbkg}", "pl");
  leg->AddEntry(h_fnp_sum, "f_{NP}", "pl");
  //leg->Draw();
  
  c->SaveAs("plots/fNP_unc.pdf");
  c->Clear();
  
  // PART 2: generate f_np 2d histo

  // get costh binning from the stored data histos
  TFile *infile = new TFile("files/histoStore.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("PRH"), hist);

  // get the binning
  int nBinsX = hist->GetNbinsX();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/nBinsX;

  infile->Close();

  // f_np(pT) but generating 2d map so it's easier to apply uncertainties
  TH2D *h_fnp2d = new TH2D("h_fnp2d", "Run 2 f_{NP}", nBinsX, minX, maxX, nPtBins, ptBins);
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    // same result for all costh bins
    for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
      h_fnp2d->SetBinContent(i_cos+1, i_pt+1, h_fnp_psi->GetBinContent(i_pt+1));
      h_fnp2d->SetBinError(i_cos+1, i_pt+1, h_fnp_psi->GetBinError(i_pt+1));
    }
  }
  
  c->Destructor();

  // scale fractions down from percentage
  h_fnp_psi->Scale(1./100.);
  h_fnp_bkg->Scale(1./100.);
  h_fnp_sum->Scale(1./100.);
  h_fnp2d->Scale(1./100.);
  
  TFile *fout = new TFile("files/NPFrac.root", "recreate");
  h_fnp_psi->SetName("fnppsi_unc");
  h_fnp_psi->Write();
  h_fnp_bkg->SetName("fnpbkg_unc");
  h_fnp_bkg->Write();
  h_fnp_sum->SetName("fnp_unc");
  h_fnp_sum->Write();
  h_fnp2d->SetName("h_fnp");
  h_fnp2d->Write();
  fout->Close();

}
