// macro to propagate the uncertainties to f_NP analytically
// then save it as a 2d histo unc band

#import "../ptbins.C"

// global definitions
TF1 *fNP;
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
  double tnp = par[7];

  double sig2 = sigR21*sig1;
  double sig3 = sigR31*sig1;

  fNP->SetParameters(N_NP, f1, f2, mu, sig1, sig2, sig3, tnp);
  return fNP->Integral(-pr_lim, pr_lim);
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
  fin->GetObject("ltH", h_d2d);
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
  fNP = new TF1("fNP", *fcNP, lowt, hit, fcNP->GetNpar());
  
  // get fit parameters from fitres
  TFile *inNP = new TFile("files/ltfitres2d.root");
  // get fitres to get cov in each pT bin
  TFitResult *fitres = new TFitResult();
  inNP->GetObject("fitres", fitres);
  inNP->Close();

  // start cycle of calculations
  TH1D *h_fnp = new TH1D("h_fnp", "Run 2 f_{NP}", nPtBins, ptBins);
  double fnp[nPtBins], efnp[nPtBins], evt_all[nPtBins];
  int di[] = {1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1}; // determines which pars are constant  
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
    
    h_fnp->SetBinContent(i_pt+1, fnp[i_pt]);
    h_fnp->SetBinError(i_pt+1, efnp[i_pt]);
  }

  // plotting in pT
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.11);
  c->SetTopMargin(0.015);

  TH1F *fr1 = c->DrawFrame(ptBins[0]-5, 0, ptBins[nPtBins]+5, 50);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f_{NP} (%)");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("");

  h_fnp->SetStats(0);
  h_fnp->SetMarkerStyle(20);
  h_fnp->SetMarkerSize(.5);
  h_fnp->SetMarkerColor(kBlack);
  h_fnp->SetLineColor(kBlack);
  h_fnp->Draw("e1 same");
  
  c->SaveAs("plots/fNP_unc.pdf");
 
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
      h_fnp2d->SetBinContent(i_cos+1, i_pt+1, h_fnp->GetBinContent(i_pt+1));
      h_fnp2d->SetBinError(i_cos+1, i_pt+1, h_fnp->GetBinError(i_pt+1));
    }
  }

  c->Destructor();

  // scale fractions down from percentage
  h_fnp->Scale(1./100.);
  h_fnp2d->Scale(1./100.);

  TFile *fout = new TFile("files/NPFrac.root", "recreate");
  h_fnp->SetName("fnp_unc");
  h_fnp->Write();
  h_fnp2d->SetName("h_fnp");
  h_fnp2d->Write();
  fout->Close();

}
