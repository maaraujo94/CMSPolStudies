// macro to propagate the uncertainties to f_NP analytically
// then save it as a 2d histo unc band

int n_pt;

// define negative exponential only for positive x
double pos_exp(double x, double ld)
{
  if(x > 0) return exp(-x/ld);
  else return 0;
}

int ent_v(int i, int i_pt)
{
  if(i == 0) return n_pt+i_pt;
  else if (i < 3) return 2*n_pt + (i-1);
  else return n_pt*(i-1)+2+i_pt;
}

void fnpProp()
{
  // PART 1: get the uncertainty of the f_NP
  
  // get all needed input
  // define aux vals for plotting
  double lowPlot = -0.1;
  double pr_lim = 0.05;
  double lowt = -0.05;
  double hit = 0.5;

  // define the resolution (=PR) function
  TF1 *fres = new TF1("fres", "[0]*([1]*TMath::Gaus(x, [2],[3]) + (1.-[1])*TMath::Gaus(x, [2], [4]))", 5*lowt, 5*hit);

  // define the NP function by convolution
  TF1 *fexp = new TF1("fexp", "pos_exp(x,[0])", 5*lowt, 5*hit);
  TF1Convolution *fcNP = new TF1Convolution(fres, fexp, 5*lowt, 5*hit);
  fcNP->SetRange(5*lowt, 5*hit);
  fcNP->SetNofPointsFFT(1000);
  TF1 *fNP = new TF1("fNP", *fcNP, lowPlot, hit, fcNP->GetNpar());
  fNP->SetParNames("N_NP", "f", "mu", "sig1", "sig2", "lambda");
  
  // get fit parameters - need to know which params are being used
  TFile *inNP = new TFile("files/ltfitres2d.root");
  int n_par = fNP->GetNpar(); 
  TGraphErrors** g_par = new TGraphErrors*[n_par];
  for(int i = 0; i < n_par; i++) {
    inNP->GetObject(Form("fit_%s", fNP->GetParName(i)), g_par[i]);
  }
  // get fitres to get cov in each pT bin
  TFitResult *fitres = new TFitResult();
  inNP->GetObject("fitres", fitres);
  inNP->Close();

  // get pT binning from the fit pars
  n_pt = g_par[0]->GetN();
  double ptBins[n_pt+1];
  for(int i = 0; i < n_pt; i++) {
    ptBins[i] = g_par[0]->GetX()[i]-g_par[0]->GetEX()[i];
  }
  ptBins[n_pt] = g_par[0]->GetX()[n_pt-1]+g_par[0]->GetEX()[n_pt-1];
  
  // prepare lt histograms
  TH1D **h_d1d = new TH1D*[n_pt];
  TFile *fin = new TFile("files/ltStore.root");
  for(int ip = 0; ip < n_pt; ip++) {
    fin->GetObject(Form("ltH%.0f", ptBins[ip]), h_d1d[ip]);
    h_d1d[ip]->SetDirectory(0);
  }
  fin->Close();

  // fNP = integral / evt_all (in prompt region)
  double epsrel = 1e-6, error = 0;

  TH1D *h_fnp = new TH1D("h_fnp", "2018 f_{NP}", n_pt, ptBins);
  double ln = 10000;
  double fit_v[n_par], dpar[n_par];
  double cov[n_par][n_par];

  for(int i_pt = 0; i_pt < n_pt; i_pt++) {
    // N_histo has no uncertainty: just histogram integral
    double min_bin = h_d1d[i_pt]->GetXaxis()->FindBin(-pr_lim+1e-6);
    double max_bin = h_d1d[i_pt]->GetXaxis()->FindBin(pr_lim-1e-6);
    double evt_all = h_d1d[i_pt]->Integral(min_bin, max_bin, "width");

    // define input parameters for integral - varies by pT bin
    for(int i = 0; i < n_par; i++) {
      fit_v[i] = g_par[i]->GetY()[i_pt];
      for(int j = 0; j < n_par; j++) {
	// covariance taken straight from fit, no correction
	int pos_i = ent_v(i,i_pt);
	int pos_j = ent_v(j,i_pt);
	cov[i][j] = fitres->GetCovarianceMatrix()[pos_i][pos_j];
      }
    }
    // integral is a function of all above params
    fNP->SetParameters(fit_v);
    double fv = fNP->IntegralOneDim(-pr_lim, pr_lim, epsrel, 1., error), fe = 0;

    // get the function deviation for each parameter
    for(int i = 0; i < n_par; i++) {
      if(cov[i][i]>0) {
	fNP->SetParameter(i, fit_v[i] + sqrt(cov[i][i])/ln);
	dpar[i] = (fNP->IntegralOneDim(-pr_lim, pr_lim, epsrel, 1., error)-fv)/(sqrt(cov[i][i])/ln);
	fNP->SetParameter(i, fit_v[i]);
      }
      else dpar[i] = 0;
    }
    
    for(int i = 0; i < n_par; i++)  
      for(int j = 0; j < n_par; j++) 
	fe += dpar[i]*dpar[j]*cov[i][j];
    fe = sqrt(fe);

    // fill pT bin
    h_fnp->SetBinContent(i_pt+1, fv/evt_all*100.);
    h_fnp->SetBinError(i_pt+1, fe/evt_all*100.);
  }

  // plotting in pT
  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1F *fr1 = c->DrawFrame(ptBins[0]-5, 0, ptBins[n_pt]+5, 50);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f_{NP} (%)");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("2018 f_{NP}");

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
  infile->GetObject(Form("dataH_ab"), hist);

  // get the binning
  int nBinsX = hist->GetNbinsX();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/nBinsX;

  infile->Close();

  // f_np(pT) but generating 2d map so it's easier to apply uncertainties
  TH2D *h_fnp2d = new TH2D("h_fnp2d", "2018 f_{NP}", nBinsX, minX, maxX, n_pt, ptBins);
  for(int i_pt = 0; i_pt < n_pt; i_pt++) {
    // same result for all costh bins
    for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
      h_fnp2d->SetBinContent(i_cos+1, i_pt+1, h_fnp->GetBinContent(i_pt+1));
      h_fnp2d->SetBinError(i_cos+1, i_pt+1, h_fnp->GetBinError(i_pt+1));
    }
  }

  // plotting the 1d projection into pT
  TH1D* h_fnppt = h_fnp2d->ProjectionY("h_fnppt", 1, 1);

  h_fnppt->SetStats(0);
  h_fnppt->SetMinimum(0);
  h_fnppt->SetMaximum(50);
  h_fnppt->GetXaxis()->SetTitle("p_{T} (GeV)");
  h_fnppt->GetYaxis()->SetTitle("f_{NP} (%)");
  h_fnppt->GetYaxis()->SetTitleOffset(1.3);
  h_fnppt->GetYaxis()->SetLabelOffset(0.01);
  h_fnppt->SetTitle("2018 f_{NP}");
  h_fnppt->SetFillColorAlpha(kBlue, 0.5);
  h_fnppt->Draw("e3");
  h_fnp->Draw("error same");

  c->SaveAs("plots/fNP_band.pdf");
  c->Clear();
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
