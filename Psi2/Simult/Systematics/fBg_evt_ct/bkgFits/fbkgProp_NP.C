#import "../../../ptbins.C"

// set mass limits for integral fit
double x[] = {3.43, 3.82};
double y[] = {3.52, 4.0};

// set vector for int values and errors
double dataInt[2];
double dataErr[2];

// function that gives mass integral
double intFunc(double mi, double mf, double A, double b)
{
  return A*b*(exp(-mi/b) - exp(-mf/b));
}

// main fit function
double myFunction(double A, double b)
{
  double intV, factor;
  double chis = 0;
  
  for(int i = 0; i < 2; i++){
    intV = intFunc(x[i], y[i], A, b);
    factor = pow((intV - dataInt[i])/dataErr[i], 2);
    chis += factor;
  }
  return chis;
}

//Minuit-specific "wrapping"
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg)
{
  result=myFunction(par[0], par[1]);
}

// main code for getting the yield ratios (properly scaled already)
void fbkgProp_NP()
{
  // PART 1 : get the pT dists of LSB and RSB data
  TH1D *h_LSB = new TH1D("lsbH_NP", "Run 2 NP Data (NP LSB)", nPtBins, ptBins);
  TH1D *h_RSB = new TH1D("rsbH_NP", "Run 2 NP Data (NP RSB)", nPtBins, ptBins);
  TH1D *h_tot = new TH1D("totH_NP", "Run 2 NP Data (NP SR)", nPtBins, ptBins);

  TFile *fin = new TFile("files/yieldHistos.root");
  h_LSB = (TH1D*)fin->Get("lsbH_NP");
  h_LSB->SetDirectory(0);
  h_RSB = (TH1D*)fin->Get("rsbH_NP");
  h_RSB->SetDirectory(0);
  h_tot = (TH1D*)fin->Get("totH_NP");
  h_tot->SetDirectory(0);
  fin->Close();
    
  // PART 2 : do the fitting bin by bin to get N_SR from N_RSB and N_LSB
  TFitter* fit = new TFitter(2);
  fit->SetFCN(minuitFunction);

  TH1D *h_SR = new TH1D("srH", "Run 2 Data (NP SR)", nPtBins, ptBins);
  
  // cycle over each pT bin
  double pt[nPtBins], ept[nPtBins], par_A[nPtBins], par_b[nPtBins], epar_A[nPtBins], epar_b[nPtBins];
  
  for(int i = 0; i < nPtBins; i++) {
    // get the pt values
    pt[i] = 0.5*(ptBins[i+1]+ptBins[i]);
    ept[i] = 0.5*(ptBins[i+1]-ptBins[i]);
    
    // fill the datasigma vector
    dataInt[0] = h_LSB->GetBinContent(i+1);
    dataInt[1] = h_RSB->GetBinContent(i+1);
    dataErr[0] = h_LSB->GetBinError(i+1);
    dataErr[1] = h_RSB->GetBinError(i+1);

    // use existing results to set initial values
    fit->SetParameter(0, "A", dataInt[0]*1e4, dataInt[0]*1e2, 0, 0);
    fit->SetParameter(1, "b", 0.3, 0.1, 0, 0);

    fit->ExecuteCommand("MIGRAD",0,0);

    // get the fit parameters
    par_A[i] = fit->GetParameter(0);
    par_b[i] = fit->GetParameter(1);
    epar_A[i] = fit->GetParError(0);
    epar_b[i] = fit->GetParError(1);

    double N_SR = intFunc(3.57, 3.81, par_A[i], par_b[i]);
    
    // now get the uncertainty by propagation of errors
    double ln = 10000;
    double cov_v = fit->GetCovarianceMatrixElement(0,1);
    
    // deviation calculations
    double dA = (intFunc(3.57, 3.81, par_A[i] + epar_A[i]/ln, par_b[i]) - N_SR)/(epar_A[i]/ln);
    double db = (intFunc(3.57, 3.81, par_A[i], par_b[i] + epar_b[i]/ln) - N_SR)/(epar_b[i]/ln);
    double v_unc = sqrt( pow(dA * epar_A[i], 2) + pow(db * epar_b[i], 2) + 2*dA*db*cov_v );
    
    h_SR->SetBinContent(i+1, N_SR);
    h_SR->SetBinError(i+1, v_unc);
  }
  // store the fit results
  TGraphErrors *g_A = new TGraphErrors(nPtBins, pt, par_A, ept, epar_A);
  TGraphErrors *g_b = new TGraphErrors(nPtBins, pt, par_b, ept, epar_b);
  TFile *fout_fr = new TFile("files/ctFitRes_NP.root", "recreate");
  g_A->SetName("graph_A");
  g_A->Write();
  g_b->SetName("graph_b");
  g_b->Write();
  fout_fr->Close();

  // now we can get the f_SB from yield ratio
  TH1D *f_SB = new TH1D("f_SB", "Run 2 NP Data f_{SB}", nPtBins, ptBins);
  f_SB = (TH1D*)h_SR->Clone("f_SB");
  f_SB->Sumw2();
  f_SB->Divide(h_tot);

  // PART 3 : plotting
  
  // now the plotting
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);
  c->SetLogy(0);

  // scaled and comparing to f_bkg from fit
  TFile *finBg = new TFile("../../../bkgFits/files/bkgFrac_NP.root");
  TH1D* fbg_fit = (TH1D*)finBg->Get("fbkg_unc");
  fbg_fit->SetDirectory(0);
  finBg->Close();

  TH1F *fr2c = c->DrawFrame(15, 0., 105, 100);
  fr2c->SetXTitle("p_{T} (GeV)");
  fr2c->SetYTitle("f_{NPBg} (%)");
  fr2c->GetYaxis()->SetTitleOffset(2.);
  fr2c->GetYaxis()->SetLabelOffset(0.01);

  f_SB->SetStats(0);
  f_SB->Scale(100.);
  f_SB->SetLineColor(kBlue);
  f_SB->SetMarkerColor(kBlue);
  f_SB->SetMarkerStyle(20);
  f_SB->SetMarkerSize(.5);
  f_SB->Draw("error same");
  
  fbg_fit->SetStats(0);
  fbg_fit->Scale(100.);
  fbg_fit->SetLineColor(kBlue);
  fbg_fit->SetMarkerColor(kBlue);
  fbg_fit->SetMarkerStyle(20);
  fbg_fit->SetMarkerSize(.5);
  fbg_fit->SetLineStyle(kDashed);
  fbg_fit->Draw("error same");

  c->SaveAs("plots/fNPSB.pdf");
  c->Clear();  
  c->Destructor();

  // get 2d version accounting for cos dimension
  // get costh binning from the stored data histos
  TFile *infile = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Simult/PR_fit/files/histoStore.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("PRH"), hist);

  // get the binning
  int nBinsX = hist->GetNbinsX();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/nBinsX;

  infile->Close();

  TH2D *h_fbkg2d = new TH2D("h_fbkg2d", "Run 2 f_{NPBg}", nBinsX, minX, maxX, nPtBins, ptBins);
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    // same result for all costh bins
    for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
      h_fbkg2d->SetBinContent(i_cos+1, i_pt+1, f_SB->GetBinContent(i_pt+1));
      h_fbkg2d->SetBinError(i_cos+1, i_pt+1, f_SB->GetBinError(i_pt+1));
    }
  }

  // scale fractions down from percentage
  f_SB->Scale(1./100.);
  h_fbkg2d->Scale(1./100.);

  TFile *fout_f = new TFile("files/bkgFrac_NP.root", "recreate");
  f_SB->SetTitle("Run 2 f_{NPBg} from yield ratio");
  f_SB->SetName("fbkg_unc");
  f_SB->Write();
  h_fbkg2d->SetName("h_fbkg");
  h_fbkg2d->Write();
  fout_f->Close();

}
