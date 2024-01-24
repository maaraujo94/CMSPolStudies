#import "../ptbins.C"
#import "mbins.C"

void store_SB()
{
  double mB_min[] = {3.4, 3.57, 3.82};
  double mB_max[] = {3.52, 3.81, 4.0};
  double m_SR = 0.5*(mB_max[1]+mB_min[1]);
  double dm_SR = mB_max[1]-mB_min[1];
    
  const int npar = 4;
  double pars[npar][nPtBins], epars[npar][nPtBins];
  double pt[nPtBins], ept[nPtBins];
  
  // get the shape parameters
  TFile *fin_SB = new TFile("files/ltfitres2d_tnp12_f.root");
  TFitResult* fitres_s = (TFitResult*)fin_SB->Get("fitres");
  double f1_1 = fitres_s->Parameter(nmBins), f1_2 = fitres_s->Parameter(nmBins+1);
  double tnp1 = fitres_s->Parameter(nmBins*2);
  double tnp2 = fitres_s->Parameter(nmBins*3);
  double ef1_1 = fitres_s->ParError(nmBins), ef1_2 = fitres_s->ParError(nmBins+1);
  double etnp1 = fitres_s->ParError(nmBins*2);
  double etnp2 = fitres_s->ParError(nmBins*3);
  fin_SB->Close();

  // f1 needs error propagation including covariance term
  double ln = 1e4;
  double cov = fitres_s->GetCovarianceMatrix()[nmBins][nmBins+1];
  double devp1 = (f1_1 + ef1_1/ln) + m_SR * f1_2;
  double devp2 = f1_1 + m_SR * (f1_2 + ef1_2/ln);
  double f_val = f1_1 + m_SR*f1_2;
  double dp_1 = (devp1-f_val)/(ef1_1/ln);
  double dp_2 = (devp2-f_val)/(ef1_2/ln);
  double f_err = pow(dp_1*ef1_1,2) + pow(dp_2*ef1_2,2) + 2*dp_1*dp_2*cov;
  f_err = sqrt(f_err);

  // set values for all pT bins - constant tnp1,2; linear f1
  for(int i = 0; i < nPtBins; i++) {
    pt[i] = 0.5*(ptBins[i+1]+ptBins[i]);
    ept[i] = 0.5*(ptBins[i+1]-ptBins[i]);
    
    pars[2][i] = tnp1;
    epars[2][i] = etnp1;
    pars[3][i] = tnp2;
    epars[3][i] = etnp2;

    // f1 needs error propagation including covariance term
    pars[1][i] = f_val;
    epars[1][i] = f_err;
  }

  // get N from exp fit
  TFitResult *fitres = new TFitResult();
  TFile *fin = new TFile("files/ltfitres_N.root");
  for(int j = 0; j < nPtBins; j++) {
    fin->GetObject(Form("fitres_%.0f", ptBins[j]), fitres);
    double dpt = ptBins[j+1]-ptBins[j];

    double p1 = fitres->Parameter(0);
    double p2 = fitres->Parameter(1);
    double sigp1 = fitres->ParError(0);
    double sigp2 = fitres->ParError(1);
    cov = fitres->GetCovarianceMatrix()[0][1];

    //evaluate N_NP
    f_val = p1*exp(-m_SR/p2);
    pars[0][j] = dm_SR*f_val;

    // propagate uncertainty
    devp1 = (p1+sigp1/ln)*exp(-m_SR/p2);
    devp2 = p1*exp(-m_SR/(p2+sigp2/ln));
    dp_1 = (devp1-f_val)/(sigp1/ln);
    dp_2 = (devp2-f_val)/(sigp2/ln);
    f_err = pow(dp_1*sigp1,2) + pow(dp_2*sigp2,2) + 2*dp_1*dp_2*cov;
    f_err = sqrt(f_err);

    epars[0][j] = dm_SR*f_err;
  }
  fin->Close();

  // store the resulting parameters vs pT
  string parsav[] = {"N_NP", "f_1", "t_NP1", "t_NP2"};
  TGraphErrors **g_par = new TGraphErrors*[npar];
  for(int i = 0; i < npar; i++) {
    g_par[i] = new TGraphErrors(nPtBins, pt, pars[i], ept, epars[i]);
    g_par[i]->SetMarkerStyle(20);
    g_par[i]->SetMarkerSize(.75);
    g_par[i]->SetMarkerColor(kBlack);
    g_par[i]->SetLineColor(kBlack);
  }

  TFile *fout = new TFile("files/store_SB.root", "recreate");
  for(int i = 0; i < npar; i++) {
    g_par[i]->SetName(Form("g_%s", parsav[i].c_str()));
    g_par[i]->Write();
  }
  fout->Close();
}
