#import "/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/ptbins.C"
#import "mbins.C"

void store_SB()
{
  double mB_min[] = {2.92, 3.0, 3.21};
  double mB_max[] = {2.95, 3.2, 3.28};

  double m_SR = 0.5*(mB_max[1]+mB_min[1]);
  double dm_SR = mB_max[1]-mB_min[1];
    
  const int npar = 2;
  double pars[npar][nPtBins], epars[npar][nPtBins];
  double pt[nPtBins], ept[nPtBins];
  
  // get the shape parameters
  TFile *fin_SB = new TFile("files/ltfitres2d.root");
  TFitResult* fitres_s = (TFitResult*)fin_SB->Get("fitres");
  double tnp = fitres_s->Parameter(nmBins*1);
  double etnp = fitres_s->ParError(nmBins*1);
  fin_SB->Close();

  // set values for all pT bins - constant tnp1
  for(int i = 0; i < nPtBins; i++) {
    pt[i] = 0.5*(ptBins[i+1]+ptBins[i]);
    ept[i] = 0.5*(ptBins[i+1]-ptBins[i]);
    
    pars[1][i] = tnp;
    epars[1][i] = etnp;
  }

  // get N from exp fit
  TFitResult *fitres = new TFitResult();
  TFile *fin = new TFile("files/ltfitres_N.root");
  for(int j = 0; j < nPtBins; j++) {
    fin->GetObject(Form("fitres_%.0f", ptBins[j]), fitres);
    double dpt = ptBins[j+1]-ptBins[j];

    double ln = 1e4;
    double p1 = fitres->Parameter(0);
    double p2 = fitres->Parameter(1);
    double sigp1 = fitres->ParError(0);
    double sigp2 = fitres->ParError(1);
    double cov = fitres->GetCovarianceMatrix()[0][1];

    //evaluate N_NP
    double f_val = p1*exp(-m_SR/p2);
    pars[0][j] = dm_SR*f_val;

    // propagate uncertainty
    double devp1 = (p1+sigp1/ln)*exp(-m_SR/p2);
    double devp2 = p1*exp(-m_SR/(p2+sigp2/ln));
    double dp_1 = (devp1-f_val)/(sigp1/ln);
    double dp_2 = (devp2-f_val)/(sigp2/ln);
    double f_err = pow(dp_1*sigp1,2) + pow(dp_2*sigp2,2) + 2*dp_1*dp_2*cov;
    f_err = sqrt(f_err);

    epars[0][j] = dm_SR*f_err;
  }
  fin->Close();

  // store the resulting parameters vs pT
  string parsav[] = {"N_NP", "t_NP"};
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
