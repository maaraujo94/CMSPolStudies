// macro to get the corrected f_bkg in each pt bin
void getfBkg()
{
  // get binning from the stored data histos
  TFile *infile = new TFile("files/histoStore.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("dataH_ab"), hist);
  hist->SetDirectory(0);
  infile->Close();
  
  int nBinsX = hist->GetNbinsX(), nBinsY = hist->GetNbinsY();
  const double *yBins = hist->GetYaxis()->GetXbins()->GetArray();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);
  
  // get fractions from full pT interval
  TFile *ffull = new TFile("../bkgFits/files/ltfit.root");
  TGraph *g_f = (TGraph*)ffull->Get("g_fNP");
  double f_SR = g_f->GetX()[0];
  double f_LSB = g_f->GetX()[1];
  double f_RSB = g_f->GetX()[2];
  ffull->Close();

  // get relative weight f_L
  TFile *ffl = new TFile("files/store_fL.root");
  TGraph *g_fL = (TGraph*)ffl->Get("g_fL");
  double fL = g_fL->GetY()[0];
  ffl->Close();

  // get f_NP in each pT bin
  TFile *fpt = new TFile("files/ltfit.root");
  TGraphErrors *g_fNP = (TGraphErrors*)fpt->Get("fit_b_fNP");
  double *f_NP = g_fNP->GetY();

  // get f_bkg as a function of pT
  TFile *ffbkg = new TFile("files/bkgFrac.root");
  TF1 *f_fbkg = (TF1*)ffbkg->Get("fit_SB");
  ffbkg->Close();
  
  // STEP 1 : calculate f_NP^SB
  double f_SB = fL * f_LSB + (1.-fL) * f_RSB;

  // STEP 2 : get f_NP^SB in each pT bin
  double f_SBpt[3][nBinsY];
  double f_prod = f_SB * f_SR, f_rat = f_SB / f_SR;

  cout << "f_SB = " << f_SB << " ; ratio = " << f_rat << " ; product = " << f_prod << endl;

  for(int i = 0; i < nBinsY; i++) {
    // constant in pT
    f_SBpt[0][i] = f_SB;
    // increasing in pT
    f_SBpt[1][i] = f_rat * f_NP[i]/100.;
    // decreasing in pT
    f_SBpt[2][i] = f_prod / (f_NP[i]/100.);
  }

  // STEP 3 : get the corrected f_bkg * ( 1 - f_SBpt ) in each pT bin
  double *ptv = g_fNP->GetX();
  double f_bkg_o[nBinsY];
  double f_bkg[3][nBinsY];
  for(int ip = 0; ip < nBinsY; ip++) {
    for(int im = 0; im < 3; im++) {
      f_bkg[im][ip] = f_fbkg->Eval(ptv[ip]) * (1.-f_SBpt[im][ip]);
      f_bkg_o[ip] = f_fbkg->Eval(ptv[ip]);
    }
  }

  TGraphErrors **g_fbkg = new TGraphErrors*[3];
  TGraphErrors **g_fNPSB = new TGraphErrors*[3];
  TFile *fout = new TFile("files/bkgFracCorr.root", "recreate");

  for(int i = 0; i < 3; i++){
    g_fbkg[i] = new TGraphErrors(nBinsY, ptv, f_bkg[i], g_fNP->GetEX(), g_fNP->GetEY());
    g_fbkg[i]->SetName(Form("g_fbkg_%d", i));
    g_fbkg[i]->Write();

    g_fNPSB[i] = new TGraphErrors(nBinsY, ptv, f_SBpt[i], g_fNP->GetEX(), g_fNP->GetEY());
    g_fNPSB[i]->SetName(Form("g_fNPSB_%d", i));
    g_fNPSB[i]->Write();
  }
  TGraphErrors* g_fbkg_o = new TGraphErrors(nBinsY, ptv, f_bkg_o, g_fNP->GetEX(), g_fNP->GetEY());
  g_fbkg_o->SetName("g_fbkg_org");
  g_fbkg_o->Write();
  fout->Close();

}
