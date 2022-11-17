#import "ltPerPt.C"
#import "ltPerPt_muFix.C"
#import "ltPerPt_bFix.C"
#import "../ptbins.C"

int do_round(double val)
{
  int valR = (int)val;
  if (val-valR > 0.5) return valR+1;
  else return valR;
}

// macro to run all the lifetime fits in the 3 models
void ltBkg()
{
  // prepare output tables
  ofstream ftable;
  ftable.open("text_output/lt_fit.txt");
  ftable << "pt_min\t pt_max\t N_PR\t eN_PR\t N_NP\t eN_NP\t f\t ef\t mu\t emu\t sigma1\t esigma1\t sigma2\t esigma2\t lambda\t elambda\t chi2\t NDF\t f_NP\n";
  ftable.close();

  ofstream ftable2;
  ftable2.open("text_output/lt_fit_mf.txt");
  ftable2 << "pt_min\t pt_max\t N_PR\t eN_PR\t N_NP\t eN_NP\t f\t ef\t mu\t emu\t sigma1\t esigma1\t sigma2\t esigma2\t lambda\t elambda\t chi2\t NDF\t f_NP\n";
  ftable2.close();

  ofstream ftable3;
  ftable3.open("text_output/lt_fit_bf.txt");
  ftable3 << "pt_min\t pt_max\t N_PR\t eN_PR\t N_NP\t eN_NP\t f\t ef\t mu\t emu\t sigma1\t esigma1\t sigma2\t esigma2\t lambda\t elambda\t chi2\t NDF\t f_NP\n";
  ftable3.close();

  // run first fits - all params free
  for(int i = 0; i < nPtBins; i++) {
    ltPerPt(ptBins[i], ptBins[i+1], i);
    cout << endl;
  }

  // open output table, get mu values
  ifstream ifile;
  string data;
  double aux;
  ifile.open("text_output/lt_fit.txt");
  getline(ifile, data);

  double mu[nPtBins], emu[nPtBins], pt_avg[nPtBins], pt_err[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    ifile >> aux  >> aux;
    for(int ip = 0; ip < 3; ip++) {
      ifile >> aux >> aux;
    }
    ifile >> mu[i] >> emu[i];
    for(int ip = 0; ip < 3; ip++) {
      ifile >> aux >> aux;
    }
    ifile >> aux >> aux >> aux;
    pt_avg[i] = 0.5*(ptBins[i+1]+ptBins[i]);
    pt_err[i] = 0.5*(ptBins[i+1]-ptBins[i]);
  }
  ifile.close();

  // fit to constant
  TGraphErrors *g_mu = new TGraphErrors(nPtBins, pt_avg, mu, pt_err, emu);
  TF1 *fc = new TF1("fc", "[0]", ptBins[0], ptBins[nPtBins]);
  fc->SetParameter(0,0);
  g_mu->Fit(fc);
  double mu_avg = fc->GetParameter(0);
  int nm = ceil(-log10(mu_avg))+2;	
  mu_avg = do_round(mu_avg*pow(10, nm))/pow(10, nm);

  cout << endl << endl << nm << " " << mu_avg << endl << endl;
  
  // run lifetime fit with mu fixed to above value
  for(int i = 0; i < nPtBins; i++) {
    ltPerPt_muFix(ptBins[i], ptBins[i+1], i, mu_avg);
    cout << endl;
  }

  // open output table, get f values
  ifile.open("text_output/lt_fit_mf.txt");
  getline(ifile, data);

  double fv[nPtBins], efv[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    ifile >> aux  >> aux;
    for(int ip = 0; ip < 2; ip++) {
      ifile >> aux >> aux;
    }
    ifile >> fv[i] >> efv[i];
    for(int ip = 0; ip < 4; ip++) {
      ifile >> aux >> aux;
    }
    ifile >> aux >> aux >> aux;
  }
  ifile.close();

  // fit to constant
  TGraphErrors *g_f = new TGraphErrors(nPtBins, pt_avg, fv, pt_err, efv);
  fc->SetParameter(0,0.8);
  g_f->Fit(fc);
  double f_avg = fc->GetParameter(0);
  nm = ceil(-log10(f_avg))+2;	
  f_avg = do_round(f_avg*pow(10, nm))/pow(10, nm);

  cout << endl << endl << nm << " " << f_avg << endl << endl;

  TFile *fout = new TFile("files/ltfitres.root", "recreate");
  fout->Close();
  
  for(int i = 0; i < nPtBins; i++) {
    ltPerPt_bFix(ptBins[i], ptBins[i+1], i, mu_avg, f_avg);
    cout << endl;
  }
}
