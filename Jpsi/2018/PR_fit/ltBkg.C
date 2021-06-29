#import "ltPerPt.C"
#import "ltPerPt_muFix.C"
#import "ltPerPt_bFix.C"
#import "plotLtPars.C"

void ltBkg()
{
  // prepare binning and histograms for plots
  const int nPtBins = 17;
  double ptBins[nPtBins+1];
  int yBins_c[nPtBins+1];
  for(int i = 0; i < 7; i++) ptBins[i] = 25 + 3.*i;
  for(int i = 0; i < 6; i++) ptBins[i+7] = 46 + 5.*i;
  for(int i = 0; i < 3; i++) ptBins[i+13] = 76 + 8.*i;
  for(int i = 0; i < 2; i++) ptBins[i+16] = 100 + 20.*i;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;
  
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

  for(int i = 0; i < nPtBins; i++) {
    ltPerPt(ptBins[i], ptBins[i+1]);
    cout << endl << endl;
    ltPerPt_muFix(ptBins[i], ptBins[i+1]);
    cout << endl << endl;
    ltPerPt_bFix(ptBins[i], ptBins[i+1]);
    cout << endl << endl;
  }
  plotLtPars();
}
