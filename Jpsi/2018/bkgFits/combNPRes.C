void combNPRes()
{
  string data;
  double aux;
  const int npar = 10;
  double par[3][npar], epar[3][npar];
  string fname[] = {"", "_hpt", "_hpt_fn"};
  string par_n[] = {"N per 1 GeV", "f (\\%)", "$\\mu$ (MeV)", "$\\sigma_1$ (MeV)", "$\\sigma_2$ (MeV)", "n", "$\\alpha$", "$f_G$ (\\%)", "$\\sigma_G$ (MeV)"};
  double mult[] = {1, 100, 1000, 1000, 1000, 1, 1, 100, 1000, 1.};
  double normv[] = {1, 1, 1, 1, 1, 2, 1, 1, 0, 1};
  
  // open NP fit output, get results
  for(int i = 0; i < 3; i++) {
    ifstream ifile;
    ifile.open(Form("text_output/NP%s_fit.txt", fname[i].c_str()));
    for(int j = 0; j < npar; j++) {
      ifile >> par[i][j]  >> epar[i][j];
      par[i][j] *= mult[j];
      epar[i][j] *= mult[j];
    }
    ifile.close();
  }

  // output fit parameters as a table
  ofstream ftexF;
  ftexF.open("text_output/mfit_MC_NP.tex");
  ftexF << "\\begin{tabular}{c||c|c|c}\n";
  ftexF << "Parameter & Full $\\pt$ & $\\pt>50$ GeV, G fixed & $\\pt>50$ GeV, G+n fixed \\\\\n";
  ftexF << "\\hline\n";

  for(int i = 0; i < npar-1; i++) {
    // first do the 9 parameters
    ftexF << par_n[i];
    for(int j = 0; j < 3; j++) {
      double val = par[j][i], unc = epar[j][i];
      int p_norm = 1.;
      if(unc > 0) {
	if(unc < 1 ) 
	  p_norm = ceil(-log10(unc))+1;	
	ftexF << " & " << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc;
      }
      else {
	p_norm = normv[i];	
	ftexF << " & " << setprecision(p_norm) << fixed << val;
      }
    }
    ftexF << "\\\\\n";
  }
  ftexF << "\\hline\n";
  // then do the chi^2/ndf
  ftexF << "$\\chi^2$/ndf";
  for(int j = 0; j < 3; j++) {
    ftexF << " & " << setprecision(0) << fixed << par[j][npar-1] << "/" << epar[j][npar-1];
  }
  ftexF << "\\\\\n";
  ftexF << "\\end{tabular}\n";
  ftexF.close();

  // open NPR fit output, get results
  for(int i = 0; i < 3; i++) {
    ifstream ifile;
    ifile.open(Form("text_output/NPR%s_fit.txt", fname[i].c_str()));
    for(int j = 0; j < npar; j++) {
      ifile >> par[i][j]  >> epar[i][j];
      par[i][j] *= mult[j];
      epar[i][j] *= mult[j];
    }
    ifile.close();
  }

  // output fit parameters as a table
  ofstream ftexNP;
  ftexNP.open("text_output/mfit_MC_NPR.tex");
  ftexNP << "\\begin{tabular}{c||c|c|c}\n";
  ftexNP << "Parameter & Full $\\pt$ & $\\pt>50$ GeV, G fixed & $\\pt>50$ GeV, G+n fixed \\\\\n";
  ftexNP << "\\hline\n";

  for(int i = 0; i < npar-1; i++) {
    // first do the 9 parameters
    ftexNP << par_n[i];
    for(int j = 0; j < 3; j++) {
      double val = par[j][i], unc = epar[j][i];
      int p_norm = 1.;
      if(unc > 0) {
	if(unc < 1 ) 
	  p_norm = ceil(-log10(unc))+1;	
	ftexNP << " & " << setprecision(p_norm) << fixed << val << " $\\pm$ " << unc;
      }
      else {
	p_norm = normv[i];	
	ftexNP << " & " << setprecision(p_norm) << fixed << val;
      }
    }
    ftexNP << "\\\\\n";
  }
  ftexNP << "\\hline\n";
  // then do the chi^2/ndf
  ftexNP << "$\\chi^2$/ndf";
  for(int j = 0; j < 3; j++) {
    ftexNP << " & " << setprecision(0) << fixed << par[j][npar-1] << "/" << epar[j][npar-1];
  }
  ftexNP << "\\\\\n";
  ftexNP << "\\end{tabular}\n";
  ftexNP.close();
 
  
}
