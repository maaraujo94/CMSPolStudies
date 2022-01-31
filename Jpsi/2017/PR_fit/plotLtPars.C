double sig_m(double f, double sig1, double sig2)
{
  double sig1_s = pow(sig1, 2);
  double sig2_s = pow(sig2, 2);
  f /= 100.;

  return sqrt(f*sig1_s + (1.-f)*sig2_s);
}

double esig_m(double f, double sig1, double sig2, double ef, double esig1, double esig2)
{
  double sig1_s = pow(sig1, 2);
  double sig2_s = pow(sig2, 2);
  f/= 100.;
  ef /= 100.;

  return sqrt(pow((sig1_s-sig2_s)*ef/2, 2) + pow(f * sig1 * esig1, 2) + pow((1.-f) * sig2 * esig2, 2))/sqrt(f*sig1_s + (1.-f)*sig2_s);
}

void plotLtPars()
{
  const int n_p = 7;

  // read the fit results
  ifstream ifile;
  string data;
  int pt_bins = 17;
  double pt_min[pt_bins], pt_max[pt_bins], pt_avg[pt_bins], pt_err[pt_bins];
  double par[7][pt_bins], epar[7][pt_bins];
  double par_f[7][pt_bins], epar_f[7][pt_bins];
  double par_b[7][pt_bins], epar_b[7][pt_bins];
  double chis[3][pt_bins], ndf[3][pt_bins], zero[pt_bins];
  double fNP[3][pt_bins], chiN[3][pt_bins], pChi[3][pt_bins];

  // read the fit results
  ifile.open("text_output/lt_fit.txt");
  getline(ifile, data);

  double mults[] = {1, 1, 100., 1e3, 1e3, 1e3, 1e3};
  for(int i = 0; i < pt_bins; i++) {
    ifile >> pt_min[i] >> pt_max[i];
    for(int ip = 0; ip < n_p; ip++) {
      ifile >> par[ip][i] >> epar[ip][i];
      par[ip][i] *= mults[ip];
      epar[ip][i] *= mults[ip];
    }
    ifile >> chis[0][i] >> ndf[0][i] >> fNP[0][i];
    fNP[0][i]*=100;

    pt_avg[i] = 0.5*(pt_max[i]+pt_min[i]);
    pt_err[i] = 0.5*(pt_max[i]-pt_min[i]);
    chiN[0][i] = chis[0][i]/ndf[0][i];
    pChi[0][i] = TMath::Prob(chis[0][i], ndf[0][i]);
     zero[i] = 0;
  }
  ifile.close();

  // read the fit results - fixed mu = 0
  ifile.open("text_output/lt_fit_mf.txt");
  getline(ifile, data);

  for(int i = 0; i < pt_bins; i++) {
    ifile >> pt_min[i] >> pt_max[i];
    for(int ip = 0; ip < n_p; ip++) {
      ifile >> par_f[ip][i] >> epar_f[ip][i];
      par_f[ip][i] *= mults[ip];
      epar_f[ip][i] *= mults[ip];
    }
    ifile >> chis[1][i] >> ndf[1][i] >> fNP[1][i];
    fNP[1][i]*=100;

    chiN[1][i] = chis[1][i]/ndf[1][i];
    pChi[1][i] = TMath::Prob(chis[1][i], ndf[1][i]);
  }
  ifile.close();

  // read the fit results - fixed mu = 0, f = 0.14
  ifile.open("text_output/lt_fit_bf.txt");
  getline(ifile, data);

  for(int i = 0; i < pt_bins; i++) {
    ifile >> pt_min[i] >> pt_max[i];
    for(int ip = 0; ip < n_p; ip++) {
      ifile >> par_b[ip][i] >> epar_b[ip][i];
      par_b[ip][i] *= mults[ip];
      epar_b[ip][i] *= mults[ip];
    }
    ifile >> chis[2][i] >> ndf[2][i] >> fNP[2][i];
    fNP[2][i]*=100;

    chiN[2][i] = chis[2][i]/ndf[2][i];
    pChi[2][i] = TMath::Prob(chis[2][i], ndf[2][i]);
  }
  ifile.close();

  // aux arrays
  string parlab[] = {"N_PR", "N_NP", "f", "mu", "sig1", "sig2", "t"};
  string partit[] = {"N_{PR}", "N_{NP}", "f", "#mu", "#sigma_{1}", "#sigma_{2}", "t_{NP}"};
  string par_unit[] = {" per 1 GeV", " per 1 GeV", " (%)", " (#mum)", " (#mum)", " (#mum)", " (#mum)"};

  double parmin[] = {4e1, 3e2, 0,  -5., 0,  0,  300};
  double parmax[] = {1e5, 5e5, 100, 5., 20, 40, 400};

  // initialize tgraphs for parameters
  TGraphErrors **g_par = new TGraphErrors*[n_p];
  TGraphErrors **g_par_f = new TGraphErrors*[n_p];
  TGraphErrors **g_par_b = new TGraphErrors*[n_p];
  TFile *fout = new TFile("files/ltfit.root", "recreate");

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);

  for(int i = 0; i < n_p; i++) {

    g_par[i] = new TGraphErrors(pt_bins, pt_avg, par[i], pt_err, epar[i]);
    g_par_f[i] = new TGraphErrors(pt_bins, pt_avg, par_f[i], pt_err, epar_f[i]);
    g_par_b[i] = new TGraphErrors(pt_bins, pt_avg, par_b[i], pt_err, epar_b[i]);

    if(i < 2) c->SetLogy();
    else c->SetLogy(0);

    TH1F *fp = c->DrawFrame(pt_min[0]-5, parmin[i], pt_max[pt_bins-1], parmax[i]);
    fp->SetXTitle("p_{T} (GeV)");
    fp->SetYTitle(Form("%s%s", partit[i].c_str(), par_unit[i].c_str()));
    fp->GetYaxis()->SetTitleOffset(1.5);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("2017 %s", partit[i].c_str()));

    g_par[i]->SetMarkerStyle(20);
    g_par[i]->SetMarkerSize(.75);
    g_par[i]->SetMarkerColor(kBlack);
    g_par[i]->SetLineColor(kBlack);
    g_par[i]->Draw("psame");

    g_par[i]->SetName(Form("fit_%s", parlab[i].c_str()));
    g_par[i]->Write();

    g_par_f[i]->SetMarkerStyle(20);
    g_par_f[i]->SetMarkerSize(.75);
    g_par_f[i]->SetMarkerColor(kBlue);
    g_par_f[i]->SetLineColor(kBlue);
    g_par_f[i]->Draw("psame");

    g_par_f[i]->SetName(Form("fit_mf_%s", parlab[i].c_str()));
    g_par_f[i]->Write();
    
    g_par_b[i]->SetMarkerStyle(20);
    g_par_b[i]->SetMarkerSize(.75);
    g_par_b[i]->SetMarkerColor(kRed);
    g_par_b[i]->SetLineColor(kRed);
    g_par_b[i]->Draw("psame");

    g_par_b[i]->SetName(Form("fit_bf_%s", parlab[i].c_str()));
    g_par_b[i]->Write();

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(g_par[i], "all free", "pl");
    leg->AddEntry(g_par_b[i], "#mu, f fixed", "pl");
    //leg->Draw();

    c->SaveAs(Form("plots/lifetime/par_%s.pdf", parlab[i].c_str()));
    c->Clear();
  }

  TGraphErrors *g_chi = new TGraphErrors(pt_bins, pt_avg, chiN[0], pt_err, zero);
  TGraphErrors *g_chi_f = new TGraphErrors(pt_bins, pt_avg, chiN[1], pt_err, zero);
  TGraphErrors *g_chi_b = new TGraphErrors(pt_bins, pt_avg, chiN[2], pt_err, zero);

  TH1F *fchi = c->DrawFrame(pt_min[0]-5, 0, pt_max[pt_bins-1]+5, 3);
  fchi->SetXTitle("p_{T} (GeV)");
  fchi->SetYTitle("#chi^{2}/ndf");
  fchi->GetYaxis()->SetTitleOffset(1.5);
  fchi->GetYaxis()->SetLabelOffset(0.01);
  fchi->SetTitle("2017 #chi^{2}/ndf");
  
  g_chi->SetMarkerStyle(20);
  g_chi->SetMarkerSize(.75);
  g_chi->SetMarkerColor(kBlack);
  g_chi->SetLineColor(kBlack);
  g_chi->Draw("psame");

  g_chi_f->SetMarkerStyle(20);
  g_chi_f->SetMarkerSize(.75);
  g_chi_f->SetMarkerColor(kBlue);
  g_chi_f->SetLineColor(kBlue);
  g_chi_f->Draw("psame");

  g_chi_b->SetMarkerStyle(20);
  g_chi_b->SetMarkerSize(.75);
  g_chi_b->SetMarkerColor(kRed);
  g_chi_b->SetLineColor(kRed);
  g_chi_b->Draw("psame");

  TF1 *fc = new TF1("fc", "[0]", pt_min[0]-5, pt_max[pt_bins-1]+5);
  fc->SetParameter(0, 1.);
  fc->SetLineColor(kBlue);
  fc->SetLineStyle(kDashed);
  fc->Draw("lsame");
  
  c->SaveAs(Form("plots/lifetime/par_chiN.pdf"));
  c->Clear();
    // plot chi prob
  TGraph *g_chiP   = new TGraph(pt_bins, pt_avg, pChi[0]);
  TGraph *g_chiP_f = new TGraph(pt_bins, pt_avg, pChi[1]);
  TGraph *g_chiP_b = new TGraph(pt_bins, pt_avg, pChi[2]);

  TH1F *fchiP = c->DrawFrame(pt_min[0]-5, 0, pt_max[pt_bins-1]+5, 1);
  fchiP->SetXTitle("p_{T} (GeV)");
  fchiP->SetYTitle("P(#chi^{2})");
  fchiP->GetYaxis()->SetTitleOffset(1.5);
  fchiP->GetYaxis()->SetLabelOffset(0.01);
  fchiP->SetTitle("2017 P(#chi^{2})");
  
  g_chiP->SetMarkerStyle(20);
  //g_chiP->SetMarkerSize(.75);
  g_chiP->SetMarkerColor(kBlack);
  g_chiP->SetLineColor(kBlack);
  g_chiP->Draw("psame");

  g_chiP_f->SetMarkerStyle(20);
  //g_chiP_f->SetMarkerSize(.75);
  g_chiP_f->SetMarkerColor(kBlue);
  g_chiP_f->SetLineColor(kBlue);
  g_chiP_f->Draw("psame");

  g_chiP_b->SetMarkerStyle(20);
  //g_chiP_b->SetMarkerSize(.75);
  g_chiP_b->SetMarkerColor(kRed);
  g_chiP_b->SetLineColor(kRed);
  g_chiP_b->Draw("psame");

  c->SaveAs(Form("plots/lifetime/par_Pchi.pdf"));
  c->Clear();

  TGraphErrors *g_frac = new TGraphErrors(pt_bins, pt_avg, fNP[0], pt_err, zero);
  TGraphErrors *g_frac_f = new TGraphErrors(pt_bins, pt_avg, fNP[1], pt_err, zero);
  TGraphErrors *g_frac_b = new TGraphErrors(pt_bins, pt_avg, fNP[2], pt_err, zero);

  TH1F *ffr = c->DrawFrame(pt_min[0]-5, 0, pt_max[pt_bins-1]+5, 50);
  ffr->SetXTitle("p_{T} (GeV)");
  ffr->SetYTitle("f_{NP} (%)");
  ffr->GetYaxis()->SetTitleOffset(1.3);
  ffr->GetYaxis()->SetLabelOffset(0.01);
  ffr->SetTitle("2017 f_{NP}");
  
  g_frac->SetMarkerStyle(20);
  g_frac->SetMarkerSize(.75);
  g_frac->SetMarkerColor(kBlack);
  g_frac->SetLineColor(kBlack);
  g_frac->Draw("psame");

  g_frac_f->SetMarkerStyle(20);
  g_frac_f->SetMarkerSize(.75);
  g_frac_f->SetMarkerColor(kBlue);
  g_frac_f->SetLineColor(kBlue);
  g_frac_f->Draw("psame");

  g_frac_b->SetMarkerStyle(20);
  g_frac_b->SetMarkerSize(.75);
  g_frac_b->SetMarkerColor(kRed);
  g_frac_b->SetLineColor(kRed);
  g_frac_b->Draw("psame");

  g_frac->SetName("fit_fNP");
  g_frac->Write();
  g_frac_f->SetName("fit_f_fNP");
  g_frac_f->Write();
  g_frac_b->SetName("fit_b_fNP");
  g_frac_b->Write();
  fout->Close();

  c->SaveAs(Form("plots/lifetime/par_fNP.pdf"));
  c->Clear();
    
  // do the sigma plots to compare with rms from previous ANs
  double sig_avg[3][pt_bins], esig_avg[3][pt_bins], rms[3][pt_bins], erms[3][pt_bins], sig_rat[3][pt_bins], esig_rat[3][pt_bins];
  for(int i_p = 0; i_p < pt_bins; i_p++){

    sig_avg[0][i_p] = sig_m(par[2][i_p], par[4][i_p], par[5][i_p]);
    esig_avg[0][i_p] = esig_m(par[2][i_p], par[4][i_p], par[5][i_p], epar[2][i_p], epar[4][i_p], epar[5][i_p]);
			      
    rms[0][i_p] = sig_avg[0][i_p] * pt_avg[i_p]/3.097;
    erms[0][i_p] = esig_avg[0][i_p] * pt_avg[i_p]/3.097;
    
    sig_rat[0][i_p] = par[5][i_p]/par[4][i_p];
    esig_rat[0][i_p] = sig_rat[0][i_p] * sqrt(pow(epar[4][i_p]/par[4][i_p], 2) + pow(epar[5][i_p]/par[5][i_p], 2));
    
    // now for second set of pars
    sig_avg[1][i_p] = sig_m(par_f[2][i_p], par_f[4][i_p], par_f[5][i_p]);
    esig_avg[1][i_p] = esig_m(par_f[2][i_p], par_f[4][i_p], par_f[5][i_p], epar_f[2][i_p], epar_f[4][i_p], epar_f[5][i_p]);
			      
    rms[1][i_p] = sig_avg[1][i_p] * pt_avg[i_p]/3.097;
    erms[1][i_p] = esig_avg[1][i_p] * pt_avg[i_p]/3.097;
    
    sig_rat[1][i_p] = par_f[5][i_p]/par_f[4][i_p];
    esig_rat[1][i_p] = sig_rat[1][i_p] * sqrt(pow(epar_f[4][i_p]/par_f[4][i_p], 2) + pow(epar_f[5][i_p]/par_f[5][i_p], 2));
    
    // now for third set of pars
    sig_avg[2][i_p] = sig_m(par_b[2][i_p], par_b[4][i_p], par_b[5][i_p]);
    esig_avg[2][i_p] = esig_m(par_b[2][i_p], par_b[4][i_p], par_b[5][i_p], epar_b[2][i_p], epar_b[4][i_p], epar_b[5][i_p]);
			      
    rms[2][i_p] = sig_avg[2][i_p] * pt_avg[i_p]/3.097;
    erms[2][i_p] = esig_avg[2][i_p] * pt_avg[i_p]/3.097;
    
    sig_rat[2][i_p] = par_b[5][i_p]/par_b[4][i_p];
    esig_rat[2][i_p] = sig_rat[2][i_p] * sqrt(pow(epar_b[4][i_p]/par_b[4][i_p], 2) + pow(epar_b[5][i_p]/par_b[5][i_p], 2));

  }

  TGraphErrors *g_sig = new TGraphErrors(pt_bins, pt_avg, sig_avg[0], pt_err, esig_avg[0]);
  TGraphErrors *g_sig_f = new TGraphErrors(pt_bins, pt_avg, sig_avg[1], pt_err, esig_avg[1]);
  TGraphErrors *g_sig_b = new TGraphErrors(pt_bins, pt_avg, sig_avg[2], pt_err, esig_avg[2]);

  TH1F *fsig = c->DrawFrame(pt_min[0]-5, 0., pt_max[pt_bins-1]+5, 20);
  fsig->SetXTitle("p_{T} (GeV)");
  fsig->SetYTitle("#sigma (#mum)");
  fsig->GetYaxis()->SetTitleOffset(1.3);
  fsig->GetYaxis()->SetLabelOffset(0.01);
  fsig->SetTitle("2017 #sigma_{m}");
  
  g_sig->SetMarkerStyle(20);
  g_sig->SetMarkerSize(.75);
  g_sig->SetMarkerColor(kBlack);
  g_sig->SetLineColor(kBlack);
  g_sig->Draw("psame");

  g_sig_f->SetMarkerStyle(20);
  g_sig_f->SetMarkerSize(.75);
  g_sig_f->SetMarkerColor(kBlue);
  g_sig_f->SetLineColor(kBlue);
  g_sig_f->Draw("psame");

  g_sig_b->SetMarkerStyle(20);
  g_sig_b->SetMarkerSize(.75);
  g_sig_b->SetMarkerColor(kRed);
  g_sig_b->SetLineColor(kRed);
  g_sig_b->Draw("psame");

  c->SaveAs(Form("plots/lifetime/par_sigma.pdf"));
  c->Clear();

  TGraphErrors *g_rms = new TGraphErrors(pt_bins, pt_avg, rms[0], pt_err, erms[0]);
  TGraphErrors *g_rms_f = new TGraphErrors(pt_bins, pt_avg, rms[1], pt_err, erms[1]);
  TGraphErrors *g_rms_b = new TGraphErrors(pt_bins, pt_avg, rms[2], pt_err, erms[2]);

  TH1F *frms = c->DrawFrame(pt_min[0]-5, 0., pt_max[pt_bins-1]+5, 500);
  frms->SetXTitle("p_{T} (GeV)");
  frms->SetYTitle("rms (Lxy) (#mum)");
  frms->GetYaxis()->SetTitleOffset(1.3);
  frms->GetYaxis()->SetLabelOffset(0.01);
  frms->SetTitle("2017 rms (Lxy)");
  
  g_rms->SetMarkerStyle(20);
  g_rms->SetMarkerSize(.75);
  g_rms->SetMarkerColor(kBlack);
  g_rms->SetLineColor(kBlack);
  g_rms->Draw("psame");

  g_rms_f->SetMarkerStyle(20);
  g_rms_f->SetMarkerSize(.75);
  g_rms_f->SetMarkerColor(kBlue);
  g_rms_f->SetLineColor(kBlue);
  g_rms_f->Draw("psame");

  g_rms_b->SetMarkerStyle(20);
  g_rms_b->SetMarkerSize(.75);
  g_rms_b->SetMarkerColor(kRed);
  g_rms_b->SetLineColor(kRed);
  g_rms_b->Draw("psame");

  c->SaveAs(Form("plots/lifetime/par_rms.pdf"));
  c->Clear();

  TGraphErrors *g_rat = new TGraphErrors(pt_bins, pt_avg, sig_rat[0], pt_err, esig_rat[0]);
  TGraphErrors *g_rat_f = new TGraphErrors(pt_bins, pt_avg, sig_rat[1], pt_err, esig_rat[1]);
  TGraphErrors *g_rat_b = new TGraphErrors(pt_bins, pt_avg, sig_rat[2], pt_err, esig_rat[2]);

  TH1F *frat = c->DrawFrame(pt_min[0]-5, 1., pt_max[pt_bins-1]+5, 3.);
  frat->SetXTitle("p_{T} (GeV)");
  frat->SetYTitle("#sigma_{2}/#sigma_{1}");
  frat->GetYaxis()->SetTitleOffset(1.3);
  frat->GetYaxis()->SetLabelOffset(0.01);
  frat->SetTitle("2017 #sigma_{2}/#sigma_{1}");
  
  g_rat->SetMarkerStyle(20);
  g_rat->SetMarkerSize(.75);
  g_rat->SetMarkerColor(kBlack);
  g_rat->SetLineColor(kBlack);
  g_rat->Draw("psame");

  g_rat_f->SetMarkerStyle(20);
  g_rat_f->SetMarkerSize(.75);
  g_rat_f->SetMarkerColor(kBlue);
  g_rat_f->SetLineColor(kBlue);
  g_rat_f->Draw("psame");
  
  g_rat_b->SetMarkerStyle(20);
  g_rat_b->SetMarkerSize(.75);
  g_rat_b->SetMarkerColor(kRed);
  g_rat_b->SetLineColor(kRed);
  g_rat_b->Draw("psame");
  
  c->SaveAs(Form("plots/lifetime/par_sigRat.pdf"));
  c->Clear();

  c->Destructor();

  // outputting tex file w par results
  ofstream ftex;
  ftex.open(Form("text_output/tfit_res.tex"));
  ftex << "\\begin{tabular}{c||c|c|c|c|c|c|c||c|c}\n";
  ftex << "$\\pt$ (GeV) & $N_{PR}$ & $N_{NP}$ & f (\\%) & $\\mu$ ($\\mu$m) & $\\sigma_1$ ($\\mu$m) & $\\sigma_2$ ($\\mu$m)  & $t_{NP}$ ($\\mu$m) & $f_{NP}$ (\\%) & $\\chi^2$/ndf \\\\\n";
  ftex << "\\hline\n";
  for(int i = 0; i < pt_bins; i++) {
    // pT bin
    ftex << Form("$[%.0f, %.0f]$", pt_min[i], pt_max[i]);
    for(int i_p = 0; i_p < n_p; i_p++) {
      double val = par[i_p][i], unc = epar[i_p][i];
      if (unc > 0) {
	int p_norm = 1.; 
	if(unc < 1 ) 
	  p_norm = ceil(-log10(unc))+1;	
	ftex << " & " <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc;
      }
      else {
	int p_norm = 1.; 
	if(val < 1 ) 
	  p_norm = ceil(-log10(val))+1;	
	ftex << " & " <<  setprecision(p_norm) << fixed << val ;
      }
    }
    ftex << " & " << setprecision(2) << fixed << fNP[0][i];
    ftex << " & " << setprecision(0) << fixed << chis[0][i] << "/" << ndf[0][i];
    ftex <<  "\\\\\n";
  }
  ftex << "\\end{tabular}\n";
  ftex.close();

  // same for fixed m
  ofstream ftex_mf;
  ftex_mf.open(Form("text_output/tfit_mf_res.tex"));
  ftex_mf << "\\begin{tabular}{c||c|c|c|c|c|c|c||c|c}\n";
  ftex_mf << "$\\pt$ (GeV) & $N_{PR}$ & $N_{NP}$ & f (\\%) & $\\mu$ ($\\mu$m) & $\\sigma_1$ ($\\mu$m) & $\\sigma_2$ ($\\mu$m)  & $t_{NP}$ ($\\mu$m) & $f_{NP}$ (\\%) & $\\chi^2$/ndf \\\\\n";
  ftex_mf << "\\hline\n";
  for(int i = 0; i < pt_bins; i++) {
    // pT bin
    ftex_mf << Form("$[%.0f, %.0f]$", pt_min[i], pt_max[i]);
    for(int i_p = 0; i_p < n_p; i_p++) {
      double val = par_f[i_p][i], unc = epar_f[i_p][i];
      if (unc > 0) {
	int p_norm = 1.; 
	if(unc < 1 ) 
	  p_norm = ceil(-log10(unc))+1;	
	ftex_mf << " & " <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc;
      }
      else {
	int p_norm = 1.; 
	if(val < 1 ) 
	  p_norm = ceil(-log10(val))+2;	
	ftex_mf << " & " <<  setprecision(p_norm) << fixed << val ;
      }
    }
    ftex_mf << " & " << setprecision(2) << fixed << fNP[1][i];
    ftex_mf << " & " << setprecision(0) << fixed << chis[1][i] << "/" << ndf[1][i];
    ftex_mf <<  "\\\\\n";
  }
  ftex_mf << "\\end{tabular}\n";
  ftex_mf.close();

  // same for fixed m, f
  ofstream ftex_bf;
  ftex_bf.open(Form("text_output/tfit_bf_res.tex"));
  ftex_bf << "\\begin{tabular}{c||c|c|c|c|c|c|c||c|c}\n";
  ftex_bf << "$\\pt$ (GeV) & $N_{PR}$ & $N_{NP}$ & f (\\%) & $\\mu$ ($\\mu$m) & $\\sigma_1$ ($\\mu$m) & $\\sigma_2$ ($\\mu$m)  & $t_{NP}$ ($\\mu$m) & $f_{NP}$ (\\%) & $\\chi^2$/ndf \\\\\n";
  ftex_bf << "\\hline\n";
  for(int i = 0; i < pt_bins; i++) {
    // pT bin
    ftex_bf << Form("$[%.0f, %.0f]$", pt_min[i], pt_max[i]);
    for(int i_p = 0; i_p < n_p; i_p++) {
      double val = par_b[i_p][i], unc = epar_b[i_p][i];
      if (unc > 0) {
	int p_norm = 1.; 
	if(unc < 1 ) 
	  p_norm = ceil(-log10(unc))+1;	
	ftex_bf << " & " <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc;
      }
      else {
	int p_norm = 1.; 
	if(val < 1 ) 
	  p_norm = ceil(-log10(val))+2;	
	ftex_bf << " & " <<  setprecision(p_norm) << fixed << val ;
      }
    }
    ftex_bf << " & " << setprecision(2) << fixed << fNP[2][i];
    ftex_bf << " & " << setprecision(0) << fixed << chis[2][i] << "/" << ndf[2][i];
    ftex_bf <<  "\\\\\n";
  }
  ftex_bf << "\\end{tabular}\n";
  ftex_bf.close();


}
