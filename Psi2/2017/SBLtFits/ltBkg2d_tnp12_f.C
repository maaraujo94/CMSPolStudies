#import "mbins.C"

// define negative exponential only for positive x
double pos_exp(double x, double f1, double ld1, double ld2)
{
  if(x > 0) return f1*exp(-x/ld1) + (1.-f1)*exp(-x/ld2);
  else return 0;
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

// define final fit function (2d)
// parameters: N_NP (exponential), f1 (free), tnp1 (constant), tnp2 (constant)
double func_sum(double *xx, double *pp)
{
  // get m bin
  double lt = xx[0], m = xx[1];
  int m_bin = 100;
  for(int i = 0; i < nmBins; i++)
    if(m_min[i] < m && m_max[i] > m)
      m_bin = i;
  if (m_bin > nmBins){
    TF2::RejectPoint();
    return 0;
  }
  double d_m = m_max[m_bin]-m_min[m_bin];
  
  double N_NP = pp[m_bin]; // free
  double f1 = pp[nmBins] + m*pp[nmBins+1]; // linear
  double tnp1 = pp[2*nmBins]; // constant
  double tnp2 = pp[3*nmBins]; // constant

  return N_NP * pos_exp(lt, f1, tnp1, tnp2);
}

// define final fit function (1d)
double sum_1d(double *xx, double *pp)
{
  double lt = xx[0];
  double N_NP = pp[0], f1 = pp[1], tnp1 = pp[2], tnp2 = pp[3];

  return N_NP * pos_exp(lt, f1, tnp1, tnp2);
}

void ltBkg2d_tnp12_f()
{
  // prepare binning and histograms for plots
  TH1D **h_d1d = new TH1D*[nmBins];
  TFile *fin = new TFile("files/ltStore.root");
  for(int i = 0; i < nmBins; i++) {
    fin->GetObject(Form("ltH_%s", lbl[i].c_str()), h_d1d[i]);
    h_d1d[i]->SetDirectory(0);
  }
  fin->Close();

  int tbins = h_d1d[0]->GetNbinsX();
  double lowt = h_d1d[0]->GetXaxis()->GetBinLowEdge(1);
  double hit = h_d1d[0]->GetXaxis()->GetBinUpEdge(tbins);
  double wbin = (hit-lowt)/(double)tbins;

  //need to define m bin edge array for the 2d histo
  double mBins[nmBins+2];
  for(int i = 0; i < 4; i++) mBins[i] = m_min[i];
  mBins[4] = m_max[3];
  mBins[5] = m_min[4];
  for(int i = 0; i < 4; i++) mBins[6+i] = m_max[i+4];

  // get 2d histo from the sb 1d histos
  TH2D *h_d2d = new TH2D("ltH", "2017 data c#tau (sidebands)", tbins, lowt, hit, nmBins+1, mBins);
  for(int i_x = 0; i_x < tbins; i_x++) {
    for(int i_y = 0; i_y < 4; i_y++) {
      h_d2d->SetBinContent(i_x+1, i_y+1, h_d1d[i_y]->GetBinContent(i_x+1));
      h_d2d->SetBinError(i_x+1, i_y+1, h_d1d[i_y]->GetBinError(i_x+1));
    }
    for(int i_y = 4; i_y < nmBins; i_y++) {
      h_d2d->SetBinContent(i_x+1, i_y+2, h_d1d[i_y]->GetBinContent(i_x+1));
      h_d2d->SetBinError(i_x+1, i_y+2, h_d1d[i_y]->GetBinError(i_x+1));
    }
  }
  
  // define aux vals for plotting
  double pr_lim = 0.05;
  double np_lim = 0.1;
  double lowPlot = -0.1;

  // the fit function
  TF2 *fitS = new TF2("fitS", func_sum, np_lim, hit, m_min[0], m_max[nmBins-1], 4*nmBins, 2);
  string par_n[] =  {"N_NP", "f1",  "t_NP1", "t_NP2"};
  double par_v[] =  {1.e6, 3.5,  3.5e-1,  5.5e-2};
  double par2_v[] = {1.,   -0.8,   1.,      1.};

  // define the parameters
  for(int i = 0; i < nmBins; i++) {
    for(int j = 0; j < 4; j++) {
      // name and initialize all parameters
      fitS->SetParName(j*nmBins+i, Form("%s_%d", par_n[j].c_str(), i));
      fitS->SetParameter(j*nmBins+i, par_v[j]);

      // setting the scale parameter
      if(j==0) fitS->SetParameter(j*nmBins+i, h_d1d[i]->GetMaximum()/3.);
      // setting the linear parameter
      else if(j==1 && i == 1) fitS->SetParameter(j*nmBins+i, par2_v[j]);
      else if(j==1 && i > 1) fitS->FixParameter(j*nmBins+i, par_v[j]);
      // setting the constant parameters t_NP
      else if(j>1 && i > 0) fitS->FixParameter(j*nmBins+i, par_v[j]);
    }
  }
  // fit the 2d function to the lifetime:pT map
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.13);
  TFitResultPtr fitres = h_d2d->Fit("fitS", "SR");

  // tf1 for plotting in the 1D bins
  TF1 *f_1d = new TF1("f_1d", sum_1d, np_lim, hit, 4);
  f_1d->SetParNames("N_NP", "f1", "tnp1", "tnp2");

  double m_val[nmBins], m_err[nmBins];
  double pars[4][nmBins], epars[4][nmBins];
 
  // cycle over all m bins
  for(int i_m = 0; i_m < nmBins; i_m++) {
    m_val[i_m] = 0.5*(m_max[i_m]+m_min[i_m]);
    m_err[i_m] = 0.5*(m_max[i_m]-m_min[i_m]);

    // storing all parameters - N_NP, f1, tnp_1, tnp_2
    for(int j = 0; j < 4; j++) {
      if(j ==0) { // N is free
	pars[j][i_m] = fitS->GetParameter(i_m);
	epars[j][i_m] = fitS->GetParError(i_m);
      }
      else if(j== 1) { // f is linear
	pars[j][i_m] = fitS->GetParameter(j*nmBins) + m_val[i_m] * fitS->GetParameter(j*nmBins+1);
	double ln = 1e4;
	double sigp1 = fitS->GetParError(j*nmBins);
	double sigp2 = fitS->GetParError(j*nmBins+1);
	double cov = fitres->GetCovarianceMatrix()[j*nmBins][j*nmBins+1];

	double devp1 = fitS->GetParameter(j*nmBins)+sigp1/ln +  m_val[i_m] * fitS->GetParameter(j*nmBins+1);
	double devp2 = fitS->GetParameter(j*nmBins) +  m_val[i_m] * (fitS->GetParameter(j*nmBins+1)+sigp2/ln);

	double dp_1 = (devp1-pars[j][i_m])/(sigp1/ln);
	double dp_2 = (devp2-pars[j][i_m])/(sigp2/ln);
	
	epars[j][i_m] = pow(dp_1*sigp1,2) + pow(dp_2*sigp2,2) + 2*dp_1*dp_2*cov;
	epars[j][i_m] = sqrt(epars[j][i_m]);
      }
      else{
	pars[j][i_m] = fitS->GetParameter(j*nmBins); // constant
	epars[j][i_m] = fitS->GetParError(j*nmBins); // constant
      }
    }

    // initializing f_1d and plotting
    f_1d->SetParameters(pars[0][i_m],
			pars[1][i_m],
			pars[2][i_m],
			pars[3][i_m]);
    //f_1d->Eval(0.2);

    c->SetLogy();
    c->SetTopMargin(0.015);
     
    h_d1d[i_m]->SetMaximum(h_d1d[i_m]->GetMaximum()*1.2);
    h_d1d[i_m]->SetMinimum(h_d1d[i_m]->GetMaximum()*5e-3);

    TH1F *fh = c->DrawFrame(lowPlot, h_d1d[i_m]->GetMinimum(), hit, h_d1d[i_m]->GetMaximum());
    fh->SetXTitle("c#tau (mm)");
    fh->SetYTitle(Form("Events per %.0f #mum", wbin*1000.));
    fh->GetYaxis()->SetTitleOffset(1.8);
    fh->GetYaxis()->SetLabelOffset(0.01);
    fh->GetXaxis()->SetLabelOffset(0.015);
    fh->GetXaxis()->SetTitleOffset(1.3);

    h_d1d[i_m]->SetMarkerStyle(20);
    h_d1d[i_m]->SetMarkerColor(kBlack);
    h_d1d[i_m]->SetLineColor(kBlack);
    h_d1d[i_m]->SetMarkerSize(0.75);
    h_d1d[i_m]->Draw("error same");
    
    f_1d->SetLineColor(kBlue);
    f_1d->Draw("lsame");

    // plot each exponential
    TF1 *fexp1 = new TF1("fexp1", "[0]*[1]*exp(-x/[2])", np_lim, hit);
    fexp1->SetParameters(pars[0][i_m], pars[1][i_m], pars[2][i_m]);
    fexp1->SetLineColor(kGreen);
    fexp1->SetLineStyle(kDashDotted);
    fexp1->Draw("lsame");
    TF1 *fexp2 = new TF1("fexp2", "[0]*(1.-[1])*exp(-x/[2])", np_lim, hit);
    fexp2->SetParameters(pars[0][i_m], pars[1][i_m], pars[3][i_m]);
    fexp2->SetLineColor(kRed);
    fexp2->SetLineStyle(kDashDotted);
    fexp2->Draw("lsame");

    // aux lines for the 2.5 sigma and 4 sigma limits
    TLine *lsig1 = new TLine(-pr_lim, h_d1d[i_m]->GetMinimum(), -pr_lim, h_d1d[i_m]->GetMaximum());
    lsig1->SetLineStyle(kDashed);
    lsig1->Draw("lsame");
    TLine *lsig2 = new TLine(pr_lim, h_d1d[i_m]->GetMinimum(), pr_lim, h_d1d[i_m]->GetMaximum());
    lsig2->SetLineStyle(kDashed);
    lsig2->Draw("lsame");
    TLine *lsig3 = new TLine(np_lim, h_d1d[i_m]->GetMinimum(), np_lim, h_d1d[i_m]->GetMaximum());
    lsig3->SetLineStyle(kDashed);
    lsig3->Draw("lsame");

    TLatex lc;
    lc.SetTextSize(0.03);

    // draw mass bin
    double xp = getPos(lowt, hit, 0.35, 0);
    double yp = getPos(h_d1d[i_m]->GetMinimum(), h_d1d[i_m]->GetMaximum(), 0.9, 1);
    if(i_m < 4)
      lc.DrawLatex(xp, yp, Form("#bf{%.2f < M(#mu#mu) < %.2f GeV}", m_min[i_m], m_max[i_m]));
    else
      lc.DrawLatex(xp, yp, Form("#bf{%.3f < M(#mu#mu) < %.3f GeV}", m_min[i_m], m_max[i_m]));
    yp = getPos(h_d1d[i_m]->GetMinimum(), h_d1d[i_m]->GetMaximum(), 0.83, 1);
    // draw the chi^2/ndf
    xp = getPos(lowt, hit, 0.35, 0);
    lc.DrawLatex(xp, yp, Form("#bf{global #chi^{2}/ndf = %.0f / %d}", fitS->GetChisquare(), fitS->GetNDF()));

    c->SaveAs(Form("plots/lifetime2d_tnp12_f/fit_%s.pdf", lbl[i_m].c_str()));
    c->Clear();
    
    // calculating pulls
    double tv[tbins], pv[tbins], dv[tbins];
    for(int i_t = 0 ; i_t < tbins; i_t++) {
      tv[i_t] = h_d1d[i_m]->GetBinCenter(i_t+1);
      double fitv = f_1d->Eval(tv[i_t]);
      double datav = h_d1d[i_m]->GetBinContent(i_t+1);
      double datau = h_d1d[i_m]->GetBinError(i_t+1);
      if(datau > 0 && tv[i_t] > np_lim) pv[i_t] = (datav-fitv)/datau;
      else pv[i_t] = 0;
      if(fitv > 0 && tv[i_t] > np_lim ) dv[i_t] = (datav-fitv)/fitv * 100.;
      else dv[i_t] = 0;
    }
    
    // plotting the pulls
    c->SetLogy(0);
    c->SetTopMargin(0.1);

    TH1F *fp = c->DrawFrame(lowPlot, -9, hit, 9);
    fp->SetXTitle("c#tau (mm)");
    fp->SetYTitle("pulls");
    fp->GetYaxis()->SetTitleOffset(1.3);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("Lifetime fit pulls"));
  
    TGraph *g_pull = new TGraph(tbins, tv, pv);
    g_pull->SetLineColor(kBlack);
    g_pull->SetMarkerColor(kBlack);
    g_pull->SetMarkerStyle(20);
    g_pull->Draw("psame");

    // aux lines - pull=0 and sigma limits
    TF1 *fcons = new TF1("fcons", "[0]", lowPlot, hit);
    fcons->SetParameter(0, 0);
    fcons->SetLineStyle(kDashed);
    fcons->SetLineColor(kBlack);
    fcons->Draw("lsame");

    TLine *psig1 = new TLine(-pr_lim, -9, -pr_lim, 9);
    psig1->SetLineStyle(kDashed);
    psig1->Draw("lsame");
    TLine *psig2 = new TLine(pr_lim, -9, pr_lim, 9);
    psig2->SetLineStyle(kDashed);
    psig2->Draw("lsame");
    TLine *psig3 = new TLine(np_lim, -9, np_lim, 9);
    psig3->SetLineStyle(kDashed);
    psig3->Draw("lsame");

    TLine *plim1 = new TLine(lowPlot, -5, hit, -5);
    plim1->SetLineStyle(kDotted);
    plim1->Draw("lsame");
    TLine *plim2 = new TLine(lowPlot, -3, hit, -3);
    plim2->SetLineStyle(kDotted);
    plim2->Draw("lsame");
    TLine *plim3 = new TLine(lowPlot, 3, hit, 3);
    plim3->SetLineStyle(kDotted);
    plim3->Draw("lsame");
    TLine *plim4 = new TLine(lowPlot, 5, hit, 5);
    plim4->SetLineStyle(kDotted);
    plim4->Draw("lsame");

    xp = getPos(lowt, hit, 0.4, 0);
    yp = getPos(-9, 9, 0.9, 0);
    if(i_m < 4)
      lc.DrawLatex(xp, yp, Form("#bf{%.2f < M(#mu#mu) < %.2f GeV}", m_min[i_m], m_max[i_m]));
    else
      lc.DrawLatex(xp, yp, Form("#bf{%.3f < M(#mu#mu) < %.3f GeV}", m_min[i_m], m_max[i_m]));
    
    c->SaveAs(Form("plots/lifetime2d_tnp12_f/pulls_%s.pdf", lbl[i_m].c_str()));
    c->Clear();

    // plotting the devs
    TH1F *fd = c->DrawFrame(lowPlot, -15, hit, 15);
    fd->SetXTitle("c#tau (mm)");
    fd->SetYTitle("relative difference (%)");
    fd->GetYaxis()->SetTitleOffset(1.3);
    fd->GetYaxis()->SetLabelOffset(0.01);
    fd->SetTitle(Form("Lifetime rel. difference (%s)", lbl[i_m].c_str()));

    TGraph *g_dev = new TGraph(tbins, tv, dv);
    g_dev->SetLineColor(kBlack);
    g_dev->SetMarkerColor(kBlack);
    g_dev->SetMarkerStyle(20);
    g_dev->Draw("psame");
  
    // aux lines - pull = 0 and sigma limits
    fcons->Draw("lsame");

    TLine *dsig1 = new TLine(-pr_lim, -15, -pr_lim, 15);
    dsig1->SetLineStyle(kDashed);
    dsig1->Draw("lsame");
    TLine *dsig2 = new TLine(pr_lim, -15, pr_lim, 15);
    dsig2->SetLineStyle(kDashed);
    dsig2->Draw("lsame");
    TLine *dsig3 = new TLine(np_lim, -15, np_lim, 15);
    dsig3->SetLineStyle(kDashed);
    dsig3->Draw("lsame");

    xp = getPos(lowt, hit, 0.4, 0);
    yp = getPos(-15, 15, 0.9, 0);
    lc.DrawLatex(xp, yp, "#bf{2d fit with fixed t_{NP1,2} and f_{1}}");

    c->SaveAs(Form("plots/lifetime2d_tnp12_f/devs_%s.pdf", lbl[i_m].c_str()));
    c->Clear();
  }
  
  // storing the free parameters
  TFile *fout = new TFile("files/ltfitres2d_tnp12_f.root", "recreate");
  string parlab[] = {"N_NP", "f1", "tnp1", "tnp2"};

  for(int i_p = 0; i_p < 4; i_p++) {
    TGraphErrors *g_par = new TGraphErrors(nmBins, m_val, pars[i_p], m_err, epars[i_p]);
    g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
  }
  
  fitres->SetName("fitres");
  fitres->Write();

  h_d2d->Write();
  fout->Close();

  cout << fitS->GetChisquare() << "/" << fitS->GetNDF() << endl;
  
  c->Destructor();

}
