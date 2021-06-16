void plotLtPars()
{
  // read the fit results
  ifstream ifile;
  string data;
  int pt_bins = 7;
  double pt_min[pt_bins], pt_max[pt_bins], pt_avg[pt_bins], pt_err[pt_bins];
  double par[7][pt_bins], epar[7][pt_bins];
  double chis[pt_bins], ndf[pt_bins], fNP[pt_bins], chiN[pt_bins], zero[pt_bins];

  // read the fit results
  ifile.open("text_output/lt_fit.txt");
  getline(ifile, data);

  double mults[] = {1./10., 1./10., 100., 1e3, 1e3, 1e3, 1e3};
  for(int i = 0; i < pt_bins; i++) {
    ifile >> pt_min[i] >> pt_max[i];
    for(int ip = 0; ip < 7; ip++) {
      ifile >> par[ip][i] >> epar[ip][i];
      par[ip][i] *= mults[ip];
      epar[ip][i] *= mults[ip];
    }
    ifile >> chis[i] >> ndf[i] >> fNP[i];
    fNP[i]*=100;

    pt_avg[i] = 0.5*(pt_max[i]+pt_min[i]);
    pt_err[i] = 0.5*(pt_max[i]-pt_min[i]);
    chiN[i] = chis[i]/ndf[i];
    zero[i] = 0;
  }
  ifile.close();

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);

  double ax_min[] = {1e2, 1e3, 0,  0, 0,  0,  300};
  double ax_max[] = {1e5, 5e5, 20, 0, 20, 40, 400};
  string par_tit[] = {"N_PR", "N_NP", "f", "mu", "sigma1", "sigma2", "lambda"};
  string par_name[] = {"N_{PR}", "N_{NP}", "f", "#mu", "#sigma_{1}", "#sigma_{2}", "#lambda"};
  string par_unit[] = {" per 5 GeV", " per 5 GeV", " (%)", " (#mum)", " (#mum)", " (#mum)", " (#mum)"};
  
  TGraphErrors **g_par = new TGraphErrors*[7];
  for(int i = 0; i < 7; i++) {
    if(i==3) continue;
    g_par[i] = new TGraphErrors(pt_bins, pt_avg, par[i], pt_err, epar[i]);

    if(i < 2) c->SetLogy();
    else c->SetLogy(0);

    TH1F *fp = c->DrawFrame(pt_min[0]-5, ax_min[i], pt_max[pt_bins-1]+5, ax_max[i]);
    fp->SetXTitle("p_{T} (GeV)");
    fp->SetYTitle(Form("%s%s", par_name[i].c_str(), par_unit[i].c_str()));
    fp->GetYaxis()->SetTitleOffset(1.5);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("%s vs p_{T}", par_name[i].c_str()));

    if(i == 2) {
      TF1 *fc = new TF1("fc", "[0]", pt_min[0], pt_max[pt_bins-1]);
      //g_par[i]->Fit(fc);
    }
    
    g_par[i]->SetMarkerStyle(20);
    g_par[i]->SetMarkerSize(.75);
    g_par[i]->SetMarkerColor(kBlack);
    g_par[i]->SetLineColor(kBlack);
    g_par[i]->Draw("psame");


    c->SaveAs(Form("plots/lifetime/par_%s.pdf", par_tit[i].c_str()));
    c->Clear();
  }

  TGraphErrors *g_chi = new TGraphErrors(pt_bins, pt_avg, chiN, pt_err, zero);

  TH1F *fchi = c->DrawFrame(pt_min[0]-5, 0, pt_max[pt_bins-1]+5, 35);
  fchi->SetXTitle("p_{T} (GeV)");
  fchi->SetYTitle("#chi^{2}/ndf");
  fchi->GetYaxis()->SetTitleOffset(1.5);
  fchi->GetYaxis()->SetLabelOffset(0.01);
  fchi->SetTitle("#chi^{2}/ndf vs p_{T}");
  
  g_chi->SetMarkerStyle(20);
  g_chi->SetMarkerSize(.75);
  g_chi->SetMarkerColor(kBlack);
  g_chi->SetLineColor(kBlack);
  g_chi->Draw("psame");

  TF1 *fc = new TF1("fc", "[0]", pt_min[0]-5, pt_max[pt_bins-1]+5);
  fc->SetParameter(0, 1.);
  fc->SetLineColor(kBlue);
  fc->SetLineStyle(kDashed);
  fc->Draw("lsame");
  
  c->SaveAs(Form("plots/lifetime/par_chiN.pdf"));
  c->Clear();
  
  TGraphErrors *g_frac = new TGraphErrors(pt_bins, pt_avg, fNP, pt_err, zero);

  TH1F *ffr = c->DrawFrame(pt_min[0]-5, 0, pt_max[pt_bins-1]+5, 50);
  ffr->SetXTitle("p_{T} (GeV)");
  ffr->SetYTitle("f_{NP} (%)");
  ffr->GetYaxis()->SetTitleOffset(1.3);
  ffr->GetYaxis()->SetLabelOffset(0.01);
  ffr->SetTitle("f_{NP} vs p_{T}");
  
  g_frac->SetMarkerStyle(20);
  g_frac->SetMarkerSize(.75);
  g_frac->SetMarkerColor(kBlack);
  g_frac->SetLineColor(kBlack);
  g_frac->Draw("psame");
  
  c->SaveAs(Form("plots/lifetime/par_fNP.pdf"));
  c->Clear();

  // do the sigma plots to compare with rms from previous ANs
  double sig_avg[pt_bins], esig_avg[pt_bins], rms[pt_bins], erms[pt_bins], sig_rat[pt_bins], esig_rat[pt_bins];
  for(int i_p = 0; i_p < pt_bins; i_p++){
    double sig1_s = pow(par[4][i_p], 2);
    double sig2_s = pow(par[5][i_p], 2);
    sig_avg[i_p] = sqrt(sig1_s+sig2_s);

    double esig1_s = pow(epar[4][i_p], 2);
    double esig2_s = pow(epar[5][i_p], 2);
    esig_avg[i_p] = sqrt(sig1_s*esig1_s+sig2_s*esig2_s)/sig_avg[i_p];

    rms[i_p] = sig_avg[i_p] * pt_avg[i_p]/3.097;
    erms[i_p] = esig_avg[i_p] * pt_avg[i_p]/3.097;

    sig_rat[i_p] = par[5][i_p]/par[4][i_p];
    esig_rat[i_p] = sig_rat[i_p] * sqrt(esig1_s/sig1_s + esig2_s/sig2_s);
  }

  TGraphErrors *g_sig = new TGraphErrors(pt_bins, pt_avg, sig_avg, pt_err, esig_avg);

  TH1F *fsig = c->DrawFrame(pt_min[0]-5, 0., pt_max[pt_bins-1]+5, 60);
  fsig->SetXTitle("p_{T} (GeV)");
  fsig->SetYTitle("#sigma (#mum)");
  fsig->GetYaxis()->SetTitleOffset(1.3);
  fsig->GetYaxis()->SetLabelOffset(0.01);
  fsig->SetTitle("#sigma vs p_{T}");
  
  g_sig->SetMarkerStyle(20);
  g_sig->SetMarkerSize(.75);
  g_sig->SetMarkerColor(kBlack);
  g_sig->SetLineColor(kBlack);
  g_sig->Draw("psame");

  c->SaveAs(Form("plots/lifetime/par_sigma.pdf"));
  c->Clear();

  TGraphErrors *g_rms = new TGraphErrors(pt_bins, pt_avg, rms, pt_err, erms);

  TH1F *frms = c->DrawFrame(pt_min[0]-5, 0., pt_max[pt_bins-1]+5, 800);
  frms->SetXTitle("p_{T} (GeV)");
  frms->SetYTitle("rms (Lxy) (#mum)");
  frms->GetYaxis()->SetTitleOffset(1.3);
  frms->GetYaxis()->SetLabelOffset(0.01);
  frms->SetTitle("rms (Lxy) vs p_{T}");
  
  g_rms->SetMarkerStyle(20);
  g_rms->SetMarkerSize(.75);
  g_rms->SetMarkerColor(kBlack);
  g_rms->SetLineColor(kBlack);
  g_rms->Draw("psame");
  
  c->SaveAs(Form("plots/lifetime/par_rms.pdf"));
  c->Clear();

  TGraphErrors *g_rat = new TGraphErrors(pt_bins, pt_avg, sig_rat, pt_err, esig_rat);

  TH1F *frat = c->DrawFrame(pt_min[0]-5, 1., pt_max[pt_bins-1]+5, 3.);
  frat->SetXTitle("p_{T} (GeV)");
  frat->SetYTitle("#sigma_{2}/#sigma_{1}");
  frat->GetYaxis()->SetTitleOffset(1.3);
  frat->GetYaxis()->SetLabelOffset(0.01);
  frat->SetTitle("#sigma_{2}/#sigma_{1} vs p_{T}");
  
  g_rat->SetMarkerStyle(20);
  g_rat->SetMarkerSize(.75);
  g_rat->SetMarkerColor(kBlack);
  g_rat->SetLineColor(kBlack);
  g_rat->Draw("psame");

  TF1 *fl = new TF1("fl", "[0]*x+[1]", pt_min[0]-5, pt_max[pt_bins-1]+5);
  fl->SetParameters(5e-2, 2.);
  fl->SetLineColor(kBlue);
  fl->SetLineStyle(kDashed);
  g_rat->Fit(fl);
  
  c->SaveAs(Form("plots/lifetime/par_sigRat.pdf"));
  c->Clear();
  c->Destructor();
}
