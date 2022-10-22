void plotLtParsBoth()
{
  const int n_p = 7;

  // read the fit results - 1d fit
  ifstream ifile;
  string data;
  int pt_bins = 19;
  double pt_min[pt_bins], pt_max[pt_bins], pt_avg[pt_bins], pt_err[pt_bins];
  double par[7][pt_bins], epar[7][pt_bins], fNP[pt_bins], zero[pt_bins];
  double aux;
  
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
    ifile >> aux >> aux >> fNP[i];
    fNP[i]*=100;
    zero[i] = 0;
    
    pt_avg[i] = 0.5*(pt_max[i]+pt_min[i]);
    pt_err[i] = 0.5*(pt_max[i]-pt_min[i]);
  }
  ifile.close();

  // aux arrays
  string parlab[] = {"N_PR", "N_NP", "f", "mu", "sig1", "sig2", "lambda"};
  string partit[] = {"N_{PR}", "N_{NP}", "f", "#mu", "#sigma", "#sigma_{2}", "t_{NP}"};
  string par_unit[] = {" per 1 GeV", " per 1 GeV", " (%)", " (#mum)", " (#mum)", " (#mum)", " (#mum)"};

  double parmin[] = {4e1, 4e2, 0,  -5., 0,  0,  300};
  double parmax[] = {2e5, 8e5, 100, 5., 30, 40, 400};

  // read the fit results - 2d fit
  // initialize tgraphs for parameters
  TGraphErrors **g_par2d = new TGraphErrors*[n_p];
  TGraphErrors *g_fNP2d = new TGraphErrors();
  TFile *fin = new TFile(Form("files/ltfitres2d.root"));
  for(int i_p = 0; i_p < n_p; i_p++) {
    fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par2d[i_p]);
  }
  fin->GetObject("fit_fNP", g_fNP2d);
  fin->Close();

  double *xv = g_par2d[0]->GetX();
  double *xe = g_par2d[0]->GetEX();
  int n = g_par2d[0]->GetN();
  // scale all graphs for plotting
  TGraphErrors **g_par_s = new TGraphErrors*[n_p];
  double mult[] = {1., 1., 100., 1e3, 1e3, 1e3, 1e3};
  for(int i = 0; i < n_p; i++) {
    double *yv = g_par2d[i]->GetY();
    double *ye = g_par2d[i]->GetEY();

    for(int j = 0; j < n; j++) {
      if(i < 2) {
	yv[j] /= (2.*xe[j]);
	ye[j] /= (2.*xe[j]);
      }
      else {
	yv[j] *= mult[i];
	ye[j] *= mult[i];
      }
    }
    g_par_s[i] = new TGraphErrors(n, xv, yv, xe, ye);
  }

  double *yv = g_fNP2d->GetY();
  for(int j = 0; j < n; j++) {
    yv[j] *= 100;
  }
  TGraphErrors *g_fNP_s = new TGraphErrors(n, xv, yv, xe, zero);
  
  // initialize tgraphs for parameters
  TGraphErrors **g_par = new TGraphErrors*[n_p];

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);

  for(int i = 0; i < n_p; i++) {

    g_par[i] = new TGraphErrors(pt_bins, pt_avg, par[i], pt_err, epar[i]);
  }
  
  for(int i = 0; i < n_p; i++) {

    if(i < 2) c->SetLogy();
    else c->SetLogy(0);

    TH1F *fp = c->DrawFrame(pt_min[0]-5, parmin[i], pt_max[pt_bins-1], parmax[i]);
    fp->SetXTitle("p_{T} (GeV)");
    fp->SetYTitle(Form("%s%s", partit[i].c_str(), par_unit[i].c_str()));
    fp->GetYaxis()->SetTitleOffset(1.5);
    fp->GetYaxis()->SetLabelOffset(0.01);
    //fp->SetTitle(Form("2017 %s", partit[i].c_str()));

    // drawing the 1d fit results
    g_par[i]->SetMarkerStyle(20);
    g_par[i]->SetMarkerSize(.75);
    g_par[i]->SetMarkerColor(kBlack);
    g_par[i]->SetLineColor(kBlack);
    g_par[i]->Draw("psame");

    // drawing the 2d fit results
    // constant pars
    if( i == 2 || i == 3) {
      g_par_s[i]->SetLineColor(kRed);
      g_par_s[i]->SetFillColorAlpha(kRed, 0.5);
      g_par_s[i]->Draw("ce3");
    }

    else {
      g_par_s[i]->SetMarkerStyle(20);
      g_par_s[i]->SetMarkerSize(.75);
      g_par_s[i]->SetMarkerColor(kRed);
      g_par_s[i]->SetLineColor(kRed);
      g_par_s[i]->Draw("psame");
    }
    
    if(i == 4) {
      g_par[i+1]->SetMarkerStyle(20);
      g_par[i+1]->SetMarkerSize(.75);
      g_par[i+1]->SetMarkerColor(kBlack);
      g_par[i+1]->SetLineColor(kBlack);
      g_par[i+1]->Draw("psame");

      g_par_s[i+1]->SetMarkerStyle(20);
      g_par_s[i+1]->SetMarkerSize(.75);
      g_par_s[i+1]->SetMarkerColor(kRed);
      g_par_s[i+1]->SetLineColor(kRed);
      g_par_s[i+1]->Draw("psame");

      TLatex ls;
      ls.DrawLatex(60, 11, "#sigma_{1}");
      ls.DrawLatex(60, 21, "#sigma_{2}");

    }

    if(i != 5)
      c->SaveAs(Form("plots/lifetime2d/parB_%s.pdf", parlab[i].c_str()));
    c->Clear();
  }

  TGraphErrors *g_frac = new TGraphErrors(pt_bins, pt_avg, fNP, pt_err, zero);

  TH1F *ffr = c->DrawFrame(pt_min[0]-5, 0, pt_max[pt_bins-1]+5, 50);
  ffr->SetXTitle("p_{T} (GeV)");
  ffr->SetYTitle("f_{NP} (%)");
  ffr->GetYaxis()->SetTitleOffset(1.3);
  ffr->GetYaxis()->SetLabelOffset(0.01);
  ffr->SetTitle("2017 f_{NP} vs p_{T}");
  
  g_frac->SetMarkerStyle(20);
  g_frac->SetMarkerSize(.75);
  g_frac->SetMarkerColor(kBlack);
  g_frac->SetLineColor(kBlack);
  g_frac->Draw("psame");

  g_fNP_s->SetMarkerStyle(20);
  g_fNP_s->SetMarkerSize(.75);
  g_fNP_s->SetMarkerColor(kRed);
  g_fNP_s->SetLineColor(kRed);
  g_fNP_s->Draw("psame");

  c->SaveAs(Form("plots/lifetime2d/parB_fNP.pdf"));
  c->Destructor();
}
