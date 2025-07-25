// macro to draw all the fit parameters

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

void plotMMG()
{
  // aux arrays
  const int n_p = 9;
  int n_m = 2;
  string fn[] = {"4", "G"};
  
  string parlab[] = {"f", "N", "mu", "sig1", "sig2", "n", "alpha", "fG" , "sigG"};
  string partit[] = {"f", "N", "#mu", "#sigma", "#sigma_{2}", "n", "#alpha", "f_{G}", "#sigma_{G}"};
  string parax[] = {"f (%)", "N per 1 GeV", "#mu (MeV)", "#sigma (MeV)", "#sigma_{2} (MeV)", "n", "#alpha", "f_{G} (%)", "#sigma_{G} (MeV)"};
  double parmin[] = {0,   3,   3650, 0.0, 30., 0.6, 1.8, 0,   0};
  double parmax[] = {100, 9e3, 3750, 110, 55., 1.8, 2.5, 100, 100};

  // initialize tgraphs for parameters
  TGraphErrors ***g_par = new TGraphErrors**[n_m];
  for(int i_m = 0; i_m < n_m; i_m++) {
    g_par[i_m] = new TGraphErrors*[n_p];
    for(int i_p = 0; i_p < n_p; i_p++) {
      g_par[i_m][i_p] = new TGraphErrors();
    }
  }

  // get params from fit result files
  for(int i_m = 0; i_m < n_m; i_m++) {
    TFile *fin = new TFile(Form("files/MCfit_%s.root", fn[i_m].c_str()));
    int n_pv = i_m == 0 ? 7 : 9;
    for(int i_p = 0; i_p < n_pv; i_p++) {
      fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_m][i_p]);
    }
    fin->Close();
  }
  
  int nv = g_par[0][0]->GetN();
  double *vx = g_par[0][0]->GetX();
  double *ex = g_par[0][0]->GetEX();
 
  // do the plotting - compare models
  double pt_min = vx[0]-ex[0]-5;
  double pt_max = vx[nv-1]+ex[nv-1]+5;

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.03);

  for(int i_p = 0; i_p < 7; i_p++) {

    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.8);
    if(i_p == 3) fl->GetYaxis()->SetTitleOffset(1.7);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("MC %s vs p_{T}", partit[i_p].c_str()));

    if (i_p == 4) continue;
    
    if(i_p == 1) c->SetLogy();
    else c->SetLogy(0);

    TLegend *leg = new TLegend(0.77, 0.7, 0.97, 0.9);
    leg->SetTextSize(0.03);

    // drawing base model
    // free pars
    if(i_p == 1 || i_p == 6) {
      g_par[0][i_p]->SetMarkerStyle(20);
      g_par[0][i_p]->SetMarkerSize(.75);
      g_par[0][i_p]->SetLineColor(kBlack);
      g_par[0][i_p]->SetMarkerColor(kBlack);
      g_par[0][i_p]->Draw("p");

      g_par[1][i_p]->SetMarkerStyle(20);
      g_par[1][i_p]->SetMarkerSize(.75);
      g_par[1][i_p]->SetLineColor(kBlue);
      g_par[1][i_p]->SetMarkerColor(kBlue);
      g_par[1][i_p]->Draw("p");	
    }
    // linear or constant pars 
    else if(i_p != 3) {
      g_par[0][i_p]->SetMarkerStyle(20);
      g_par[0][i_p]->SetMarkerSize(.75);
      g_par[0][i_p]->SetMarkerColor(kBlack);
      g_par[0][i_p]->SetLineColor(kBlack);
      g_par[0][i_p]->SetFillColorAlpha(kBlack, 0.5);
      g_par[0][i_p]->Draw("pce3");

      g_par[1][i_p]->SetMarkerStyle(20);
      g_par[1][i_p]->SetMarkerSize(.75);
      g_par[1][i_p]->SetMarkerColor(kBlue);
      g_par[1][i_p]->SetLineColor(kBlue);
      g_par[1][i_p]->SetFillColorAlpha(kBlue, 0.5);
      g_par[1][i_p]->Draw("le3");
    }

    // SPECIAL CASES - COMBINED PLOTS
    // if we're plotting f, add fG
    if( i_p == 0) {      
      g_par[1][7]->SetMarkerStyle(20);
      g_par[1][7]->SetMarkerSize(.75);
      g_par[1][7]->SetLineColor(kBlue);
      g_par[1][7]->SetMarkerColor(kBlue);
      g_par[1][7]->Draw("p");

      // fit f_G with constant function
      TF1 *f_fg = new TF1("f_fg", "[0]+[1]*x", 20, 60);
      f_fg->FixParameter(1,0);
      g_par[1][7]->Fit(f_fg, "R0");
      
      leg->AddEntry(g_par[1][0], "f_{CB1}", "l");
      leg->AddEntry(g_par[1][7], "f_{G}", "pl");
      leg->Draw();
    }

    // if we're plotting par sig1, add sig2 and sigG
    else if( i_p == 3) {
      g_par[0][i_p]->SetMarkerStyle(20);
      g_par[0][i_p]->SetMarkerSize(.75);
      g_par[0][i_p]->SetMarkerColor(kBlack);
      g_par[0][i_p]->SetLineColor(kBlack);
      g_par[0][i_p]->SetFillColorAlpha(kBlack, 0.5);
      g_par[0][i_p]->Draw("pce3");

      g_par[1][i_p]->SetMarkerStyle(20);
      g_par[1][i_p]->SetMarkerSize(.75);
      g_par[1][i_p]->SetMarkerColor(kBlue);
      g_par[1][i_p]->SetLineColor(kBlue);
      g_par[1][i_p]->SetFillColorAlpha(kBlue, 0.5);
      g_par[1][i_p]->Draw("pce3");

      g_par[0][4]->SetMarkerStyle(22);
      g_par[0][4]->SetMarkerSize(.75);
      g_par[0][4]->SetMarkerColor(kBlack);
      g_par[0][4]->SetLineColor(kBlack);
      g_par[0][4]->SetFillColorAlpha(kBlack, 0.5);
      g_par[0][4]->Draw("pce3");

      g_par[1][4]->SetMarkerStyle(22);
      g_par[1][4]->SetMarkerColor(kBlue);
      g_par[1][4]->SetLineColor(kBlue);
      g_par[1][4]->SetFillColorAlpha(kBlue, 0.5);
      g_par[1][4]->Draw("pce3");

      g_par[1][8]->SetMarkerStyle(29);
      g_par[1][8]->SetMarkerColor(kBlue);
      g_par[1][8]->SetLineColor(kBlue);
      g_par[1][8]->Draw("p");

      // fit sigma_G with linear function
      TF1 *f_sg = new TF1("f_sg", "[0]+[1]*x", 20, 60);
      f_sg->SetParameters(60, 0.1);
      g_par[1][8]->Fit(f_sg, "R0");

      leg->AddEntry(g_par[1][3], "#sigma_{CB_{1}}", "p");
      leg->AddEntry(g_par[1][4], "#sigma_{CB_{2}}", "p");
      leg->AddEntry(g_par[1][8], "#sigma_{G}", "p");
      leg->Draw();

    }
    
    int isLog = 0;
    if(i_p == 1 ) isLog = 1;
    
    TLatex lc;
    lc.SetTextSize(0.03);
    lc.DrawLatex(70, getPos(parmin[i_p], parmax[i_p], 0.9, isLog), "no G");
    lc.SetTextColor(kBlue);
    lc.DrawLatex(70, getPos(parmin[i_p], parmax[i_p], 0.83, isLog), "with G");
    
    c->SaveAs(Form("plots/MCMass/parG_%s.pdf", parlab[i_p].c_str()));
    c->Clear();
  }
  
  c->Destructor();
}
