// macro to draw all the fit parameters
void plotMMPars()
{
  // aux arrays
  const int n_p = 7;

  n_m = 2;
  
  string mode[] = {"free", "dep"};

  string parlab[] = {"f", "N", "mu", "sig1", "sig2", "n", "alpha", "chiN"};
  string partit[] = {"f", "N", "#mu", "#sigma_{1}", "#sigma_{2}", "n", "#alpha"};
  string parax[] = {"f", "N", "#mu (GeV)", "#sigma_{1} (GeV)", "#sigma_{2} (GeV)", "n", "#alpha"};
  double parmin[] = {0, 4.7e-3, 3.096, 0.018, 0.03,  0.6, 1.8};
  double parmax[] = {1, 5e-3,   3.098, 0.03,  0.055, 1.8, 2.5};

  // initialize tgraphs for parameters
  TGraphErrors ***g_par = new TGraphErrors**[n_p];
  for(int i_p = 0; i_p < n_p; i_p++) {
    g_par[i_p] = new TGraphErrors**[n_m];
    for(int i_m = 0; i_m < n_m; i_m++) {
      g_par[i_p][i_m] = new TGraphErrors();
    }
  }
  TGraphErrors *g_chi = new TGraphErrors();
  TLine *l_chi = new TLine();
  
  for(int i_m = 0; i_m < n_m; i_m++) {
    TFile *fin = new TFile(Form("files/mfit_MC_%s.root", mode[i_m].c_str()));

    for(int i_p = 0; i_p < n_p; i_p++) {
      fin->GetObject(Form("%s_%s", parlab[i_p].c_str(), mode[i_m].c_str()), g_par[i_p][i_m]);
      
      if(i_m == 0)
	fin->GetObject(Form("chiN_%s", mode[i_m].c_str()), g_chi);
      else
	fin->GetObject(Form("chiN_%s", mode[i_m].c_str()), l_chi);
    }
    fin->Close();
  }

  double pt_min = g_par[0][0]->GetX()[0]-5;
  double pt_max = g_par[0][0]->GetX()[g_par[0][0]->GetN()-1]+5;
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  for(int i_p = 0; i_p < n_p; i_p++) {

    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("2018 %s", partit[i_p].c_str()));
      
    g_par[i_p][0]->SetLineColor(kBlack);
    g_par[i_p][0]->SetMarkerColor(kBlack);
    g_par[i_p][0]->Draw("p");
    if(i_p == 1 || i_p == 5 || i_p == 6) {
      g_par[i_p][1]->SetLineColor(kBlue);
      g_par[i_p][1]->SetMarkerColor(kBlue);
      if(i_p ==  6 ){
	TF1 *fc = new TF1("fc", "[0]*x+[1]", pt_min, pt_max);
	fc->SetLineColor(kViolet);
	fc->SetLineStyle(kDashed);
	//g_par[i_p][1]->Fit(fc);
      }
      g_par[i_p][1]->Draw("p");

    }
    else {
      g_par[i_p][1]->SetLineColor(kBlue);
      g_par[i_p][1]->SetFillColorAlpha(kBlue, 0.5);
      g_par[i_p][1]->Draw("ce3");
    }
    
    TLine *l1 = new TLine(46, parmin[i_p], 46, parmax[i_p]);
    l1->SetLineColor(kBlack);
    l1->SetLineStyle(kDashed);
    l1->Draw();
    TLine *l2 = new TLine(66, parmin[i_p], 66, parmax[i_p]);
    l2->SetLineColor(kBlack);
    l2->SetLineStyle(kDashed);
    l2->Draw();
    
    c->SaveAs(Form("plots/MCMass/%s.pdf", parlab[i_p].c_str()));
    c->Clear();
  }

  TH1F *fl = c->DrawFrame(pt_min, 0, pt_max, 10);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#chi^{2}/ndf");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle(Form("2018 #chi^{2}/ndf"));
  
  g_chi->SetLineColor(kBlack);
  g_chi->SetMarkerColor(kBlack);
  g_chi->SetMarkerStyle(20);
  g_chi->Draw("p");
  
  l_chi->SetLineColor(kBlue);
  l_chi->Draw("lsame");
  
  c->SaveAs(Form("plots/MCMass/chiN.pdf"));
  c->Clear();
  
  c->Destructor();
}
