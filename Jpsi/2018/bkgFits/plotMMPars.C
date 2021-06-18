// macro to draw all the fit parameters
void plotMMPars()
{
  // aux arrays
  const int n_p = 7;

  int n_m = 5;
  
  string mode[] = {"free", "mix", "dep0", "dep1", "dep2"};
  string legT[] = {"n, #alpha free", "n fixed, #alpha free", "n, #alpha fixed"};

  string parlab[] = {"f", "N", "mu", "sig1", "sig2", "n", "alpha", "chiN"};
  string partit[] = {"f", "N", "#mu", "#sigma_{1}", "#sigma_{2}", "n", "#alpha"};
  string parax[] = {"f (%)", "N", "#mu (MeV)", "#sigma_{1} (MeV)", "#sigma_{2} (MeV)", "n", "#alpha"};
  double parmin[] = {0,   8e1, 3096,   0.0, 30., 0.6, 1.8};
  double parmax[] = {100, 2e3, 3098,   50., 55., 1.8, 2.5};
  double parmind[] = {60, 8e1, 3096.4, 20., 36., 0.6, 1.8};
  double parmaxd[] = {80, 2e3, 3097.0, 30., 50., 1.8, 2.5};


  // initialize tgraphs for parameters
  TGraphErrors ***g_par = new TGraphErrors**[n_p];
  for(int i_p = 0; i_p < n_p; i_p++) {
    g_par[i_p] = new TGraphErrors*[n_m];
    for(int i_m = 0; i_m < n_m; i_m++) {
      g_par[i_p][i_m] = new TGraphErrors();
    }
  }
  TGraphErrors *g_chi = new TGraphErrors();
  TLine **l_chi = new TLine*[4];
  for(int i_m = 0; i_m < 4; i_m++) {
    l_chi[i_m] = new TLine();
  }
  
  for(int i_m = 0; i_m < n_m; i_m++) {
    if(i_m == 0) {
      TFile *fin = new TFile(Form("files/mfit_MC_%s.root", mode[i_m].c_str()));
      fin->GetObject(Form("chiN_%s", mode[i_m].c_str()), g_chi);

      for(int i_p = 0; i_p < n_p; i_p++)
	fin->GetObject(Form("%s_%s", parlab[i_p].c_str(), mode[i_m].c_str()), g_par[i_p][i_m]);
      fin->Close();
    }
    else if (i_m == 1) {
      TFile *fin = new TFile(Form("files/mfit_MC_%s.root", mode[i_m].c_str()));
      fin->GetObject(Form("chiN_%s", mode[i_m].c_str()), l_chi[0]);
      
      for(int i_p = 0; i_p < n_p; i_p++)
	fin->GetObject(Form("%s_%s", parlab[i_p].c_str(), mode[i_m].c_str()), g_par[i_p][i_m]);
      fin->Close();
      
    }
    
    else {	
      TFile *fin = new TFile(Form("files/mfit_MC_dep.root"));
      fin->GetObject(Form("chiN_%s", mode[i_m].c_str()), l_chi[i_m-1]);

      for(int i_p = 0; i_p < n_p; i_p++)
	fin->GetObject(Form("%s_%s", parlab[i_p].c_str(), mode[i_m].c_str()), g_par[i_p][i_m]);
      fin->Close();
    }
  }

  // get sigma2/sigma1 ratio
  int nv = g_par[0][0]->GetN();
  TGraphErrors **g_sig = new TGraphErrors*[3];
  double vx[2][nv], vy[2][nv], ex[2][nv], ey[2][nv];
  for(int j = 0; j < 3; j++) {
    for(int i = 0; i < nv; i++) {
      vx[j][i] = g_par[0][0]->GetX()[i];
      ex[j][i] = g_par[0][0]->GetEX()[i];
      
      vy[j][i] = g_par[4][j]->GetY()[i]/g_par[3][j]->GetY()[i];
      
      ey[j][i] = vy[j][i]*sqrt(pow(g_par[3][j]->GetEY()[i]/g_par[3][j]->GetY()[i], 2) + pow(g_par[4][j]->GetEY()[i]/g_par[4][j]->GetY()[i], 2));    
    }
    g_sig[j] = new TGraphErrors(nv, vx[j], vy[j], ex[j], ey[j]);
  }
 
  // do the plotting - free vs dep
  double pt_min = g_par[0][0]->GetX()[0]-g_par[0][0]->GetEX()[0]-5;
  double pt_max = g_par[0][0]->GetX()[nv-1]+g_par[0][0]->GetEX()[nv-1]+5;

  TCanvas *c = new TCanvas("", "", 900, 900);
  for(int i_p = 0; i_p < n_p; i_p++) {

    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("2018 %s", partit[i_p].c_str()));

    if(i_p == 1) c->SetLogy();
    else c->SetLogy(0);

    // free mode always plots points
    g_par[i_p][0]->SetLineColor(kBlack);
    g_par[i_p][0]->SetMarkerColor(kBlack);
    g_par[i_p][0]->SetMarkerStyle(20);
    g_par[i_p][0]->SetMarkerSize(.75);
    g_par[i_p][0]->Draw("p");

    // mixed mode only plots lines for (0,2) = (f,mu)
    if(i_p == 0 || i_p == 2) {
      g_par[i_p][1]->SetLineColor(kRed);
      g_par[i_p][1]->SetFillColorAlpha(kRed, 0.5);
      g_par[i_p][1]->Draw("ce3");
    }
    else {
      g_par[i_p][1]->SetLineColor(kRed);
      g_par[i_p][1]->SetMarkerColor(kRed);
      g_par[i_p][1]->SetMarkerStyle(20);
      g_par[i_p][1]->SetMarkerSize(.75);
      g_par[i_p][1]->Draw("p");
    }

    //dep mode plots points for (1,5,6) = (N,n,alpha)
      if(i_p == 1 || i_p == 5 || i_p == 6) {
	g_par[i_p][2]->SetLineColor(kBlue);
	g_par[i_p][2]->SetMarkerColor(kBlue);
	g_par[i_p][2]->SetMarkerStyle(20);
	g_par[i_p][2]->SetMarkerSize(.75);
	g_par[i_p][2]->Draw("p");
      }
      else {
	g_par[i_p][2]->SetLineColor(kBlue);
	g_par[i_p][2]->SetFillColorAlpha(kBlue, 0.5);
	g_par[i_p][2]->Draw("ce3");
      }

      // if we're plotting par 3, add par 4 (both sigmas in 1)
      if( i_p == 3) {      
      g_par[i_p+1][0]->SetLineColor(kBlack);
      g_par[i_p+1][0]->SetMarkerColor(kBlack);
      g_par[i_p+1][0]->SetMarkerStyle(24);
      g_par[i_p+1][0]->SetMarkerSize(.75);
      g_par[i_p+1][0]->Draw("p");

      g_par[i_p+1][1]->SetLineColor(kRed);
      g_par[i_p+1][1]->SetMarkerColor(kRed);
      g_par[i_p+1][1]->SetMarkerStyle(24);
      g_par[i_p+1][1]->SetMarkerSize(.75);
      g_par[i_p+1][1]->Draw("p");

      g_par[i_p+1][2]->SetLineColor(kGreen);
      g_par[i_p+1][2]->SetFillColorAlpha(kGreen, 0.5);
      g_par[i_p+1][2]->Draw("ce3");

      TLegend *leg = new TLegend(0.75, 0.3, 0.9, 0.45);
      leg->SetTextSize(0.03);
      leg->AddEntry(g_par[i_p][0], "#sigma_{1}", "pl");
      leg->AddEntry(g_par[i_p+1][0], "#sigma_{2}", "pl");
      leg->Draw();
    }
    else if (i_p == 4) continue;
    
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

  g_chi->SetMarkerColor(kBlack);
  g_chi->SetLineColor(kBlack);
  g_chi->Draw("p");

  
  l_chi[0]->SetLineColor(kRed);
  l_chi[0]->Draw("lsame");
  l_chi[1]->SetLineColor(kBlue);
  l_chi[1]->Draw("lsame");
  
  c->SaveAs(Form("plots/MCMass/chiN.pdf"));
  c->Clear();

  TH1F *fr = c->DrawFrame(pt_min, 1.6, pt_max, 2);
  fr->SetXTitle("p_{T} (GeV)");
  fr->SetYTitle("#sigma_{2}/#sigma_{1}");
  fr->GetYaxis()->SetTitleOffset(1.3);
  fr->GetYaxis()->SetLabelOffset(0.01);
  fr->SetTitle(Form("2018 #sigma_{2}/#sigma_{1}"));
  
  g_sig[0]->SetLineColor(kBlack);
  g_sig[0]->SetMarkerColor(kBlack);
  g_sig[0]->Draw("p");

  TF1 *fc = new TF1("fc", "[0]*x+[1]", pt_min, pt_max);
  fc->SetLineColor(kBlack);
  fc->SetLineStyle(kDashed);
  g_sig[0]->Fit(fc);

  g_sig[1]->SetLineColor(kRed);
  g_sig[1]->SetMarkerColor(kRed);
  g_sig[1]->Draw("p");
  
  TF1 *fcm = new TF1("fcm", "[0]*x+[1]", pt_min, pt_max);
  fcm->SetLineColor(kRed);
  fcm->SetLineStyle(kDashed);
  g_sig[1]->Fit(fcm);

  
  g_sig[2]->SetLineColor(kBlue);
  g_sig[2]->SetFillColorAlpha(kBlue, 0.5);
  g_sig[2]->Draw("ce3");
  
  c->SaveAs(Form("plots/MCMass/sigR.pdf"));
  c->Clear();

  // do the plotting - the 3 dep models


  for(int i_p = 0; i_p < n_p; i_p++) {

    TH1F *fl = c->DrawFrame(pt_min, parmind[i_p], pt_max, parmaxd[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("2018 %s", partit[i_p].c_str()));

    if(i_p == 1) c->SetLogy();
    else c->SetLogy(0);
    
    TLegend *leg = new TLegend(0.65, 0.65, 0.9, 0.9);
    leg->SetTextSize(0.03);
    
    for(int i_m = 2; i_m < 5; i_m++) {
      if(i_p == 1 || i_p == 5 || i_p == 6) {
	g_par[i_p][i_m]->SetLineColor(i_m);
	g_par[i_p][i_m]->SetMarkerColor(i_m);
	g_par[i_p][i_m]->SetMarkerStyle(20);
	g_par[i_p][i_m]->SetMarkerSize(.75);
	g_par[i_p][i_m]->Draw("p");
      }
      else {
	g_par[i_p][i_m]->SetLineColor(i_m);
	g_par[i_p][i_m]->SetMarkerColor(i_m);
	g_par[i_p][i_m]->SetFillColorAlpha(i_m, 0.5);
	g_par[i_p][i_m]->Draw("ce3");
      }

      leg->AddEntry(g_par[i_p][i_m], legT[i_m-2].c_str(), "pl");
    }
    
    TLine *l1 = new TLine(46, parmind[i_p], 46, parmaxd[i_p]);
    l1->SetLineColor(kBlack);
    l1->SetLineStyle(kDashed);
    l1->Draw();
    TLine *l2 = new TLine(66, parmind[i_p], 66, parmaxd[i_p]);
    l2->SetLineColor(kBlack);
    l2->SetLineStyle(kDashed);
    l2->Draw();

    leg->Draw();
    
    c->SaveAs(Form("plots/MCMass/%s_dep.pdf", parlab[i_p].c_str()));
    c->Clear();
  }

  TF1 *fcc = new TF1("fcc", "[0]", pt_min, pt_max);
  g_par[5][2]->Fit(fcc);

  g_par[6][3]->Fit(fcc);

  
  TH1F *fld = c->DrawFrame(pt_min, 3, pt_max, 7);
  fld->SetXTitle("p_{T} (GeV)");
  fld->SetYTitle("#chi^{2}/ndf");
  fld->GetYaxis()->SetTitleOffset(1.3);
  fld->GetYaxis()->SetLabelOffset(0.01);
  fld->SetTitle(Form("2018 #chi^{2}/ndf"));

  TLegend *leg = new TLegend(0.65, 0.65, 0.9, 0.9);
  leg->SetTextSize(0.03);
  
  for(int i_m = 1; i_m < 4; i_m++) {
    l_chi[i_m]->SetLineColor(i_m+1);
    l_chi[i_m]->Draw("lsame");
    leg->AddEntry(l_chi[i_m], legT[i_m-1].c_str(), "l");
  }

  leg->Draw();
  
  c->SaveAs(Form("plots/MCMass/chiN_dep.pdf"));
  c->Clear();
  
  c->Destructor();
}
