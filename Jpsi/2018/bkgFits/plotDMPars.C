// macro to draw all the fit parameters
void plotDMPars(int DO_EXP)
{
  // aux arrays
  const int n_p = 10;

  string parlab[] = {"f", "NS", "mu", "sig1", "sig2", "n", "alpha", "p1", "p2", "fBG", "chiN"};
  string partit[] = {"f", "N_{SR}", "#mu", "#sigma_{1}", "#sigma_{2}", "n", "#alpha", "p_{1}", "p_{2}", "f_{bkg}"};
  string parax[] = {"f (%)", "N_{SR}", "#mu (GeV)", "#sigma_{1} (GeV)", "#sigma_{2} (GeV)", "n", "#alpha", "p_{1} (GeV^{-1})", "p_{2} (GeV)", "f_{bkg} (%)"};
  double parmin[] = {0,    0.008, 3.0934, 0.015, 0.03, 1.0, 2.0, 1e2, 3.6, 0.};
  double parmax[] = {100., 0.01, 3.0944, 0.035, 0.05, 1.4, 2.3, 5e4, 3.8, 10.};

  if(DO_EXP == 1) {
    parmin[7] = 0.;
    parmax[7] = 1.5;
    parmin[8] = 0.5;
    parmax[8] = 0.7;
    parlab[7] = "NB";
    parlab[8] = "lambda";
    partit[7] = "N_{BG}";
    partit[8] = "#lambda";
    parax[7] = "N_{BG}";
    parax[8] = "#lambda (GeV)";
  }
  
  // initialize tgraphs for parameters
  TGraphErrors **g_par = new TGraphErrors*[n_p];
  for(int i_p = 0; i_p < n_p; i_p++) {
    g_par[i_p] = new TGraphErrors();
  }
  TGraphErrors *g_chi = new TGraphErrors();
  TLine *l_chi = new TLine();
  
  TFile *fin = new TFile(Form("files/mfit_data.root"));

  for(int i_p = 0; i_p < n_p; i_p++) {
    fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_p]);
  }

  fin->GetObject(Form("fit_chiN"), l_chi);
      
  fin->Close();

  double pt_min = l_chi->GetX1()-5;
  double pt_max = l_chi->GetX2()+5;
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  for(int i_p = 0; i_p < n_p; i_p++) {
    
    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("2018 %s", partit[i_p].c_str()));
      
    if(i_p == 1 || i_p >= 5) {
      g_par[i_p]->SetLineColor(kBlue);
      g_par[i_p]->SetMarkerColor(kBlue);
      g_par[i_p]->Draw("p");

      }
    else {
      g_par[i_p]->SetLineColor(kBlue);
      g_par[i_p]->SetFillColorAlpha(kBlue, 0.5);
      g_par[i_p]->Draw("ce3");
    }
      
    c->SaveAs(Form("plots/dataMass/%s.pdf", parlab[i_p].c_str()));
    c->Clear();
  }

  TH1F *fl = c->DrawFrame(pt_min, 0, pt_max, 50);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#chi^{2}/ndf");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle(Form("2018 #chi^{2}/ndf"));

  l_chi->SetLineColor(kBlue);
  l_chi->Draw("lsame");
  
  c->SaveAs(Form("plots/dataMass/chiN.pdf"));
  c->Clear();
  
  c->Destructor();
}
