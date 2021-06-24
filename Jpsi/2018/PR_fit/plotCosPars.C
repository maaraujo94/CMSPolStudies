void plotCosPars()
{
  // aux arrays
  const int n_inp = 3;
  int n_p = 2;
  string lbl[] = {"NP", "LSB", "RSB"};

  string parlab[] = {"N", "l2", "l4"};
  string partit[] = {"N", "#lambda_{NP}", "#lambda_{4}"};
  double parmin[] = {3e-2, -1, -1};
  double parmax[] = {5e0, 1, 1};

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);
  
  for(int i_inp = 0; i_inp < n_inp; i_inp++) {
    
    if(i_inp > 0) {
      n_p = 3; 
      parmin[0] = 5e-4;
      parmax[0] = 1e-1;
      partit[1] = "#lambda_{2}";
    }
    // initialize tgraphs for parameters
    TGraphErrors **g_par = new TGraphErrors*[n_p];
    TFile *fin = new TFile(Form("files/%s2d_fitres.root", lbl[i_inp].c_str()));
    for(int i_p = 0; i_p < n_p; i_p++) {
      fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_p]);
    }
    fin->Close();

    double *xv = g_par[0]->GetX();
    double *xe = g_par[0]->GetEX();
    int n = g_par[0]->GetN();
    // scale all graphs for plotting
    TGraphErrors **g_par_s = new TGraphErrors*[n_p];
    double mult[] = {1., 1., 1.};
    for(int i = 0; i < n_p; i++) {
      double *yv = g_par[i]->GetY();
      double *ye = g_par[i]->GetEY();
      
      for(int j = 0; j < n; j++) {
	yv[j] *= mult[i];
	ye[j] *= mult[i];
      }
      g_par_s[i] = new TGraphErrors(n, xv, yv, xe, ye);
    }
    
    double pt_min = xv[0]-xe[0]-5;
    double pt_max = xv[n-1]+xe[n-1]+5;
  
    for(int i = 0; i < n_p; i++) {

      if(i < 1) c->SetLogy();
      else c->SetLogy(0);

      TH1F *fp = c->DrawFrame(pt_min, parmin[i], pt_max, parmax[i]);
      fp->SetXTitle("p_{T} (GeV)");
      fp->SetYTitle(Form("%s", partit[i].c_str()));
      fp->GetYaxis()->SetTitleOffset(1.5);
      fp->GetYaxis()->SetLabelOffset(0.01);
      fp->SetTitle(Form("2018 %s %s", lbl[i_inp].c_str(), partit[i].c_str()));

      // constant pars
      if( i > 0) {
	g_par_s[i]->SetLineColor(kBlack);
	g_par_s[i]->SetFillColorAlpha(kBlack, 0.5);
	g_par_s[i]->Draw("ce3");
      }

      else {
	g_par_s[i]->SetMarkerStyle(20);
	g_par_s[i]->SetMarkerSize(.75);
	g_par_s[i]->SetMarkerColor(kBlack);
	g_par_s[i]->SetLineColor(kBlack);
	g_par_s[i]->Draw("psame");
      }
  
      c->SaveAs(Form("plots/%s2d/par_%s.pdf", lbl[i_inp].c_str(), parlab[i].c_str()));
      c->Clear();
    }
  }
  c->Destructor();

}
