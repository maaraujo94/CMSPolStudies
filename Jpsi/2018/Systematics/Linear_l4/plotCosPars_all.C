int col_s(int inp)
{
  if(inp == 0) return kBlack;
  if(inp == 1) return kBlue;
  if(inp == 2) return kRed;
  if(inp == 3) return kGreen;
  else return kWhite;
}

void plotCosPars_all()
{
  // aux arrays
  const int n_inp = 2;
  int n_p = 3;
  int n_fit = 3;
  string lbl[] = {"LSB", "RSB"};
  string loc[] = {"", "../../PR_fit/", ""};
  string filename[] = {"", "2d", "lin"};

  string parlab[] = {"N", "l2", "l4", "chiP"};
  string partit[] = {"N", "#lambda_{2}", "#lambda_{4}", "P(#chi^{2})"};
  double parmin[] = {6e-4, -2, -2, 0};
  double parmax[] = {6e-3, 2, 2, 1};

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);

  // initialize tgraphs for parameters
  TGraphErrors ****g_par = new TGraphErrors***[n_fit];
  TGraph *g_chi = new TGraph();
  TLine **l_chi =  new TLine*[n_fit-1];

  // cycle to fill TGraphs
  for(int i_f = 0; i_f < n_fit; i_f++) { // cycle through fit format
    g_par[i_f] = new TGraphErrors**[n_inp];
    for(int i_inp = 0; i_inp < n_inp; i_inp++) { // cycle through SB
      // open file
      TFile *fin = new TFile(Form("%sfiles/%s%s_fitres.root", loc[i_f].c_str(), lbl[i_inp].c_str(), filename[i_f].c_str()));	

      // get fit pars
      g_par[i_f][i_inp] = new TGraphErrors*[n_p];
      for(int i_p = 0; i_p < n_p; i_p++) {
	fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_f][i_inp][i_p]);	  
      }

      // get fit chi^2
      if(i_f == 0)
	fin->GetObject("fit_chiP", g_chi);
      else
	fin->GetObject("fit_chiP", l_chi[i_f-1]);

      fin->Close();
    }
  }
  
  double pt_min = g_par[0][0][0]->GetX()[0]-g_par[0][0][0]->GetEX()[0]-5;
  double pt_max = g_par[0][0][0]->GetX()[g_par[0][0][0]->GetN()-1]+g_par[0][0][0]->GetEX()[g_par[0][0][0]->GetN()-1]+5;
  
  for(int i_inp = 0; i_inp < n_inp; i_inp++) {
    for(int i = 0; i < n_p; i++) {

      if(i < 1) c->SetLogy();
      else c->SetLogy(0);

      TH1F *fp = c->DrawFrame(pt_min, parmin[i], pt_max, parmax[i]);
      fp->SetXTitle("p_{T} (GeV)");
      fp->SetYTitle(Form("%s", partit[i].c_str()));
      fp->GetYaxis()->SetTitleOffset(1.5);
      fp->GetYaxis()->SetLabelOffset(0.01);
      fp->SetTitle(Form("2018 %s %s", lbl[i_inp].c_str(), partit[i].c_str()));

      for(int i_f = 0; i_f < n_fit; i_f++) {
	// plot free parameters (all for fit 1, N for fit 2, 3)
	if(i_f == 0 || i == 0) {
	  g_par[i_f][i_inp][i]->SetMarkerStyle(20);
	  g_par[i_f][i_inp][i]->SetMarkerSize(.75);
	  g_par[i_f][i_inp][i]->SetMarkerColor(col_s(i_f));
	  g_par[i_f][i_inp][i]->SetLineColor(col_s(i_f));
	  g_par[i_f][i_inp][i]->Draw("psame");
	}

	else {
	g_par[i_f][i_inp][i]->SetLineColor(col_s(i_f));
	g_par[i_f][i_inp][i]->SetFillColorAlpha(col_s(i_f), 0.5);
	g_par[i_f][i_inp][i]->Draw("ce3");
	}
      
      }
      
      c->SaveAs(Form("plots/%s/parA_%s.pdf", lbl[i_inp].c_str(), parlab[i].c_str()));
      c->Clear();
    }
    
    // draw chi^2 prob
    TH1F *fp = c->DrawFrame(pt_min, parmin[n_p], pt_max, parmax[n_p]);
    fp->SetXTitle("p_{T} (GeV)");
    fp->SetYTitle(Form("%s", partit[n_p].c_str()));
    fp->GetYaxis()->SetTitleOffset(1.5);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("2018 %s %s", lbl[i_inp].c_str(), partit[n_p].c_str()));
    
    g_chi->SetMarkerStyle(20);
    g_chi->SetMarkerSize(1.5);
    g_chi->SetMarkerColor(col_s(0));
    g_chi->SetLineColor(col_s(0));
    g_chi->Draw("psame");

    for(int i = 0; i < n_fit-1; i++) {
      l_chi[i]->SetLineColor(col_s(i+1));
      l_chi[i]->Draw("lsame");
    }
    
    c->SaveAs(Form("plots/%s/parA_%s.pdf", lbl[i_inp].c_str(), parlab[n_p].c_str()));
    c->Clear();
    
  }
  c->Destructor();
  
}
