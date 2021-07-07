int col_s(int inp)
{
  if(inp == 0) return kBlack;
  if(inp == 1) return kBlue;
  if(inp == 2) return kRed;
  if(inp == 3) return kGreen;
  else return kWhite;
}

void plotCosPars_both()
{
  // aux arrays
  const int n_inp = 2;
  int n_p = 4;
  int n_f = 3;
  string lbl[] = {"LSB", "RSB"};

  string parlab[] = {"N", "l2", "l4", "chiP"};
  string partit[] = {"N", "#lambda_{2}", "#lambda_{4}", "P(#chi^{2})"};
  double parmin[] = {6e-4, -2, -2, 0};
  double parmax[] = {6e-3, 2, 2, 1};

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);
  
  for(int i_inp = 0; i_inp < n_inp; i_inp++) {
    TFile *fin = new TFile(Form("files/%s_fitres.root", lbl[i_inp].c_str()));	
    // initialize tgraphs for parameters
    TGraphErrors ***g_par = new TGraphErrors**[n_f+1];
    TGraph **g_chi = new TGraph*[n_f];
     for(int i_f = 0; i_f < n_f; i_f++) {
      g_par[i_f] = new TGraphErrors*[n_p];
      for(int i_p = 0; i_p < n_p; i_p++) {
	fin->GetObject(Form("fit_%d_%s", i_f, parlab[i_p].c_str()), g_par[i_f][i_p]);
      }
      fin->GetObject(Form("fit_%d_%s", i_f, parlab[n_p-1].c_str()), g_chi[i_f]);
    }
    fin->Close();

    TFile *fin2 = new TFile(Form("files/%s2d_fitres.root", lbl[i_inp].c_str()));
    g_par[n_f] = new TGraphErrors*[n_p];
    for(int i_p = 0; i_p < n_p; i_p++) {
      fin2->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[n_f][i_p]);
    }
    TLine *l_chi =  new TLine();
    fin2->GetObject(Form("fit_%s", parlab[n_p-1].c_str()), l_chi);
    fin2->Close();

    double pt_min = g_par[0][0]->GetX()[0]-g_par[0][0]->GetEX()[0]-5;
    double pt_max = g_par[0][0]->GetX()[g_par[0][0]->GetN()-1]+g_par[0][0]->GetEX()[g_par[0][0]->GetN()-1]+5;
  
    for(int i = 0; i < n_p-1; i++) {

      if(i < 1) c->SetLogy();
      else c->SetLogy(0);

      TH1F *fp = c->DrawFrame(pt_min, parmin[i], pt_max, parmax[i]);
      fp->SetXTitle("p_{T} (GeV)");
      fp->SetYTitle(Form("%s", partit[i].c_str()));
      fp->GetYaxis()->SetTitleOffset(1.5);
      fp->GetYaxis()->SetLabelOffset(0.01);
      fp->SetTitle(Form("2017 %s %s", lbl[i_inp].c_str(), partit[i].c_str()));

      for(int i_f = 0; i_f < 1; i_f++) {
	g_par[i_f][i]->SetMarkerStyle(20);
	g_par[i_f][i]->SetMarkerSize(.75);
	g_par[i_f][i]->SetMarkerColor(col_s(i_f));
	g_par[i_f][i]->SetLineColor(col_s(i_f));
	g_par[i_f][i]->Draw("psame");
      }

      if(i < 1) {
	g_par[n_f][i]->SetMarkerStyle(20);
	g_par[n_f][i]->SetMarkerSize(.75);
	g_par[n_f][i]->SetMarkerColor(col_s(n_f));
	g_par[n_f][i]->SetLineColor(col_s(n_f));
	g_par[n_f][i]->Draw("psame");
      }

      else {
	g_par[n_f][i]->SetLineColor(col_s(n_f));
	g_par[n_f][i]->SetFillColorAlpha(col_s(n_f), 0.5);
	g_par[n_f][i]->Draw("ce3");
      }
	
      c->SaveAs(Form("plots/%s/parB_%s.pdf", lbl[i_inp].c_str(), parlab[i].c_str()));
      c->Clear();
    }

    // draw chi^2 prob
    TH1F *fp = c->DrawFrame(pt_min, parmin[n_p-1], pt_max, parmax[n_p-1]);
    fp->SetXTitle("p_{T} (GeV)");
    fp->SetYTitle(Form("%s", partit[n_p-1].c_str()));
    fp->GetYaxis()->SetTitleOffset(1.5);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("2017 %s %s", lbl[i_inp].c_str(), partit[n_p-1].c_str()));

    for(int i_f = 0; i_f < 1; i_f++) {
      g_chi[i_f]->SetMarkerStyle(20);
      g_chi[i_f]->SetMarkerSize(1.5);
      g_chi[i_f]->SetMarkerColor(col_s(i_f));
      g_chi[i_f]->SetLineColor(col_s(i_f));
      g_chi[i_f]->Draw("psame");
    }
    l_chi->SetLineColor(col_s(n_f));
    l_chi->Draw("lsame");
    
    c->SaveAs(Form("plots/%s/parB_%s.pdf", lbl[i_inp].c_str(), parlab[n_p-1].c_str()));
    c->Clear();
     
  }
  c->Destructor();
  
}
