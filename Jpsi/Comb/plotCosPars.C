int col_s(int inp)
{
  if(inp == 0) return kBlack;
  if(inp == 1) return kBlue;
  if(inp == 2) return kRed;
  else return kWhite;
}

void plotCosPars()
{
  // aux arrays
  const int n_inp = 2;
  int n_p = 2;
  int n_f = 2;
  string lbl[] = {"LSB", "RSB"};
  string st[] = {"free", "fixed"};

  string parlab[] = {"l2", "l4"};
  string partit[] = {"#lambda_{2}", "#lambda_{4}"};
  double parmin[] = {-3, -3};
  double parmax[] = {4, 4};

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);

  // open files for each plot to make less tgraphs
  // each plot has one lambda and one sb
  // each plot has both years and both fit models
  for(int i_inp = 0; i_inp < n_inp; i_inp++) {
    for(int i_p = 0; i_p < n_p; i_p++) {
      TFile *fin17 = new TFile(Form("../2017/PR_fit/files/%s_fitres.root", lbl[i_inp].c_str()));
      TGraphErrors **g_par17_b = new TGraphErrors*[n_f];
      for(int i_f = 0; i_f < n_f; i_f++) 
	fin17->GetObject(Form("fit_%d_%s", i_f, parlab[i_p].c_str()), g_par17_b[i_f]);
      fin17->Close();

      TFile *fin18 = new TFile(Form("../2018/PR_fit/files/%s_fitres.root", lbl[i_inp].c_str()));
      TGraphErrors **g_par18_b = new TGraphErrors*[n_f];
      for(int i_f = 0; i_f < n_f; i_f++) 
	fin18->GetObject(Form("fit_%d_%s", i_f, parlab[i_p].c_str()), g_par18_b[i_f]);
      fin18->Close();

      // make lambda plots slightly moved in x
      // i_f = 1 moves over 0.5; 17->18 moves over 1.0
      TGraphAsymmErrors **g_par17 = new TGraphAsymmErrors*[n_f];
      TGraphAsymmErrors **g_par18 = new TGraphAsymmErrors*[n_f];
      for(int i_f = 0; i_f < n_f; i_f++) {
	int nv = g_par17_b[i_f]->GetN();
	double xv[nv], exl[nv], exh[nv], yv[nv], ey[nv];
	for(int i_v = 0; i_v < nv; i_v++) {
	  xv[i_v] = g_par17_b[i_f]->GetX()[i_v]-1.5+0.5*i_f;
	  exl[i_v] = g_par17_b[i_f]->GetEX()[i_v]-1.5+0.5*i_f;
	  exh[i_v] = g_par17_b[i_f]->GetEX()[i_v]+1.5-0.5*i_f;
	  yv[i_v] = g_par17_b[i_f]->GetY()[i_v];
	  ey[i_v] = g_par17_b[i_f]->GetEY()[i_v];
	}
	g_par17[i_f] = new TGraphAsymmErrors(nv, xv, yv, exl, exh, ey, ey);

	for(int i_v = 0; i_v < nv; i_v++) {
	  xv[i_v] = g_par18_b[i_f]->GetX()[i_v]+1. + 0.5*i_f;
	  exl[i_v] = g_par18_b[i_f]->GetEX()[i_v]+1. + 0.5*i_f;
	  exh[i_v] = g_par18_b[i_f]->GetEX()[i_v]-1. - 0.5*i_f;
	  yv[i_v] = g_par18_b[i_f]->GetY()[i_v];
	  ey[i_v] = g_par18_b[i_f]->GetEY()[i_v];
	}
	g_par18[i_f] = new TGraphAsymmErrors(nv, xv, yv, exl, exh, ey, ey);
      }

      double pt_min = g_par17_b[0]->GetX()[0]-g_par17_b[0]->GetEX()[0]-5;
      double pt_max = g_par17_b[0]->GetX()[g_par17_b[0]->GetN()-1]+g_par17_b[0]->GetEX()[g_par17_b[0]->GetN()-1]+5;
  
      c->SetLogy(0);


      TH1F *fp = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
      fp->SetXTitle("p_{T} (GeV)");
      fp->SetYTitle(Form("%s", partit[i_p].c_str()));
      fp->GetYaxis()->SetTitleOffset(1.5);
      fp->GetYaxis()->SetLabelOffset(0.01);
      fp->SetTitle(Form("%s %s", lbl[i_inp].c_str(), partit[i_p].c_str()));

      
      for(int i_f = 0; i_f < n_f; i_f++) {
	g_par17[i_f]->SetMarkerStyle(25-i_f);
	//g_par17[i_f]->SetMarkerSize(.75);
	g_par17[i_f]->SetMarkerColor(col_s(i_f));
	g_par17[i_f]->SetLineColor(col_s(i_f));
	g_par17[i_f]->Draw("psame");
	
	g_par18[i_f]->SetMarkerStyle(21-i_f);
	//g_par18[i_f]->SetMarkerSize(.75);
	g_par18[i_f]->SetMarkerColor(col_s(i_f));
	g_par18[i_f]->SetLineColor(col_s(i_f));
	g_par18[i_f]->Draw("psame");
     } 

      TLegend *leg = new TLegend(0.65, 0.7, 0.9, 0.9);
      leg->SetTextSize(0.03);
      for(int i_f = 0; i_f < n_f; i_f++) {
	leg->AddEntry(g_par17[i_f], Form("2017 %s #lambda_{2}", st[i_f].c_str()), "pl");
	leg->AddEntry(g_par18[i_f], Form("2018 %s #lambda_{2}", st[i_f].c_str()), "pl");
      }
      if(i_p == 1)
	leg->Draw();
      
      c->SaveAs(Form("par_%s%s.pdf", lbl[i_inp].c_str(), parlab[i_p].c_str()));
      c->Clear();
    }
    
  }
  c->Destructor();
  
}
