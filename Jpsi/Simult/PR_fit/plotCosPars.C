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
  const int n_inp = 2; // LSB, RSB
  int n_p = 4; // N, l2, l4, chiP
  int n_f = 2; // free, l2 fixed
  string lbl[] = {"LSB", "RSB"};
  string lbl_f[] = {"free", "fixed"};
  
  string parlab[] = {"N", "l2", "l4", "chiP"};
  string partit[] = {"N", "#lambda_{2}", "#lambda_{4}", "P(#chi^{2})"};
  double parmin[] = {6e-4, -3, -3, 0};
  double parmax[] = {6e-3, 4, 4, 1};

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);

  // get the parameter TGraphs
  TGraphErrors ****g_par = new TGraphErrors***[n_inp];
  for(int i_inp = 0; i_inp < n_inp; i_inp++) {
    TFile *fin = new TFile(Form("files/%s_fitres.root", lbl[i_inp].c_str()));	
    // initialize tgraphs for parameters
    g_par[i_inp] = new TGraphErrors**[n_f];
    for(int i_f = 0; i_f < n_f; i_f++) {
      g_par[i_inp][i_f] = new TGraphErrors*[n_p-1];
      for(int i_p = 0; i_p < n_p-1; i_p++) {
	fin->GetObject(Form("fit_%d_%s", i_f, parlab[i_p].c_str()), g_par[i_inp][i_f][i_p]);
      }
    }
    fin->Close();
  }
  
  double pt_min = g_par[0][0][0]->GetX()[0]-g_par[0][0][0]->GetEX()[0]-5;
  double pt_max = g_par[0][0][0]->GetX()[g_par[0][0][0]->GetN()-1]+g_par[0][0][0]->GetEX()[g_par[0][0][0]->GetN()-1]+5;
  
  for(int i = 0; i < n_p-1; i++) {
    if(i < 1) c->SetLogy();
    else c->SetLogy(0);
    for(int i_f = 0; i_f < n_f; i_f++) {
      TH1F *fp = c->DrawFrame(pt_min, parmin[i], pt_max, parmax[i]);
      fp->SetXTitle("p_{T} (GeV)");
      fp->SetYTitle(Form("%s", partit[i].c_str()));
      fp->GetYaxis()->SetTitleOffset(1.5);
      fp->GetYaxis()->SetLabelOffset(0.01);
      fp->SetTitle(Form("%s (%s #lambda_{2})", partit[i].c_str(), lbl_f[i_f].c_str()));

      TLegend *leg = new TLegend(0.75, 0.7, 0.9, 0.9);
      leg->SetTextSize(0.03);
      
      for(int i_inp = 0; i_inp < n_inp; i_inp++) {
	g_par[i_inp][i_f][i]->SetMarkerStyle(20);
	g_par[i_inp][i_f][i]->SetMarkerSize(.75);
	g_par[i_inp][i_f][i]->SetMarkerColor(col_s(i_inp));
	g_par[i_inp][i_f][i]->SetLineColor(col_s(i_inp));
	g_par[i_inp][i_f][i]->Draw("psame");

	leg->AddEntry(g_par[i_inp][i_f][i], lbl[i_inp].c_str(), "pl");
      }

      leg->Draw();
      
      c->SaveAs(Form("plots/SB/par_%s_%s.pdf", parlab[i].c_str(), lbl_f[i_f].c_str()));
      c->Clear();
    }
  }
  c->Destructor();
  
}
