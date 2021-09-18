void plotCosPars2d()
{
  // aux arrays
  const int n_inp = 2;
  int n_p = 3;
  string lbl[] = {"LSB", "RSB"};

  string parlab[] = {"N", "l2", "l4"};
  string partit[] = {"N", "#lambda_{2}", "#lambda_{4}"};
  double parmin[] = {6e-4, -2, -2};
  double parmax[] = {6e-3, 2, 2};

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);
  double pt_min, pt_max;
  
  for(int i_inp = 0; i_inp < n_inp; i_inp++) {
    
    // initialize tgraphs for parameters
    TGraphErrors **g_par = new TGraphErrors*[n_p];
    TFile *fin = new TFile(Form("files/%s2d_fitres.root", lbl[i_inp].c_str()));
    for(int i_p = 0; i_p < n_p; i_p++) {
      fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_p]);
    }
    fin->Close();

    pt_min = g_par[0]->GetX()[0]-g_par[0]->GetEX()[0]-5;
    pt_max = g_par[0]->GetX()[g_par[0]->GetN()-1]+g_par[0]->GetEX()[g_par[0]->GetN()-1]+5;
  
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
	g_par[i]->SetLineColor(kBlack);
	g_par[i]->SetFillColorAlpha(kBlack, 0.5);
	g_par[i]->Draw("ce3");
      }

      else {
	g_par[i]->SetMarkerStyle(20);
	g_par[i]->SetMarkerSize(.75);
	g_par[i]->SetMarkerColor(kBlack);
	g_par[i]->SetLineColor(kBlack);
	g_par[i]->Draw("psame");
      }
  
      c->SaveAs(Form("plots/%s2d/par_%s.pdf", lbl[i_inp].c_str(), parlab[i].c_str()));
      c->Clear();
    }
  }


  // plot everything together
  TH1F *fp = c->DrawFrame(pt_min, 0, pt_max, 1);
  fp->SetXTitle("p_{T} (GeV)");
  fp->SetYTitle(Form("#lambda"));
  fp->GetYaxis()->SetTitleOffset(1.5);
  fp->GetYaxis()->SetLabelOffset(0.01);
  fp->SetTitle(Form("2018 #lambda"));

  for(int i_inp = 0; i_inp < n_inp; i_inp++) {
    
    // initialize tgraphs for parameters
    TGraphErrors **g_par = new TGraphErrors*[n_p];
    TFile *fin = new TFile(Form("files/%s2d_fitres.root", lbl[i_inp].c_str()));
    for(int i_p = 0; i_p < n_p; i_p++) {
      fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_p]);
    }
    fin->Close();
    
    
    for(int i = 1; i < n_p; i++) {

      if(i < 1) c->SetLogy();
      else c->SetLogy(0);

      // constant pars
      if( i > 0) {
	cout << i_inp << " " << i << " " << i+i_inp*(n_p-1) << endl;
	g_par[i]->SetLineColor(i+i_inp*(n_p-1));
	g_par[i]->SetFillColorAlpha(i+i_inp*(n_p-1), 0.5);
	g_par[i]->Draw("ce3 same");
      }

      else {
	g_par[i]->SetMarkerStyle(20);
	g_par[i]->SetMarkerSize(.75);
	g_par[i]->SetMarkerColor(kBlack);
	g_par[i]->SetLineColor(kBlack);
	g_par[i]->Draw("psame");
      }
  
    }
  }
  c->SaveAs(Form("plots/fitCosLbd.pdf"));
  c->Clear();

  c->Destructor();

}
