// macro to draw all the fit parameters
void plotnalpha()
{
  // aux arrays
  int n_m = 4;
  int n_p = 2;
  
  string parlab[] = {"n", "alpha"};
  string partit[] = {"n", "#alpha"};
  string parax[] = {"n", "#alpha"};
  double parmin[] = {0.6, 1.9};
  double parmax[] = {1.8, 2.5};

  // initialize tgraphs for parameters
  TGraphErrors ***g_par = new TGraphErrors**[n_p];
  for(int i_p = 0; i_p < n_p; i_p++) {
    g_par[i_p] = new TGraphErrors*[n_m];
    for(int i_m = 0; i_m < n_m; i_m++) {
      g_par[i_p][i_m] = new TGraphErrors();
    }
  }

  // get params from fit result files
  for(int i_m = 0; i_m < n_m; i_m++) {
    TFile *fin = new TFile(Form("files/MCfit_%d.root", i_m));    
    for(int i_p = 0; i_p < n_p; i_p++)
      fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_p][i_m]);
    fin->Close();
  }
  
  int nv = g_par[0][0]->GetN();
  double *vx = g_par[0][0]->GetX();
  double *ex = g_par[0][0]->GetEX();

  // make a TGraph n vs alpha
  TGraphErrors **g_both = new TGraphErrors*[n_m];
  for(int i_m = 0; i_m < n_m; i_m++) {
    double *vn = g_par[0][i_m]->GetY();
    double *en = g_par[0][i_m]->GetEY();
    double *va = g_par[1][i_m]->GetY();
    double *ea = g_par[1][i_m]->GetEY();

    g_both[i_m] = new TGraphErrors(nv, vn, va, en, ea);
  }
  
  // do the plotting
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.12);
  int col[] = {kBlack, kBlue, kViolet, kRed};
  string cond[] = {"free", "#mu", "#mu, f", "#mu, f, #sigma_{1,2}"};

  TH1F *fl = c->DrawFrame(parmin[0], parmin[1], parmax[0], parmax[1]);
  fl->SetXTitle(parax[0].c_str());
  fl->SetYTitle(parax[1].c_str());
  fl->GetYaxis()->SetTitleOffset(1.8);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle(Form("2017 %s vs %s", partit[0].c_str(), partit[1].c_str()));

  TLegend *leg = new TLegend(0.65, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  for(int i_m = 0; i_m < n_m; i_m++) {

    // free mode always plots points
    g_both[i_m]->SetLineColor(col[i_m]);
    g_both[i_m]->SetMarkerColor(col[i_m]);
    g_both[i_m]->SetMarkerStyle(20);
    g_both[i_m]->SetMarkerSize(.75);
    g_both[i_m]->Draw("p");

    leg->AddEntry(g_both[i_m], cond[i_m].c_str(), "pl");

    TF1 *fc = new TF1("fc", "[0]*x+[1]", 0.6, 1.8);
    fc->SetParameters(-1, 5);
    fc->SetLineColor(col[i_m]);
    //g_both[i_m]->Fit(fc);
  }
  leg->Draw();
  
  c->SaveAs(Form("plots/MCMass/na_comp.pdf"));
  c->Clear();
  

}
