// macro to draw all the fit parameters
void plotMMseq()
{
  // aux arrays
  const int n_p = 7;

  int n_m = 6;
  
  string parlab[] = {"f", "N", "mu", "sig1", "sig2", "n", "alpha", "chiN"};
  string partit[] = {"f", "N", "#mu", "#sigma", "#sigma_{2}", "n", "#alpha"};
  string parax[] = {"f (%)", "N per 1 GeV", "#mu (MeV)", "#sigma (MeV)", "#sigma_{2} (MeV)", "n", "#alpha"};
  double parmin[] = {0,   7e1, 3090,   0.0, 30., 0.6, 1.8};
  double parmax[] = {100, 4e3, 3100,   60., 55., 1.8, 2.5};

  // initialize tgraphs for parameters
  TGraphErrors ***g_par = new TGraphErrors**[n_p];
  for(int i_p = 0; i_p < n_p; i_p++) {
    g_par[i_p] = new TGraphErrors*[n_m];
    for(int i_m = 0; i_m < n_m; i_m++) {
      g_par[i_p][i_m] = new TGraphErrors();
    }
  }
  TLine **l_chiN = new TLine*[n_m];
  for(int i_m = 0; i_m < n_m; i_m++) {
    l_chiN[i_m] = new TLine();
  }
  TLine **l_chi= new TLine*[n_m];
  for(int i_m = 0; i_m < n_m; i_m++) {
    l_chi[i_m] = new TLine();
  }
  TLine **l_ndf = new TLine*[n_m];
  for(int i_m = 0; i_m < n_m; i_m++) {
    l_ndf[i_m] = new TLine();
  }

  // get params from fit result files
  for(int i_m = 0; i_m < n_m; i_m++) {
    TFile *fin = new TFile(Form("files/MCfit_%d.root", i_m));
    fin->GetObject("fit_chiN", l_chiN[i_m]);
    fin->GetObject("fit_chi", l_chi[i_m]);
    fin->GetObject("fit_ndf", l_ndf[i_m]);
    
    for(int i_p = 0; i_p < n_p; i_p++)
      fin->GetObject(Form("fit_%s", parlab[i_p].c_str()), g_par[i_p][i_m]);
    fin->Close();
  }
  
  int nv = g_par[0][0]->GetN();
  double *vx = g_par[0][0]->GetX();
  double *ex = g_par[0][0]->GetEX();
 
  // do the plotting - free param fit (MODEL 0)
  double pt_min = vx[0]-ex[0]-5;
  double pt_max = vx[nv-1]+ex[nv-1]+5;

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.12);
  for(int i_p = 0; i_p < n_p; i_p++) {

    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.8);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("Full %s", partit[i_p].c_str()));

    if(i_p == 1) c->SetLogy();
    else c->SetLogy(0);

    // free mode always plots points
    g_par[i_p][0]->SetLineColor(kBlack);
    g_par[i_p][0]->SetMarkerColor(kBlack);
    g_par[i_p][0]->SetMarkerStyle(20);
    g_par[i_p][0]->SetMarkerSize(.75);
    g_par[i_p][0]->Draw("p");

    // if we're plotting par 3, add par 4 (both sigmas in 1)
    if( i_p == 3) {      
      g_par[i_p+1][0]->SetLineColor(kBlue);
      g_par[i_p+1][0]->SetMarkerColor(kBlue);
      g_par[i_p+1][0]->SetMarkerStyle(24);
      g_par[i_p+1][0]->SetMarkerSize(.75);
      g_par[i_p+1][0]->Draw("p");
      
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
    
    c->SaveAs(Form("plots/MCMass/par_%s.pdf", parlab[i_p].c_str()));

    if(i_p == 2) {
      g_par[i_p][1]->SetLineColor(kRed);
      g_par[i_p][1]->SetFillColorAlpha(kRed, 0.5);
      g_par[i_p][1]->Draw("ce3");
 
      c->SaveAs(Form("plots/MCMass/par_%s_next.pdf", parlab[i_p].c_str()));
    }
   
    c->Clear();
  }
  
  // do the plotting - f, sigmas for fixed mu (MODEL 1)
  for(int i_p = 0; i_p < 4; i_p++) {
    if(i_p == 1 || i_p == 2) continue;
    
    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.8);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("Full %s (#mu constant)", partit[i_p].c_str()));

    c->SetLogy(0);

    // free mode always plots points
    g_par[i_p][1]->SetLineColor(kBlack);
    g_par[i_p][1]->SetMarkerColor(kBlack);
    g_par[i_p][1]->SetMarkerStyle(20);
    g_par[i_p][1]->SetMarkerSize(.75);
    g_par[i_p][1]->Draw("p");

    if(i_p == 0) {
      g_par[i_p][2]->SetLineColor(kRed);
      g_par[i_p][2]->SetFillColorAlpha(kRed, 0.5);
      g_par[i_p][2]->Draw("ce3");
    }

    // if we're plotting par 3, add par 4 (both sigmas in 1)
    if( i_p == 3) {      
      g_par[i_p+1][1]->SetLineColor(kBlue);
      g_par[i_p+1][1]->SetMarkerColor(kBlue);
      g_par[i_p+1][1]->SetMarkerStyle(24);
      g_par[i_p+1][1]->SetMarkerSize(.75);
      g_par[i_p+1][1]->Draw("p");
      
      TLegend *leg = new TLegend(0.75, 0.3, 0.9, 0.45);
      leg->SetTextSize(0.03);
      leg->AddEntry(g_par[i_p][1], "#sigma_{1}", "pl");
      leg->AddEntry(g_par[i_p+1][1], "#sigma_{2}", "pl");
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
    
    c->SaveAs(Form("plots/MCMass/par_1_%s.pdf", parlab[i_p].c_str()));
    c->Clear();
  }

  // do the plotting - sigmas for fixed mu, f (MODEL 2)
  TH1F *f2 = c->DrawFrame(pt_min, parmin[3], pt_max, parmax[3]);
  f2->SetXTitle("p_{T} (GeV)");
  f2->SetYTitle(parax[3].c_str());
  f2->GetYaxis()->SetTitleOffset(1.8);
  f2->GetYaxis()->SetLabelOffset(0.01);
  f2->SetTitle(Form("Full %s (#mu, f constant)", partit[3].c_str()));

  c->SetLogy(0);
  
  // free mode always plots points
  g_par[3][2]->SetLineColor(kBlack);
  g_par[3][2]->SetMarkerColor(kBlack);
  g_par[3][2]->SetMarkerStyle(20);
  g_par[3][2]->SetMarkerSize(.75);
  g_par[3][2]->Draw("p");

  g_par[3+1][2]->SetLineColor(kBlue);
  g_par[3+1][2]->SetMarkerColor(kBlue);
  g_par[3+1][2]->SetMarkerStyle(24);
  g_par[3+1][2]->SetMarkerSize(.75);
  g_par[3+1][2]->Draw("p");

  g_par[3][3]->SetLineColor(kRed);
  g_par[3][3]->SetFillColorAlpha(kRed, 0.5);
  g_par[3][3]->Draw("ce3");

  g_par[3+1][3]->SetLineColor(kRed);
  g_par[3+1][3]->SetFillColorAlpha(kRed, 0.5);
  g_par[3+1][3]->Draw("ce3");

  TLegend *leg2 = new TLegend(0.75, 0.3, 0.9, 0.45);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(g_par[3][2], "#sigma_{1}", "pl");
  leg2->AddEntry(g_par[3+1][2], "#sigma_{2}", "pl");
  leg2->Draw();
  
  TLine *l21 = new TLine(46, parmin[3], 46, parmax[3]);
  l21->SetLineColor(kBlack);
  l21->SetLineStyle(kDashed);
  l21->Draw();
  TLine *l22 = new TLine(66, parmin[3], 66, parmax[3]);
  l22->SetLineColor(kBlack);
  l22->SetLineStyle(kDashed);
  l22->Draw();
    
  c->SaveAs("plots/MCMass/par_2.pdf");
  c->Clear();

  // also draw n, alpha
  TH1F *f2na = c->DrawFrame(pt_min, parmin[5], pt_max, parmax[6]);
  f2na->SetXTitle("p_{T} (GeV)");
  f2na->SetYTitle("n, #alpha");
  f2na->GetYaxis()->SetTitleOffset(1.8);
  f2na->GetYaxis()->SetLabelOffset(0.01);
  f2na->SetTitle(Form("Full n, #alpha (#mu, f constant)"));

  c->SetLogy(0);
  
  // free mode always plots points
  g_par[5][2]->SetLineColor(kBlack);
  g_par[5][2]->SetMarkerColor(kBlack);
  g_par[5][2]->SetMarkerStyle(20);
  g_par[5][2]->SetMarkerSize(.75);
  g_par[5][2]->Draw("p");

  g_par[5+1][2]->SetLineColor(kBlue);
  g_par[5+1][2]->SetMarkerColor(kBlue);
  g_par[5+1][2]->SetMarkerStyle(24);
  g_par[5+1][2]->SetMarkerSize(.75);
  g_par[5+1][2]->Draw("p");

  TLegend *leg2na = new TLegend(0.75, 0.5, 0.9, 0.65);
  leg2na->SetTextSize(0.03);
  leg2na->AddEntry(g_par[5][2], "n", "pl");
  leg2na->AddEntry(g_par[5+1][2], "#alpha", "pl");
  leg2na->Draw();
  
  TLine *lna21 = new TLine(46, parmin[5], 46, parmax[6]);
  lna21->SetLineColor(kBlack);
  lna21->SetLineStyle(kDashed);
  lna21->Draw();
  TLine *lna22 = new TLine(66, parmin[5], 66, parmax[6]);
  lna22->SetLineColor(kBlack);
  lna22->SetLineStyle(kDashed);
  lna22->Draw();
    
  c->SaveAs("plots/MCMass/par_2_na.pdf");
  c->Clear();

  // do the plotting - n, alpha for fixed mu, f, sigmas (MODEL 3)
  for(int i_p = 5; i_p < n_p; i_p++) {

    TH1F *fl = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle(parax[i_p].c_str());
    fl->GetYaxis()->SetTitleOffset(1.8);
    fl->GetYaxis()->SetLabelOffset(0.01);
    fl->SetTitle(Form("Full %s (#mu, f constant, #sigma_{1,2} linear)", partit[i_p].c_str()));

    c->SetLogy(0);

    // free mode always plots points
    g_par[i_p][3]->SetLineColor(kBlack);
    g_par[i_p][3]->SetMarkerColor(kBlack);
    g_par[i_p][3]->SetMarkerStyle(20);
    g_par[i_p][3]->SetMarkerSize(.75);
    g_par[i_p][3]->Draw("p");

    if(i_p == 5) {
      g_par[i_p][4]->SetLineColor(kRed);
      g_par[i_p][4]->SetFillColorAlpha(kRed, 0.5);
      g_par[i_p][4]->Draw("ce3");
    }

    TLine *l1 = new TLine(46, parmin[i_p], 46, parmax[i_p]);
    l1->SetLineColor(kBlack);
    l1->SetLineStyle(kDashed);
    l1->Draw();
    TLine *l2 = new TLine(66, parmin[i_p], 66, parmax[i_p]);
    l2->SetLineColor(kBlack);
    l2->SetLineStyle(kDashed);
    l2->Draw();
    
    c->SaveAs(Form("plots/MCMass/par_3_%s.pdf", parlab[i_p].c_str()));
    c->Clear();
  }

  // do the plotting - alpha for fixed mu, n, f, sigmas (MODEL 4)
  TH1F *f4 = c->DrawFrame(pt_min, parmin[6], pt_max, parmax[6]);
  f4->SetXTitle("p_{T} (GeV)");
  f4->SetYTitle(parax[6].c_str());
  f4->GetYaxis()->SetTitleOffset(1.8);
  f4->GetYaxis()->SetLabelOffset(0.01);
  f4->SetTitle("Full #alpha (#mu, f, n constant, #sigma_{1,2} linear)");

  c->SetLogy(0);
  
  // plot the points
  g_par[6][4]->SetLineColor(kBlack);
  g_par[6][4]->SetMarkerColor(kBlack);
  g_par[6][4]->SetMarkerStyle(20);
  g_par[6][4]->SetMarkerSize(.75);
  g_par[6][4]->Draw("p");
  
  g_par[6][5]->SetLineColor(kRed);
  g_par[6][5]->SetFillColorAlpha(kRed, 0.5);
  g_par[6][5]->Draw("ce3");

  TLine *l41 = new TLine(46, parmin[6], 46, parmax[6]);
  l41->SetLineColor(kBlack);
  l41->SetLineStyle(kDashed);
  l41->Draw();
  TLine *l42 = new TLine(66, parmin[6], 66, parmax[6]);
  l42->SetLineColor(kBlack);
  l42->SetLineStyle(kDashed);
  l42->Draw();
    
  c->SaveAs("plots/MCMass/par_4.pdf");
  c->Clear();
  
  // plot all the chi^2 / ndf at the end
  TH1F *fc = c->DrawFrame(pt_min, 0, pt_max, 15);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("#chi^{2}/ndf");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);
  fc->SetTitle(Form("Full #chi^{2}/ndf for different fit models"));

  for(int i = 0; i < n_m; i++) {
    l_chiN[i]->SetLineColor(i+1);
    l_chiN[i]->Draw("lsame");
  }

  TLegend *legc = new TLegend(0.65, 0.7, 0.9, 0.9);
  legc->SetTextSize(0.03);
  legc->AddEntry(l_chiN[0], "All free", "l");
  legc->AddEntry(l_chiN[1], "#mu", "l");
  legc->AddEntry(l_chiN[2], "#mu, f", "l");
  legc->AddEntry(l_chiN[3], "#mu, f, #sigma_{1,2}", "l");
  legc->AddEntry(l_chiN[4], "#mu, f, #sigma_{1,2}, n", "l");
  legc->AddEntry(l_chiN[5], "#mu, f, #sigma_{1,2}, n, #alpha", "l");
  legc->Draw();
  
  c->SaveAs(Form("plots/MCMass/par_chiN.pdf"));
  c->Clear();

  // plot all n and alpha
  for(int i_p = 5; i_p < 7; i_p++) {

    TH1F *flna = c->DrawFrame(pt_min, parmin[i_p], pt_max, parmax[i_p]);
    flna->SetXTitle("p_{T} (GeV)");
    flna->SetYTitle(parax[i_p].c_str());
    flna->GetYaxis()->SetTitleOffset(1.8);
    flna->GetYaxis()->SetLabelOffset(0.01);
    flna->SetTitle(Form("Full %s", partit[i_p].c_str()));

    int c_val = 1;
    for(int i_m = 0; i_m < 4; i_m++) {
      if(i_m !=1) {
	g_par[i_p][i_m]->SetLineColor(c_val);
	g_par[i_p][i_m]->SetMarkerColor(c_val);
	g_par[i_p][i_m]->SetMarkerStyle(20);
	g_par[i_p][i_m]->SetMarkerSize(.75);
	g_par[i_p][i_m]->Draw("p");
	c_val++;
	if(c_val ==3) c_val++;
      }
    }

    TLine *vl = new TLine(pt_min, g_par[i_p][5]->GetY()[0], pt_max, g_par[i_p][5]->GetY()[0]);
    vl->SetLineStyle(kDashed);
    vl->SetLineColor(kBlack);
    vl->Draw();
   
    c->SaveAs(Form("plots/MCMass/na_%s.pdf", parlab[i_p].c_str()));
    c->Clear();
  }

  // also storing chi^2 as table
  ofstream ftex;
  ftex.open("text_output/mfit_chi.tex");
  ftex << "\\begin{tabular}{c|c}\n";
  ftex << "Fit model & $\\chi^2/$ndf \\\\\n";
  ftex << "\\hline\n";
  for(int i = 0; i < n_m; i++) {
    ftex << i << " & " << setprecision(0) << fixed << l_chi[i]->GetY1() << "/" << l_ndf[i]->GetY1()  <<  "\\\\\n";
  }
  ftex << "\\end{tabular}\n";
  ftex.close();
  

}
