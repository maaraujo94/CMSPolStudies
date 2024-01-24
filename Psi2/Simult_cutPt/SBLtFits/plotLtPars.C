#import "mbins.C"

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

void plotLtPars()
{
  const int n_p = 4;
  string parlab[] = {"N_NP", "f1", "tnp1", "tnp2"};
  string partit[] = {"N_{NP}", "f_{1}", "t_{NP}", "t_{NP_{2}}"};
  string par_unit[] = {" per 1 GeV", " (%)", " (#mum)", " (#mum)"};

  double chisq[4];
  int ndf[4];

  // get the 2d fit pars - fix tnp2
  TGraphErrors **g1_par_2d = new TGraphErrors*[n_p];
  TFile *fin1_2d = new TFile("files/ltfitres2d_tnp2.root");
  for(int i = 0; i < n_p; i++) {
    g1_par_2d[i] = (TGraphErrors*)fin1_2d->Get(Form("fit_%s", parlab[i].c_str()));
  }
  chisq[0] = ((TFitResult*)fin1_2d->Get("fitres"))->Chi2();
  ndf[0] = ((TFitResult*)fin1_2d->Get("fitres"))->Ndf();
  fin1_2d->Close();
  
  // get the 2d fit pars - fix tnp1, tnp2
  TGraphErrors **g2_par_2d = new TGraphErrors*[n_p];
  TFile *fin2_2d = new TFile("files/ltfitres2d_tnp12.root");
  for(int i = 0; i < n_p; i++) {
    g2_par_2d[i] = (TGraphErrors*)fin2_2d->Get(Form("fit_%s", parlab[i].c_str()));
  }
  chisq[1] = ((TFitResult*)fin2_2d->Get("fitres"))->Chi2();
  ndf[1] = ((TFitResult*)fin2_2d->Get("fitres"))->Ndf();
  fin2_2d->Close();
  
  // get the 2d fit pars - fix tnp1, tnp2, f
  TGraphErrors **g3_par_2d = new TGraphErrors*[n_p];
  TFile *fin3_2d = new TFile("files/ltfitres2d_tnp12_f.root");
  for(int i = 0; i < n_p; i++) {
    g3_par_2d[i] = (TGraphErrors*)fin3_2d->Get(Form("fit_%s", parlab[i].c_str()));
  }
  chisq[2] = ((TFitResult*)fin3_2d->Get("fitres"))->Chi2();
  ndf[2] = ((TFitResult*)fin3_2d->Get("fitres"))->Ndf();
  fin3_2d->Close();
  
  // get the 2d fit pars - fix tnp1, tnp2, N, f
  TGraphErrors **g4_par_2d = new TGraphErrors*[n_p];
  TFile *fin4_2d = new TFile("files/ltfitres2d_tnp12_N_f.root");
  for(int i = 0; i < n_p; i++) {
    g4_par_2d[i] = (TGraphErrors*)fin4_2d->Get(Form("fit_%s", parlab[i].c_str()));
  }
  chisq[3] = ((TFitResult*)fin4_2d->Get("fitres"))->Chi2();
  ndf[3] = ((TFitResult*)fin4_2d->Get("fitres"))->Ndf();
  fin4_2d->Close();
  
  // get the 1d fit pars
  double par[n_p][nmBins], epar[n_p][nmBins];
  double m[nmBins], em[nmBins];
  int stat[nmBins];
  
  TFile *fin = new TFile("files/ltfitres.root");
  double mults[] = {1, 100, 1e3, 1e3};
  TFitResult *fitres = new TFitResult();
  for(int i = 0; i < nmBins; i++) {
    fin->GetObject(Form("fitres_%s", lbl[i].c_str()), fitres);
    m[i] = 0.5*(m_max[i]+m_min[i]);
    em[i] = 0.5*(m_max[i]-m_min[i]);
    mults[0] = 1./(2.*em[i]);
    for(int j = 0; j < n_p; j++) {
      par[j][i] = fitres->GetParams()[j]*mults[j];
      epar[j][i] = fitres->GetErrors()[j]*mults[j];
      // also scaling the 2d fit pars
      g1_par_2d[j]->GetY()[i]*=mults[j];
      g1_par_2d[j]->GetEY()[i]*=mults[j];
      g2_par_2d[j]->GetY()[i]*=mults[j];
      g2_par_2d[j]->GetEY()[i]*=mults[j];
      g3_par_2d[j]->GetY()[i]*=mults[j];
      g3_par_2d[j]->GetEY()[i]*=mults[j];
      g4_par_2d[j]->GetY()[i]*=mults[j];
      g4_par_2d[j]->GetEY()[i]*=mults[j];
    }
    stat[i] = fitres->Status();
  }
  fin->Close();

  TGraphErrors **g_par = new TGraphErrors*[n_p];

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.11);

  double parmin[] = {6e4, 0,   0,    0};
  double parmax[] = {6e5, 100, 500, 1000};

  // compare 1d to 2d with fixed tnp2
  for(int j = 0; j < n_p; j++) {

    if(j == 0 ) c->SetLogy();
    else c->SetLogy(0);
    if(j == 3) continue;

    TH1F *fp = c->DrawFrame(m_min[0], parmin[j], m_max[nmBins-1], parmax[j]);
    fp->SetXTitle("m(#mu#mu) (GeV)");
    fp->SetYTitle(Form("%s%s", partit[j].c_str(), par_unit[j].c_str()));
    fp->GetYaxis()->SetTitleOffset(1.5);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("%s vs m(#mu#mu)", partit[j].c_str()));

    g_par[j] = new TGraphErrors(nmBins, m, par[j], em, epar[j]);
    g_par[j]->SetMarkerStyle(24);
    g_par[j]->SetMarkerSize(.75);
    g_par[j]->SetMarkerColor(kBlack);
    g_par[j]->SetLineColor(kBlack);
    g_par[j]->Draw("psame");

    g1_par_2d[j]->SetMarkerStyle(20);
    g1_par_2d[j]->SetMarkerSize(.75);
    g1_par_2d[j]->SetMarkerColor(kGreen+3);
    g1_par_2d[j]->SetLineColor(kGreen+3);
    if(j == 1) g1_par_2d[j]->Draw("psame");
    else {
      g1_par_2d[j]->SetFillColorAlpha(kGreen+3,0.5);
      g1_par_2d[j]->Draw("psame");
    }

    if(j == 2) {
      int k = j+1;
      g_par[k] = new TGraphErrors(nmBins, m, par[k], em, epar[k]);
      g_par[k]->SetMarkerStyle(25);
      g_par[k]->SetMarkerSize(.75);
      g_par[k]->SetMarkerColor(kBlack);
      g_par[k]->SetLineColor(kBlack);
      g_par[k]->Draw("psame");

      g1_par_2d[k]->SetMarkerStyle(21);
      g1_par_2d[k]->SetMarkerSize(.75);
      g1_par_2d[k]->SetFillColorAlpha(kGreen+3,0.5);
      g1_par_2d[k]->SetMarkerColor(kGreen+3);
      g1_par_2d[k]->SetLineColor(kGreen+3);
      g1_par_2d[k]->Draw("le3");
      
      TF1 *fcon = new TF1("fcon", "[0]", m_min[0], m_max[nmBins-1]);
      fcon->SetParameter(0, 50);
      fcon->SetLineColor(kBlack);
      fcon->SetLineStyle(kDashed);
      g_par[k]->Fit(fcon);

      TLatex lct;
      lct.SetTextSize(0.04);
      double xp = getPos(m_min[0], m_max[nmBins-1], 0.05, 0);
      double yp = getPos(parmin[j], parmax[j], 0.9, 0);
      lct.DrawLatex(xp, yp, "#bf{t_{NP_{1}}}");
      yp = getPos(parmin[j], parmax[j], 0.2, 0);
      lct.DrawLatex(xp, yp, "#bf{t_{NP_{2}}}");
    }
    
    int isLog = 0;
    if(j == 0) isLog = 1;
    
    TLatex lc;
    lc.SetTextSize(0.03);
    double xp = getPos(m_min[0], m_max[nmBins-1], 0.75, 0);
    double yp = getPos(parmin[j], parmax[j], 0.95, isLog);
    lc.DrawLatex(xp, yp, "1d fit");
    lc.SetTextColor(kGreen+3);
    yp = getPos(parmin[j], parmax[j], 0.95-0.07, isLog);
    lc.DrawLatex(xp, yp, "2d fit (fixed t_{NP2})");

    c->SaveAs(Form("plots/par1_%s.pdf", parlab[j].c_str()));
    c->Clear();
  }

  // compare 2d with fixed tnp2 to 2d with fixed tnp1, tnp2
  for(int j = 0; j < n_p; j++) {

    if(j == 0 ) c->SetLogy();
    else c->SetLogy(0);
    if(j == 3) continue;

    TH1F *fp = c->DrawFrame(m_min[0], parmin[j], m_max[nmBins-1], parmax[j]);
    fp->SetXTitle("m(#mu#mu) (GeV)");
    fp->SetYTitle(Form("%s%s", partit[j].c_str(), par_unit[j].c_str()));
    fp->GetYaxis()->SetTitleOffset(1.5);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("%s vs m(#mu#mu)", partit[j].c_str()));

    g1_par_2d[j]->SetMarkerStyle(20);
    g1_par_2d[j]->SetMarkerSize(.75);
    g1_par_2d[j]->SetMarkerColor(kGreen+3);
    g1_par_2d[j]->SetLineColor(kGreen+3);
    g1_par_2d[j]->Draw("psame");

    g2_par_2d[j]->SetMarkerStyle(20);
    g2_par_2d[j]->SetMarkerSize(.75);
    g2_par_2d[j]->SetMarkerColor(kRed);
    g2_par_2d[j]->SetLineColor(kRed);
    if(j == 2) {
      g2_par_2d[j]->SetFillColorAlpha(kRed,0.5);
      g2_par_2d[j]->Draw("le3");
    }
    else g2_par_2d[j]->Draw("psame");

    if(j == 2) {
      int k = j+1;
      g1_par_2d[k]->SetMarkerStyle(21);
      g1_par_2d[k]->SetMarkerSize(.75);
      g1_par_2d[k]->SetFillColorAlpha(kGreen+3,0.5);
      g1_par_2d[k]->SetMarkerColor(kGreen+3);
      g1_par_2d[k]->SetLineColor(kGreen+3);
      g1_par_2d[k]->Draw("le3");

      g2_par_2d[k]->SetMarkerStyle(21);
      g2_par_2d[k]->SetMarkerSize(.75);
      g2_par_2d[k]->SetFillColorAlpha(kRed,0.5);
      g2_par_2d[k]->SetMarkerColor(kRed);
      g2_par_2d[k]->SetLineColor(kRed);
      g2_par_2d[k]->Draw("le3");
 
      TF1 *fcon = new TF1("fcon", "[0]", m_min[0], m_max[nmBins-1]);
      fcon->SetParameter(0, 400);
      fcon->SetLineColor(kGreen+3);
      fcon->SetLineStyle(kDashed);
      g1_par_2d[j]->Fit(fcon);

      TLatex lct;
      lct.SetTextSize(0.04);
      double xp = getPos(m_min[0], m_max[nmBins-1], 0.05, 0);
      double yp = getPos(parmin[j], parmax[j], 0.9, 0);
      lct.DrawLatex(xp, yp, "#bf{t_{NP_{1}}}");
      yp = getPos(parmin[j], parmax[j], 0.2, 0);
      lct.DrawLatex(xp, yp, "#bf{t_{NP_{2}}}");
    }
    
    int isLog = 0;
    if(j == 0) isLog = 1;
    
    TLatex lc;
    lc.SetTextSize(0.03);
    double xp = getPos(m_min[0], m_max[nmBins-1], 0.75, 0);
    double yp = getPos(parmin[j], parmax[j], 0.95, isLog);
    lc.SetTextColor(kGreen+3);
    lc.DrawLatex(xp, yp, "2d fit (free t_{NP1})");
    yp = getPos(parmin[j], parmax[j], 0.95-0.07, isLog);
    lc.SetTextColor(kRed);
    lc.DrawLatex(xp, yp, "2d fit (fixed t_{NP1})");
    
    c->SaveAs(Form("plots/par2_%s.pdf", parlab[j].c_str()));
    c->Clear();
  }
  
  // compare 2d with fixed tnp1, tnp2 to 2d with fixed tnp1, tnp2, f
  for(int j = 0; j < n_p; j++) {

    if(j == 0 ) c->SetLogy();
    else c->SetLogy(0);
    if(j == 3) continue;

    TH1F *fp = c->DrawFrame(m_min[0], parmin[j], m_max[nmBins-1], parmax[j]);
    fp->SetXTitle("m(#mu#mu) (GeV)");
    fp->SetYTitle(Form("%s%s", partit[j].c_str(), par_unit[j].c_str()));
    fp->GetYaxis()->SetTitleOffset(1.5);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("%s vs m(#mu#mu)", partit[j].c_str()));

    g3_par_2d[j]->SetMarkerStyle(20);
    g3_par_2d[j]->SetMarkerSize(.75);
    g3_par_2d[j]->SetMarkerColor(kBlue);
    g3_par_2d[j]->SetLineColor(kBlue);
    if(j != 0) {
      g3_par_2d[j]->SetFillColorAlpha(kBlue,0.5);
      g3_par_2d[j]->Draw("le3");
    }
    else g3_par_2d[j]->Draw("psame");

    g2_par_2d[j]->SetMarkerStyle(20);
    g2_par_2d[j]->SetMarkerSize(.75);
    g2_par_2d[j]->SetMarkerColor(kRed);
    g2_par_2d[j]->SetLineColor(kRed);
    //   if(j==0 )  g2_par_2d[j]->GetFunction("fexp")->SetBit(TF1::kNotDraw);
    if(j == 2) {
      g2_par_2d[j]->SetFillColorAlpha(kRed,0.5);
      g2_par_2d[j]->Draw("le3");
    }
    else g2_par_2d[j]->Draw("psame");

    if(j == 2) {
      int k = j+1;
      g3_par_2d[k]->SetMarkerStyle(21);
      g3_par_2d[k]->SetMarkerSize(.75);
      g3_par_2d[k]->SetFillColorAlpha(kBlue,0.5);
      g3_par_2d[k]->SetMarkerColor(kBlue);
      g3_par_2d[k]->SetLineColor(kBlue);
      g3_par_2d[k]->Draw("le3");

      g2_par_2d[k]->SetMarkerStyle(21);
      g2_par_2d[k]->SetMarkerSize(.75);
      g2_par_2d[k]->SetFillColorAlpha(kRed,0.5);
      g2_par_2d[k]->SetMarkerColor(kRed);
      g2_par_2d[k]->SetLineColor(kRed);
      g2_par_2d[k]->Draw("le3");
 
      TLatex lct;
      lct.SetTextSize(0.04);
      double xp = getPos(m_min[0], m_max[nmBins-1], 0.05, 0);
      double yp = getPos(parmin[j], parmax[j], 0.9, 0);
      lct.DrawLatex(xp, yp, "#bf{t_{NP_{1}}}");
      yp = getPos(parmin[j], parmax[j], 0.2, 0);
      lct.DrawLatex(xp, yp, "#bf{t_{NP_{2}}}");
    }
    
    if(j==1) {
      TF1 *flin = new TF1("flin", "[0]+x*[1]", m_min[0], m_max[nmBins-1]);
      flin->SetParameters(200, -70);
      flin->SetLineColor(kRed);
      flin->SetLineStyle(kDashed);
      g2_par_2d[j]->Fit(flin);
    }

    int isLog = 0;
    if(j == 0) isLog = 1;
    
    TLatex lc;
    lc.SetTextSize(0.03);
    double xp = getPos(m_min[0], m_max[nmBins-1], 0.75, 0);
    double yp = getPos(parmin[j], parmax[j], 0.95, isLog);
    lc.SetTextColor(kRed);
    lc.DrawLatex(xp, yp, "2d fit (free f_{1})");
    yp = getPos(parmin[j], parmax[j], 0.95-0.07, isLog);
    lc.SetTextColor(kBlue);
    lc.DrawLatex(xp, yp, "2d fit (fixed f_{1})");
    
    c->SaveAs(Form("plots/par3_%s.pdf", parlab[j].c_str()));
    c->Clear();
  }

  // compare 2d with fixed tnp1, tnp2, f to 2d with fixed tnp1, tnp2, N, f
  for(int j = 0; j < n_p; j++) {

    if(j == 0 ) c->SetLogy();
    else c->SetLogy(0);
    if(j == 3) continue;

    TH1F *fp = c->DrawFrame(m_min[0], parmin[j], m_max[nmBins-1], parmax[j]);
    fp->SetXTitle("m(#mu#mu) (GeV)");
    fp->SetYTitle(Form("%s%s", partit[j].c_str(), par_unit[j].c_str()));
    fp->GetYaxis()->SetTitleOffset(1.5);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("%s vs m(#mu#mu)", partit[j].c_str()));

    g3_par_2d[j]->SetMarkerStyle(20);
    g3_par_2d[j]->SetMarkerSize(.75);
    g3_par_2d[j]->SetMarkerColor(kBlue);
    g3_par_2d[j]->SetLineColor(kBlue);
    if(j != 0) {
      g3_par_2d[j]->SetFillColorAlpha(kBlue,0.5);
      g3_par_2d[j]->Draw("le3");
    }
    else g3_par_2d[j]->Draw("psame");

    g4_par_2d[j]->SetMarkerStyle(20);
    g4_par_2d[j]->SetMarkerSize(.75);
    g4_par_2d[j]->SetMarkerColor(kBlack);
    g4_par_2d[j]->SetLineColor(kBlack);
    g4_par_2d[j]->SetFillColorAlpha(kBlack,0.5);
    g4_par_2d[j]->Draw("le3");
    
    if(j == 2) {
      int k = j+1;
      g3_par_2d[k]->SetMarkerStyle(21);
      g3_par_2d[k]->SetMarkerSize(.75);
      g3_par_2d[k]->SetFillColorAlpha(kBlue,0.5);
      g3_par_2d[k]->SetMarkerColor(kBlue);
      g3_par_2d[k]->SetLineColor(kBlue);
      g3_par_2d[k]->Draw("le3");

      g4_par_2d[k]->SetMarkerStyle(21);
      g4_par_2d[k]->SetMarkerSize(.75);
      g4_par_2d[k]->SetFillColorAlpha(kBlack,0.5);
      g4_par_2d[k]->SetMarkerColor(kBlack);
      g4_par_2d[k]->SetLineColor(kBlack);
      g4_par_2d[k]->Draw("le3");
 
      TLatex lct;
      lct.SetTextSize(0.04);
      double xp = getPos(m_min[0], m_max[nmBins-1], 0.05, 0);
      double yp = getPos(parmin[j], parmax[j], 0.9, 0);
      lct.DrawLatex(xp, yp, "#bf{t_{NP_{1}}}");
      yp = getPos(parmin[j], parmax[j], 0.2, 0);
      lct.DrawLatex(xp, yp, "#bf{t_{NP_{2}}}");
    }

    if(j==0) {
      TF1 *fexp = new TF1("fexp", "[0]*exp(-x/[1])", m_min[0], m_max[nmBins-1]);
      fexp->SetParameters(1e7, 1.);
      fexp->SetLineColor(kBlue);
      fexp->SetLineStyle(kDashed);
      g3_par_2d[j]->Fit(fexp);
    }


    int isLog = 0;
    if(j == 0) isLog = 1;
    
    TLatex lc;
    lc.SetTextSize(0.03);
    double xp = getPos(m_min[0], m_max[nmBins-1], 0.75, 0);
    double yp = getPos(parmin[j], parmax[j], 0.95, isLog);
    lc.SetTextColor(kBlue);
    lc.DrawLatex(xp, yp, "2d fit (free N_{NP})");
    yp = getPos(parmin[j], parmax[j], 0.95-0.07, isLog);
    lc.SetTextColor(kBlack);
    lc.DrawLatex(xp, yp, "2d fit (fixed N_{NP})");
    
    c->SaveAs(Form("plots/par4_%s.pdf", parlab[j].c_str()));
    c->Clear();
  }

  // finally compare chi^2 results
  double chiN[4], chiP[4], px[4], dx[4], dy[4];
  for(int i = 0; i < 4; i++){
    px[i] = i+1;
    dx[i] = 0.5;
    dy[i] = 0;
    chiN[i] = chisq[i]/(double)ndf[i];
    chiP[i] = TMath::Prob(chisq[i], ndf[i])*100.;
  }
  TGraphErrors *g_chiN = new TGraphErrors(4, px, chiN, dx, dy);
  TGraphErrors *g_chiP = new TGraphErrors(4, px, chiP, dx, dy);
  
  TH1F *fcN = c->DrawFrame(0.1, 0, 4.9, 2);
  fcN->SetXTitle("fit model");
  fcN->SetYTitle(Form("#chi^{2}/ndf"));
  fcN->GetYaxis()->SetTitleOffset(1.5);
  fcN->GetYaxis()->SetLabelOffset(0.01);
  fcN->GetXaxis()->SetNdivisions(5);
  fcN->SetTitle(Form("#chi^{2}/ndf per fit"));

  g_chiN->SetMarkerStyle(20);
  g_chiN->SetMarkerSize(.75);
  g_chiN->Draw("psame");

  c->SaveAs("plots/chiN.pdf");
  c->Clear();

  TH1F *fcP = c->DrawFrame(0.1, 0, 4.9, 100);
  fcP->SetXTitle("fit model");
  fcP->SetYTitle(Form("P(#chi^{2}, ndf) (%%)"));
  fcP->GetYaxis()->SetTitleOffset(1.5);
  fcP->GetYaxis()->SetLabelOffset(0.01);
  fcP->GetXaxis()->SetNdivisions(5);
  fcP->SetTitle(Form("P(#chi^{2}, ndf) per fit"));

  g_chiP->SetMarkerStyle(20);
  g_chiP->SetMarkerSize(.75);
  g_chiP->Draw("psame");

  c->SaveAs("plots/chiP.pdf");
  c->Clear();

  c->Destructor();
}
