#import "mbins.C"
#import "../../Simult/ptbins.C"

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

void plotLtPars_N()
{
  string parlab[] = {"N_NP"};
  string partit[] = {"N_{NP}"};
  string par_unit[] = {" per 1 GeV"};

  // get the 1d fit pars
  double parN[nPtBins][nmBins], eparN[nPtBins][nmBins];
  double m[nmBins], em[nmBins];
  double chi[nPtBins], chiN[nPtBins], chiP[nPtBins];
  int ndf[nPtBins];
  double mults[] = {1};
    
  TFile *fin = new TFile("files/ltfitres_N.root");
  TFitResult *fitres = new TFitResult();
  for(int j = 0; j < nPtBins; j++) {
    fin->GetObject(Form("fitres_%.0f", ptBins[j]), fitres);
    double dpt = ptBins[j+1]-ptBins[j];
    for(int i = 0; i < nmBins; i++) {
      m[i] = 0.5*(m_max[i]+m_min[i]);
      em[i] = 0.5*(m_max[i]-m_min[i]);
      mults[0] = 1./(2.*em[i]);
      mults[0] /= dpt; // scale by pT bin width as well!!
      
      parN[j][i] = 2.*em[i]*fitres->Parameter(0)*exp(-m[i]/fitres->Parameter(1));

      double ln = 1e4;
      double sigp1 = fitres->ParError(0);
      double sigp2 = fitres->ParError(1);
      double cov = fitres->GetCovarianceMatrix()[0][1];

      double devp1 = 2.*em[i]*(fitres->Parameter(0)+sigp1/ln)*exp(-m[i]/fitres->Parameter(1));
      double devp2 = 2.*em[i]*fitres->Parameter(0)*exp(-m[i]/(fitres->Parameter(1)+sigp2/ln));

      double dp_1 = (devp1-parN[j][i])/(sigp1/ln);
      double dp_2 = (devp2-parN[j][i])/(sigp2/ln);
	
      eparN[j][i] = pow(dp_1*sigp1,2) + pow(dp_2*sigp2,2) + 2*dp_1*dp_2*cov;
      eparN[j][i] = sqrt(eparN[j][i]);

      parN[j][i]*=mults[0];
      eparN[j][i]*=mults[0];
    }
    cout << fitres->Status() << " ";
    chi[j] = fitres->Chi2();
    ndf[j] = fitres->Ndf();
    chiN[j] = chi[j]/(double)ndf[j];
    chiP[j] = TMath::Prob(chi[j], ndf[j])*100.;
  }
  cout << endl << endl;
  fin->Close();

  TGraphErrors **g_parN = new TGraphErrors*[nPtBins];

  // plotting all pt bins
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.02);
  c->SetLeftMargin(0.11);
  c->SetTopMargin(0.02);

  c->SetLogy();

  double parmin[] = {1e1};
  double parmax[] = {4e4};
  int col[] = {kBlack, kBlue, kRed+1, kGreen+2};

  TH1F *fpN = c->DrawFrame(m_min[0], parmin[0], m_max[nmBins-1], parmax[0]);
  fpN->SetXTitle("m(#mu#mu) (GeV)");
  fpN->SetYTitle(Form("%s%s", partit[0].c_str(), par_unit[0].c_str()));
  fpN->GetYaxis()->SetTitleOffset(1.5);
  fpN->GetYaxis()->SetLabelOffset(0.01);
  fpN->SetTitle("");

  TLegend *leg = new TLegend(0.4, 0.975-0.04*nPtBins, 0.7, 0.975);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);

  for(int i = 0; i < nPtBins; i++) {
 
    g_parN[i] = new TGraphErrors(nmBins, m, parN[i], em, eparN[i]);
    g_parN[i]->SetMarkerStyle(20);
    if(i > 3) g_parN[i]->SetMarkerStyle(24);
    g_parN[i]->SetMarkerSize(.75);
    g_parN[i]->SetMarkerColor(col[i%4]);
    g_parN[i]->SetLineColor(col[i%4]);
    g_parN[i]->Draw("psame");
    
    leg->AddEntry(g_parN[i], Form("[%.0f,%.0f] GeV", ptBins[i], ptBins[i+1]), "pl");
  }
  int isLog =  1;
   
  leg->Draw();
  c->SaveAs(Form("plots/lifetime_N/par_%s.pdf", parlab[0].c_str()));
  c->Clear();
  
  // finally compare chi^2 results
  double pt[nPtBins], ept[nPtBins], dy[nPtBins];
  for(int i = 0; i < nPtBins; i++){
    pt[i] = 0.5*(ptBins[i+1]+ptBins[i]);
    ept[i] = 0.5*(ptBins[i+1]-ptBins[i]);
    dy[i] = 0;
  }
  TGraphErrors *g_chiN = new TGraphErrors(nPtBins, pt, chiN, ept, dy);
  TGraphErrors *g_chiP = new TGraphErrors(nPtBins, pt, chiP, ept, dy);

  c->SetLogy(0);
  
  TH1F *fcN = c->DrawFrame(ptBins[0]-5, 0, ptBins[nPtBins]+5, 2);
  fcN->SetXTitle("p_{T} (GeV)");
  fcN->SetYTitle(Form("#chi^{2}/ndf"));
  fcN->GetYaxis()->SetTitleOffset(1.5);
  fcN->GetYaxis()->SetLabelOffset(0.01);
  fcN->SetTitle("");
  
  g_chiN->SetMarkerStyle(20);
  g_chiN->SetMarkerSize(.75);
  g_chiN->Draw("psame");

  TLine *l1 = new TLine(ptBins[0]-5, 1, ptBins[nPtBins]+5, 1);
  l1->SetLineStyle(kDashed);
  l1->Draw();
  
  c->SaveAs("plots/lifetime_N/chiN.pdf");
  c->Clear();
  
  TH1F *fcP = c->DrawFrame(ptBins[0]-5, 0, ptBins[nPtBins]+5, 100);
  fcP->SetXTitle("p_{T} (GeV)");
  fcP->SetYTitle(Form("P(#chi^{2}, ndf) (%%)"));
  fcP->GetYaxis()->SetTitleOffset(1.5);
  fcP->GetYaxis()->SetLabelOffset(0.01);
  fcP->SetTitle("");

  g_chiP->SetMarkerStyle(20);
  g_chiP->SetMarkerSize(.75);
  g_chiP->Draw("psame");

  c->SaveAs("plots/lifetime_N/chiP.pdf");
  c->Clear();

  c->Destructor();

}
