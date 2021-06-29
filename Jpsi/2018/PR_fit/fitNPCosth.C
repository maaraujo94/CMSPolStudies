#import "../cosMax/imp_jumpF.C"

// macro to fit the NP/MC distributions
void fitNPCosth()
{
  string lbl = "NP";
  double dev_min = 30, pull_min = 3;

  double N_min = 3e-2, N_max = 5e-1;
  double l2_min = 1.0;
  
  // get the NP/MC 2d maps from the repository
  TH2D *h_cth = new TH2D();
  TFile *fIn = new TFile("files/histoStore.root");
  h_cth = (TH2D*)fIn->Get("ratNPH_ab");
  h_cth->SetDirectory(0);
  fIn->Close();
  
  // get the binning
  int nBinsX = h_cth->GetNbinsX(), nBinsY = h_cth->GetNbinsY();
  const double *yBins = h_cth->GetYaxis()->GetXbins()->GetArray();
  double minX = h_cth->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_cth->GetXaxis()->GetBinUpEdge(nBinsX);

  // get the fit range from our cosmax(pT)
  ifstream in;
  string dataS;
  in.open("../cosMax/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", yBins[0]-10, yBins[nBinsY]+10);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);

  TCanvas *c = new TCanvas("", "", 700, 700);
  
  // define 1d histos for fitting and par arrays
  TH1D *pHist[nBinsY];
  double par[2][nBinsY], epar[2][nBinsY];
  double chi2[nBinsY], ndf[nBinsY], chiP[nBinsY];
  double cMaxVal[nBinsY], pt[nBinsY], ept[nBinsY];

  // cycle over all pT bins
  for(int i = 0; i < nBinsY; i++) {
    
    pt[i] = 0.5*(yBins[i+1]+yBins[i]);
    ept[i] = 0.5*(yBins[i+1]-yBins[i]);
    
    // getting the pT bin projections of Data/MC
    pHist[i] = h_cth->ProjectionX(Form("%s_bin%d_1d", lbl.c_str(), i+1), i+1, i+1);
    pHist[i]->SetTitle(Form("2018 %s/MC bin %d: [%.0f, %.0f] GeV", lbl.c_str(), i+1, yBins[i], yBins[i+1]));
    
    // getting the max costh value for the fit
    double pMin = h_cth->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_cth->GetYaxis()->GetBinUpEdge(i+1);
    cMaxVal[i] = jumpF(cosMax->Integral(pMin, pMax)/(pMax-pMin));

    // define the ratio fit function
    TF1 *f_fit = new TF1("f_fit", "[0]*(1+[1]*x*x)", 0, cMaxVal[i]);
    f_fit->SetParNames("N", "l_2");
    f_fit->SetParameters(pHist[i]->GetBinContent(1)*1.1, -0.1);

    pHist[i]->Fit("f_fit", "R");
    for(int j = 0; j < 2; j++) {
      par[j][i] = f_fit->GetParameter(j);
      epar[j][i] = f_fit->GetParError(j);
    }
    chi2[i] = f_fit->GetChisquare();
    ndf[i] = f_fit->GetNDF();
    chiP[i] = TMath::Prob(chi2[i], ndf[i]);
    
    // plot the fit
    pHist[i]->SetStats(0);
    pHist[i]->SetMinimum(0);
    pHist[i]->SetMaximum(par[0][i]*2.5);
    pHist[i]->SetMarkerColor(kBlack);
    pHist[i]->SetLineColor(kBlack);
    pHist[i]->GetXaxis()->SetTitle("|cos#theta|");
    pHist[i]->Draw();
    f_fit->SetLineColor(kRed);
    f_fit->Draw("same");
    
    TLatex lc;
    lc.SetTextSize(0.03);
    lc.DrawLatex(0.1, par[0][i]*0.3, Form("#chi^{2}/ndf = %.0f/%.0f", chi2[i], ndf[i]));
    lc.DrawLatex(0.1, par[0][i]*0.2, Form("P(#chi^{2},ndf) = %.1f%%", 100*chiP[i]));
    
    c->SaveAs(Form("plots/%s/cth_pt%d.pdf", lbl.c_str(), i+1));
    c->Clear();
      
    // calculating pulls
    double cv[nBinsX], pv[nBinsX], dv[nBinsX];
    for(int i_cos = 0 ; i_cos < nBinsX; i_cos++) {
      double cMin = h_cth->GetXaxis()->GetBinLowEdge(i_cos+1);
      double cMax = h_cth->GetXaxis()->GetBinUpEdge(i_cos+1);
      cv[i_cos] = 0.5 * (cMax+cMin);
      
      double fitv = f_fit->Eval(cv[i_cos]);
      double datav = pHist[i]->GetBinContent(i_cos+1);
      double datau = pHist[i]->GetBinError(i_cos+1);
      if(cv[i_cos] < cMaxVal[i]) {
	pv[i_cos] = (datav-fitv)/datau;
	dv[i_cos] = (datav-fitv)/fitv * 100.;
      }
      else {
	pv[i_cos] = 0;
	dv[i_cos] = 0;
      }
    }
    
    // plot the pulls
    TH1F *fp = c->DrawFrame(minX, -pull_min, maxX, pull_min);
    fp->SetXTitle("|cos#theta|");
    fp->SetYTitle("pulls");
    fp->GetYaxis()->SetTitleOffset(1.3);
    fp->GetYaxis()->SetLabelOffset(0.01);
    fp->SetTitle(Form("2018 %s pulls (%.0f < p_{T} < %.0f GeV)", lbl.c_str(), yBins[i], yBins[i+1]));
    
    TGraph *g_pull = new TGraph(nBinsX, cv, pv);
    g_pull->SetLineColor(kBlack);
    g_pull->SetMarkerColor(kBlack);
    g_pull->SetMarkerStyle(20);
    g_pull->Draw("p");
    
    TLine *clim = new TLine(cMaxVal[i], -pull_min, cMaxVal[i], pull_min);
    clim->SetLineColor(kRed);
    clim->SetLineStyle(kDashed);
    clim->Draw();

    TLine *zero = new TLine(minX, 0, maxX, 0);
    zero->SetLineColor(kBlack);
    zero->SetLineStyle(kDashed);
    zero->Draw();
    
    c->SaveAs(Form("plots/%s/pulls_pt%d.pdf", lbl.c_str(), i+1));
    c->Clear();

    // plot the deviations
    TH1F *fd = c->DrawFrame(minX, -dev_min, maxX, dev_min);
    fd->SetXTitle("|cos#theta|");
    fd->SetYTitle("relative difference (%)");
    fd->GetYaxis()->SetTitleOffset(1.3);
    fd->GetYaxis()->SetLabelOffset(0.01);
    fd->SetTitle(Form("2018 %s rel. difference (%.0f < p_{T} < %.0f GeV)", lbl.c_str(), yBins[i], yBins[i+1]));

    TGraph *g_dev = new TGraph(nBinsX, cv, dv);
    g_dev->SetLineColor(kBlack);
    g_dev->SetMarkerColor(kBlack);
    g_dev->SetMarkerStyle(20);
    g_dev->Draw("p");
    
    TLine *climd = new TLine(cMaxVal[i], -dev_min, cMaxVal[i], dev_min);
    climd->SetLineColor(kRed);
    climd->SetLineStyle(kDashed);
    climd->Draw();
    
    zero->Draw();
    
    c->SaveAs(Form("plots/%s/dev_pt%d.pdf", lbl.c_str(), i+1));
    c->Clear();
  }

  // making TGraphs for all the fit parameters
  TGraphErrors *graph_N  = new TGraphErrors(nBinsY, pt, par[0], ept, epar[0]);
  TGraphErrors *graph_l2 = new TGraphErrors(nBinsY, pt, par[1], ept, epar[1]);
  
  TGraph *graph_C  = new TGraph(nBinsY, pt, chi2);
  TGraph *graph_ND = new TGraph(nBinsY, pt, ndf);
  TGraph *graph_P  = new TGraph(nBinsY, pt, chiP);
      
  c->SetLogy(0);
  TH1F *fc = c->DrawFrame(yBins[0]-5, 0, yBins[nBinsY]+5, 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("P(#chi^{2}, ndf)");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);
  fc->SetTitle(Form("2018 %s P(#chi^{2}, ndf)", lbl.c_str()));
      
  graph_P->SetLineColor(kBlack);
  graph_P->SetMarkerColor(kBlack);
  graph_P->SetMarkerStyle(20);
  graph_P->SetMarkerSize(1.5);
  graph_P->Draw("p");
      
  c->SaveAs(Form("plots/%s/par_chiP.pdf", lbl.c_str()));
  c->Clear();
      
  c->SetLogy();
  TH1F *fN = c->DrawFrame(yBins[0]-5, N_min, yBins[nBinsY]+5, N_max);
  fN->SetXTitle("p_{T} (GeV)");
  fN->SetYTitle("N");
  fN->GetYaxis()->SetTitleOffset(1.3);
  fN->GetYaxis()->SetLabelOffset(0.01);
  fN->SetTitle(Form("2018 %s N", lbl.c_str()));
      
  graph_N->SetLineColor(kBlack);
  graph_N->SetMarkerColor(kBlack);
  graph_N->SetMarkerStyle(20);
  graph_N->SetMarkerSize(0.75);
  graph_N->Draw("p");
      
  c->SaveAs(Form("plots/%s/par_N.pdf", lbl.c_str()));
  c->Clear();
      
  c->SetLogy(0);
  TH1F *fl2 = c->DrawFrame(yBins[0]-5, -l2_min, yBins[nBinsY]+5, l2_min);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#lambda_{NP}");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle(Form("2018 %s #lambda_{NP}", lbl.c_str()));
      
  graph_l2->SetLineColor(kBlack);
  graph_l2->SetMarkerColor(kBlack);
  graph_l2->SetMarkerStyle(20);
  graph_l2->SetMarkerSize(0.75);
  graph_l2->Draw("p");
            
  TLine *zero = new TLine(yBins[0]-5, 0, yBins[nBinsY]+5, 0);
  zero->SetLineStyle(kDashed);
  zero->SetLineColor(kBlack);
  zero->Draw();

  c->SaveAs(Form("plots/%s/par_l2.pdf", lbl.c_str()));
  c->Clear();
  
  // storing the fit results in tex format
  ofstream fout;
  fout.open(Form("text_output/cos_%s.tex", lbl.c_str()));
  fout << "\\begin{tabular}{c||c|c|c}\n";
  fout << "$\\pt$ (GeV) & $N$ & $\\lambda_{NP}$ & $\\chi^2$/ndf  \\\\\n";
  fout << "\\hline\n";
  for(int i = 0; i < nBinsY; i++) {
    double pMin = h_cth->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_cth->GetYaxis()->GetBinUpEdge(i+1);
    fout << Form("$[%.0f, %.0f]$", pMin, pMax);
    for(int ip = 0; ip < 2; ip++){
      if(epar[ip][i] > 0) {
	int p_norm = 1.;
	if(epar[ip][i] < 1 ) 
	  p_norm = ceil(-log10(epar[ip][i]))+1;
	fout << " & " <<  setprecision(p_norm) << fixed << par[ip][i] << " $\\pm$ " << epar[ip][i];
      }
      else
	fout << " & " <<  setprecision(2) << fixed << par[ip][i];
    }
    fout << Form(" & %.0f", chi2[i]) << "/" << Form("%.0f", ndf[i]) << "\\\\\n";
  }
  fout << "\\end{tabular}\n";
  fout.close();

  // storing the TGraphs
  TFile *fOut = new TFile(Form("files/%s_fitres.root", lbl.c_str()), "recreate");
  graph_N->SetName(Form("graph_N"));
  graph_N->Write();
  graph_l2->SetName(Form("graph_l2"));
  graph_l2->Write();
  graph_C->SetName(Form("graph_chi"));
  graph_C->Write();
  graph_ND->SetName(Form("graph_ndf"));
  graph_ND->Write();
  graph_P->SetName(Form("graph_chiP"));
  graph_P->Write();
  fOut->Close();
  
  cout << endl << "finished fitting " << lbl << endl << endl;
  c->Destructor();
}
