// macro to fit the sideband/MC distributions
void fitBkg()
{
  // get the MC and SB 2d maps from the repository
  TFile *fIn = new TFile("files/ratioHist.root");
  TH2D* bL_R = (TH2D*)fIn->Get("ratioHist_ab_L");
  TH2D* bR_R = (TH2D*)fIn->Get("ratioHist_ab_R");
  bL_R->SetDirectory(0);
  bR_R->SetDirectory(0);
  fIn->Close();
  
  // define the ratio fit function
  TF1 *f_fit = new TF1("f_fit", "[0]*(1+[1]*x*x+[2]*pow(x,4))", 0, 0.9);
  f_fit->SetParNames("N", "l_2", "l_4");
  f_fit->SetParameters(0.1, 1., 0.5);

  // get the binning
  int nBinsX = bL_R->GetNbinsX(), nBinsY = bL_R->GetNbinsY();
  const double *yBins = bL_R->GetYaxis()->GetXbins()->GetArray();
  double minX = bL_R->GetXaxis()->GetBinLowEdge(1);
  double maxX = bL_R->GetXaxis()->GetBinUpEdge(nBinsX);

  // get the fit range from our cosmax(pT)
  ifstream in;
  string dataS;
  in.open("text_output/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", yBins[0]-10, yBins[nBinsY]+10);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);

  // define 1d histos for fitting and par arrays
  TH1D *pHistL[nBinsY];
  TH1D *pHistR[nBinsY];
  double par_L[3][nBinsY], par_R[3][nBinsY];
  double epar_L[3][nBinsY], epar_R[3][nBinsY];
  double chi2[2][nBinsY], ndf[2][nBinsY], chiP[2][nBinsY], cMaxVal[nBinsY];
  double pt[nBinsY], ept[nBinsY];

  TCanvas *c = new TCanvas("", "", 700, 700);
  int n_fit = 3;
  for(int i_fit = 0; i_fit < n_fit; i_fit++)
    {
      TH1D **pullL = new TH1D*[nBinsY];
      TH1D **pullR = new TH1D*[nBinsY];
      TH1D **devL = new TH1D*[nBinsY];
      TH1D **devR = new TH1D*[nBinsY];
  
      // cycle over all pT bins
      for(int i = 0; i < nBinsY; i++) {
	pt[i] = 0.5*(yBins[i+1]+yBins[i]);
	ept[i] = 0.5*(yBins[i+1]-yBins[i]);

	// getting the pT bin projections of SB/MC
	pHistL[i] = bL_R->ProjectionX(Form("LSB_bin%d_1d", i+1), i+1, i+1);
	pHistL[i]->SetTitle(Form("2018 LSB/MC bin %d: [%.0f, %.0f] GeV", i+1, yBins[i], yBins[i+1]));
	pHistR[i] = bR_R->ProjectionX(Form("RSB_bin%d_1d", i+1), i+1, i+1);
	pHistR[i]->SetTitle(Form("2018 RSB/MC bin %d: [%.0f, %.0f] GeV", i+1, yBins[i], yBins[i+1]));

	// initializing the pull and dev histos
	pullL[i] = new TH1D(Form("pullL_%d_%d", i_fit, i), Form("2018 fit pulls (%.0f < p_{T} < %.0f GeV)", yBins[i], yBins[i+1]), nBinsX, minX, maxX);
	pullR[i] = new TH1D(Form("pullR_%d_%d", i_fit, i), Form("2018 fit pulls (%.0f < p_{T} < %.0f GeV)", yBins[i], yBins[i+1]), nBinsX, minX, maxX);
	devL[i] = new TH1D(Form("devL_%d_%d", i_fit, i), Form("2018 fit deviation (%.0f < p_{T} < %.0f GeV)", yBins[i], yBins[i+1]), nBinsX, minX, maxX);
	devR[i] = new TH1D(Form("devR_%d_%d", i_fit, i), Form("2018 fit deviation (%.0f < p_{T} < %.0f GeV)", yBins[i], yBins[i+1]), nBinsX, minX, maxX);

	// getting the max costh value for the fit
	double pMin = bL_R->GetYaxis()->GetBinLowEdge(i+1);
	double pMax = bL_R->GetYaxis()->GetBinUpEdge(i+1);
    
	cMaxVal[i] = cosMax->Integral(pMin, pMax)/(pMax-pMin);
	double cR = floor(cMaxVal[i]*10.)/10.;
	if(cMaxVal[i]-cR>0.05) cR += 0.05;
	f_fit->SetRange(0, cR);

	// fitting the LSB
	f_fit->SetParameter(0, pHistL[i]->GetMaximum()/2.);
	if(i_fit > 0)
	  f_fit->FixParameter(2,1);
	if(i_fit > 1)
	  f_fit->FixParameter(1, 0.8);
	pHistL[i]->Fit("f_fit", "R");
	for(int j = 0; j < 3; j++) {
	  par_L[j][i] = f_fit->GetParameter(j);
	  epar_L[j][i] = f_fit->GetParError(j);
	  chi2[0][i] = f_fit->GetChisquare();
	  ndf[0][i] = f_fit->GetNDF();
	  chiP[0][i] = TMath::Prob(chi2[0][i], ndf[0][i]);
	}
	
	for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
	  double cMin = bL_R->GetXaxis()->GetBinLowEdge(i_cos+1);
	  double cMax = bL_R->GetXaxis()->GetBinUpEdge(i_cos+1);
	  double cos = 0.5 * (cMax+cMin);
      
	  double data_val = pHistL[i]->GetBinContent(i_cos+1);
	  double data_err = pHistL[i]->GetBinError(i_cos+1);
	  double fit_val = f_fit->Eval(cos);

	  if(cos < cR) {
	    pullL[i]->SetBinContent(i_cos+1, (data_val-fit_val)/data_err);
	    pullL[i]->SetBinError(i_cos+1, 0);
	    devL[i]->SetBinContent(i_cos+1, (data_val-fit_val)/fit_val);
	    devL[i]->SetBinError(i_cos+1, 0);
	  }
	}
    
	pHistL[i]->SetStats(0);
	pHistL[i]->SetMinimum(0);
	pHistL[i]->SetMaximum(par_L[0][i]*2.5);
	pHistL[i]->SetMarkerColor(kBlue);
	pHistL[i]->SetLineColor(kBlue);
	pHistL[i]->GetXaxis()->SetTitle("|cos#theta|");
	pHistL[i]->Draw();
	f_fit->SetLineColor(kRed);
	f_fit->Draw("same");
    
	TLatex lcL;
	lcL.SetTextSize(0.03);
	lcL.DrawLatex(0.1, par_L[0][i]*0.3, Form("#chi^{2}/ndf = %.0f/%.0f", chi2[0][i], ndf[0][i]));
	lcL.DrawLatex(0.1, par_L[0][i]*0.2, Form("P(#chi^{2},ndf) = %.1f%%", 100*chiP[0][i]));
    
	c->SaveAs(Form("plots/bkgSub/LSB_%d_%d.pdf", i_fit+1, i+1));
	c->Clear();

	// fitting the RSB
	f_fit->SetParameters(0.1, 1., 0.5);
	f_fit->SetParameter(0, pHistR[i]->GetMaximum()/2.);
	if(i_fit > 0)
	  f_fit->FixParameter(2,1);
	if(i_fit > 1)
	  f_fit->FixParameter(1, 1);
	pHistR[i]->Fit("f_fit", "R");
	for(int j = 0; j < 3; j++) {
	  par_R[j][i] = f_fit->GetParameter(j);
	  epar_R[j][i] = f_fit->GetParError(j);
	  chi2[1][i] = f_fit->GetChisquare();
	  ndf[1][i] = f_fit->GetNDF();
	  chiP[1][i] = TMath::Prob(chi2[1][i], ndf[1][i]);
	}
	
	for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
	  double cMin = bL_R->GetXaxis()->GetBinLowEdge(i_cos+1);
	  double cMax = bL_R->GetXaxis()->GetBinUpEdge(i_cos+1);
	  double cos = 0.5 * (cMax+cMin);
	  
	  double data_val = pHistR[i]->GetBinContent(i_cos+1);
	  double data_err = pHistR[i]->GetBinError(i_cos+1);
	  double fit_val = f_fit->Eval(cos);
	  
	  if(cos < cR) {
	    pullR[i]->SetBinContent(i_cos+1, (data_val-fit_val)/data_err);
	    pullR[i]->SetBinError(i_cos+1, 0);
	    devR[i]->SetBinContent(i_cos+1, (data_val-fit_val)/fit_val);
	    devR[i]->SetBinError(i_cos+1, 0);
	  }
	}

	pHistR[i]->SetStats(0);
	pHistR[i]->SetMinimum(0);
	pHistR[i]->SetMaximum(par_R[0][i]*2.5);
	pHistR[i]->SetMarkerColor(kBlack);
	pHistR[i]->SetLineColor(kBlack);
	pHistR[i]->GetXaxis()->SetTitle("|cos#theta|");
	pHistR[i]->Draw();
	f_fit->SetLineColor(kRed);
	f_fit->Draw("same");
	
	TLatex lcR;
	lcR.SetTextSize(0.03);
	lcR.DrawLatex(0.1, par_R[0][i]*0.3, Form("#chi^{2}/ndf = %.0f/%.0f", chi2[1][i], ndf[1][i]));
	lcR.DrawLatex(0.1, par_R[0][i]*0.2, Form("P(#chi^{2},ndf) = %.1f%%", 100*chiP[1][i]));
	
	c->SaveAs(Form("plots/bkgSub/RSB_%d_%d.pdf", i_fit+1, i+1));
	c->Clear();

	pullL[i]->SetMinimum(-3);
	pullL[i]->SetMaximum(3);
	pullL[i]->SetStats(0);
	pullL[i]->SetLineColor(kBlue);
	pullL[i]->SetMarkerColor(kBlue);
	pullL[i]->SetMarkerStyle(20);
	pullL[i]->SetMarkerSize(0.75);
	pullL[i]->GetXaxis()->SetTitle("|cos#theta|");
	pullL[i]->GetYaxis()->SetTitle("pulls");
	pullL[i]->Draw("pl");

	pullR[i]->SetLineColor(kBlack);
	pullR[i]->SetMarkerColor(kBlack);
	pullR[i]->SetMarkerStyle(20);
	pullR[i]->SetMarkerSize(0.75);
	pullR[i]->Draw("pl same");

	TLine *clim = new TLine(cR, -3, cR, 3);
	clim->SetLineColor(kRed);
	clim->SetLineStyle(kDashed);
	clim->Draw();
	
	c->SaveAs(Form("plots/bkgSub/pulls_%d_%d.pdf", i_fit+1, i+1));
	c->Clear();
    
	devL[i]->SetMinimum(-0.5);
	devL[i]->SetMaximum(0.5);
	devL[i]->SetStats(0);
	devL[i]->SetLineColor(kBlue);
	devL[i]->SetMarkerColor(kBlue);
	devL[i]->SetMarkerStyle(20);
	devL[i]->SetMarkerSize(0.75);
	devL[i]->GetXaxis()->SetTitle("|cos#theta|");
	devL[i]->GetYaxis()->SetTitle("Fit deviation");
	devL[i]->Draw("pl");
	
	devR[i]->SetLineColor(kBlack);
	devR[i]->SetMarkerColor(kBlack);
	devR[i]->SetMarkerStyle(20);
	devR[i]->SetMarkerSize(0.75);
	devR[i]->Draw("pl same");
	
	TLine *climd = new TLine(cR, -0.5, cR, 0.5);
	climd->SetLineColor(kRed);
	climd->SetLineStyle(kDashed);
	climd->Draw();
	
	c->SaveAs(Form("plots/bkgSub/dev_%d_%d.pdf", i_fit+1, i+1));
	c->Clear();
      }

      // making TGraphs for all the fit parameters
      TGraphErrors *graph_NL  = new TGraphErrors(nBinsY, pt, par_L[0], ept, epar_L[0]);
      TGraphErrors *graph_NR  = new TGraphErrors(nBinsY, pt, par_R[0], ept, epar_R[0]);
      TGraphErrors *graph_l2L = new TGraphErrors(nBinsY, pt, par_L[1], ept, epar_L[1]);
      TGraphErrors *graph_l2R = new TGraphErrors(nBinsY, pt, par_R[1], ept, epar_R[1]);
      TGraphErrors *graph_l4L = new TGraphErrors(nBinsY, pt, par_L[2], ept, epar_L[2]);
      TGraphErrors *graph_l4R = new TGraphErrors(nBinsY, pt, par_R[2], ept, epar_R[2]);
  
      TGraph *graph_CL  = new TGraph(nBinsY, pt, chi2[0]);
      TGraph *graph_NDL = new TGraph(nBinsY, pt, ndf[0]);
      TGraph *graph_PL  = new TGraph(nBinsY, pt, chiP[0]);
      TGraph *graph_CR  = new TGraph(nBinsY, pt, chi2[1]);
      TGraph *graph_NDR = new TGraph(nBinsY, pt, ndf[1]);
      TGraph *graph_PR  = new TGraph(nBinsY, pt, chiP[1]);
      TGraph *graph_CM  = new TGraph(nBinsY, pt, cMaxVal);
      
      c->SetLogy(0);
      TH1F *fc = c->DrawFrame(yBins[0]-5, 0, yBins[nBinsY]+5, 1);
      fc->SetXTitle("p_{T} (GeV)");
      fc->SetYTitle("P(#chi^{2}, ndf)");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->GetYaxis()->SetLabelOffset(0.01);
      fc->SetTitle("2018 P(#chi^{2}, ndf)");
      
      graph_PL->SetLineColor(kBlue);
      graph_PL->SetMarkerColor(kBlue);
      graph_PL->SetMarkerStyle(20);
      graph_PL->SetMarkerSize(1.5);
      graph_PL->Draw("p");
      
      graph_PR->SetLineColor(kBlack);
      graph_PR->SetMarkerColor(kBlack);
      graph_PR->SetMarkerStyle(20);
      graph_PR->SetMarkerSize(1.5);
      graph_PR->Draw("p");
      
      c->SaveAs(Form("plots/bkgSub/fit_chiP_%d.pdf", i_fit+1));
      c->Clear();
      
      c->SetLogy(0);
      TH1F *fN = c->DrawFrame(yBins[0]-5, 0.005, yBins[nBinsY]+5, 0.065);
      fN->SetXTitle("p_{T} (GeV)");
      fN->SetYTitle("N");
      fN->GetYaxis()->SetTitleOffset(1.3);
      fN->GetYaxis()->SetLabelOffset(0.01);
      fN->SetTitle("2018 N");
      
      graph_NL->SetLineColor(kBlue);
      graph_NL->SetMarkerColor(kBlue);
      graph_NL->SetMarkerStyle(20);
      graph_NL->SetMarkerSize(0.75);
      graph_NL->Draw("p");
      
      graph_NR->SetLineColor(kBlack);
      graph_NR->SetMarkerColor(kBlack);
      graph_NR->SetMarkerStyle(20);
      graph_NR->SetMarkerSize(0.75);
      graph_NR->Draw("p");
      
      c->SaveAs(Form("plots/bkgSub/fit_N_%d.pdf", i_fit+1));
      c->Clear();
      
      TH1F *fl2 = c->DrawFrame(yBins[0]-5, -4, yBins[nBinsY]+5, 4);
      fl2->SetXTitle("p_{T} (GeV)");
      fl2->SetYTitle("#lambda_{2}");
      fl2->GetYaxis()->SetTitleOffset(1.3);
      fl2->GetYaxis()->SetLabelOffset(0.01);
      fl2->SetTitle("2018 #lambda_{2}");
      
      graph_l2L->SetLineColor(kBlue);
      graph_l2L->SetMarkerColor(kBlue);
      graph_l2L->SetMarkerStyle(20);
      graph_l2L->SetMarkerSize(0.75);
      graph_l2L->Draw("p");
      
      graph_l2R->SetLineColor(kBlack);
      graph_l2R->SetMarkerColor(kBlack);
      graph_l2R->SetMarkerStyle(20);
      graph_l2R->SetMarkerSize(0.75);
      graph_l2R->Draw("p");
      
      TLine *zero = new TLine(yBins[0]-5, 0, yBins[nBinsY]+5, 0);
      zero->SetLineStyle(kDashed);
      zero->SetLineColor(kBlack);
      zero->Draw();
      
      c->SaveAs(Form("plots/bkgSub/fit_l2_%d.pdf", i_fit+1));
      c->Clear();
      
      TH1F *fl4 = c->DrawFrame(yBins[0]-5, -4, yBins[nBinsY]+5, 4);
      fl4->SetXTitle("p_{T} (GeV)");
      fl4->SetYTitle("#lambda_{4}");
      fl4->GetYaxis()->SetTitleOffset(1.3);
      fl4->GetYaxis()->SetLabelOffset(0.01);
      fl4->SetTitle("2018 #lambda_{4}");
      
      graph_l4L->SetLineColor(kBlue);
      graph_l4L->SetMarkerColor(kBlue);
      graph_l4L->SetMarkerStyle(20);
      graph_l4L->SetMarkerSize(0.75);
      graph_l4L->Draw("p");
      
      graph_l4R->SetLineColor(kBlack);
      graph_l4R->SetMarkerColor(kBlack);
      graph_l4R->SetMarkerStyle(20);
      graph_l4R->SetMarkerSize(0.75);
      graph_l4R->Draw("p");
      
      zero->Draw();
      
      c->SaveAs(Form("plots/bkgSub/fit_l4_%d.pdf", i_fit+1));
      c->Clear();
      
      // storing the fit results in tex format
      ofstream fout;
      fout.open(Form("text_output/mbkg_cos_%d.tex", i_fit+1));
      fout << "\\begin{tabular}{c||c|c|c|c||c|c|c|c}\n";
      fout << "$\\pt$ (GeV) & $N_L$ (1e-3) & $\\lambda_{2,L}$ & $\\lambda_{4,L}$  & $\\chi^2_L$/ndf & $N_R$ (1e-3) & $\\lambda_{2,R}$ & $\\lambda_{4,R}$  & $\\chi^2_R$/ndf \\\\\n";
      fout << "\\hline\n";
      for(int i = 0; i < nBinsY; i++) {
	double pMin = bL_R->GetYaxis()->GetBinLowEdge(i+1);
	double pMax = bL_R->GetYaxis()->GetBinUpEdge(i+1);
	fout << Form("$[%.0f, %.0f]$", pMin, pMax) << " & ";

	fout << Form("%.2f", par_L[0][i]*1e3) << "$\\pm$" << Form("%.2f", epar_L[0][i]*1e3) << " & ";
	if(i_fit <= 1)
	  fout << Form("%.3f", par_L[1][i]) << "$\\pm$" << Form("%.3f", epar_L[1][i]) << " & ";
	else
	  fout << Form("%.1f", par_L[1][i]) << " & ";
	if(i_fit==0)
	  fout << Form("%.2f", par_L[2][i]) << "$\\pm$" << Form("%.2f", epar_L[2][i]) << " & ";
	else
	  fout << Form("%.1f", par_L[2][i]) << " & ";
	fout << Form("%.0f", chi2[0][i]) << "/" << Form("%.0f", ndf[0][i]) << " & ";
	
	fout << Form("%.2f", par_R[0][i]*1e3) << "$\\pm$" << Form("%.2f", epar_R[0][i]*1e3) << " & ";
	if(i_fit <= 1)
	  fout << Form("%.3f", par_R[1][i]) << "$\\pm$" << Form("%.3f", epar_R[1][i]) << " & ";
	else
	  fout << Form("%.1f", par_R[1][i]) << " & ";
	if(i_fit == 0)
	  fout << Form("%.2f", par_R[2][i]) << "$\\pm$" << Form("%.2f", epar_R[2][i]) << " & ";
	else
	  fout << Form("%.1f", par_R[2][i]) << " & ";
	fout << Form("%.0f", chi2[1][i]) << "/" << Form("%.0f", ndf[1][i]) << "\\\\\n";
      }
      fout << "\\end{tabular}\n";
      fout.close();

      // storing the TGraphs
      TFile *fOut3 = new TFile("files/mbSub.root", "UPDATE");
      graph_NL->SetName(Form("graph_NL_%d", i_fit+1));
      graph_NL->Write();
      graph_NR->SetName(Form("graph_NR_%d", i_fit+1));
      graph_NR->Write();
      graph_l2L->SetName(Form("graph_l2L_%d", i_fit+1));
      graph_l2L->Write();
      graph_l2R->SetName(Form("graph_l2R_%d", i_fit+1));
      graph_l2R->Write();
      graph_l4L->SetName(Form("graph_l4L_%d", i_fit+1));
      graph_l4L->Write();
      graph_l4R->SetName(Form("graph_l4R_%d", i_fit+1));
      graph_l4R->Write();
      graph_CL->SetName(Form("graph_chiL_%d", i_fit+1));
      graph_CL->Write();
      graph_CR->SetName(Form("graph_chiR_%d", i_fit+1));
      graph_CR->Write();
      graph_NDL->SetName(Form("graph_ndfL_%d", i_fit+1));
      graph_NDL->Write();
      graph_NDR->SetName(Form("graph_ndfR_%d", i_fit+1));
      graph_NDR->Write();
      graph_PL->SetName(Form("graph_chiPL_%d", i_fit+1));
      graph_PL->Write();
      graph_PR->SetName(Form("graph_chiPR_%d", i_fit+1));
      graph_PR->Write();
      graph_CM->SetName(Form("graph_cosMax_%d", i_fit+1));
      graph_CM->Write();
      fOut3->Close();
    }
  c->Destructor();
}
