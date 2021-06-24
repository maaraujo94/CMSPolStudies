// macro to fit the sideband/MC distributions
void fitBkgCosth()
{
  TH2D **h_cth = new TH2D*[3];
  string lbl[] = {"NP", "LSB", "RSB"};
  string l2n[] = {"NP", "2", "2"};
  double dev_min[] = {30, 30, 30};
  double pull_min[] = {3, 3, 3};

  double N_min[] = {3e-2, 5e-4, 5e-4};
  double N_max[] = {5e0, 1e-1, 1e-1};
  double l2_min[] = {0.5, 1, 1};
  double l4_min[] = {0.5, 2, 2};
  
  // get the MC and SB 2d maps from the repository
  TFile *fIn = new TFile("files/bkgHist.root");
  for(int i = 0; i < 3; i++) {
    h_cth[i] = (TH2D*)fIn->Get(Form("ratioH%d_ab", i));
    h_cth[i]->SetDirectory(0);
  }
  fIn->Close();
  
  // get the binning
  int nBinsX = h_cth[0]->GetNbinsX(), nBinsY = h_cth[0]->GetNbinsY();
  const double *yBins = h_cth[0]->GetYaxis()->GetXbins()->GetArray();
  double minX = h_cth[0]->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_cth[0]->GetXaxis()->GetBinUpEdge(nBinsX);

  // define the ratio fit function
  TF1 *f_fit = new TF1("f_fit", "[0]*(1+[1]*x*x+[2]*pow(x,4))", minX, maxX);
  f_fit->SetParNames("N", "l_2", "l_4");
  
  // get the fit range from our cosmax(pT)
  ifstream in;
  string dataS;
  in.open("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/2018/PR_fit/text_output/cosMaxFitRes.txt");
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
  double par[3][nBinsY], epar[3][nBinsY];
  double chi2[nBinsY], ndf[nBinsY], chiP[nBinsY];
  double cMaxVal[nBinsY], pt[nBinsY], ept[nBinsY];

  // cycle over 3 background types
  for(int i_inp = 0; i_inp < 3; i_inp++) {

    // cycle over all pT bins
    for(int i = 0; i < nBinsY; i++) {

      pt[i] = 0.5*(yBins[i+1]+yBins[i]);
      ept[i] = 0.5*(yBins[i+1]-yBins[i]);

      // getting the pT bin projections of Data/MC
      pHist[i] = h_cth[i_inp]->ProjectionX(Form("%s_bin%d_1d", lbl[i_inp].c_str(), i+1), i+1, i+1);
      pHist[i]->SetTitle(Form("2018 %s/MC bin %d: [%.0f, %.0f] GeV", lbl[i_inp].c_str(), i+1, yBins[i], yBins[i+1]));
    
      // getting the max costh value for the fit
      double pMin = h_cth[i_inp]->GetYaxis()->GetBinLowEdge(i+1);
      double pMax = h_cth[i_inp]->GetYaxis()->GetBinUpEdge(i+1);
    
      cMaxVal[i] = cosMax->Integral(pMin, pMax)/(pMax-pMin);
      double cR = floor(cMaxVal[i]*10.)/10.;
      if(cMaxVal[i]-cR>0.05) cR += 0.05;
      f_fit->SetRange(0, cR);

      // make sure each of the 3 inputs has same initialization
      f_fit->SetParameters(1., 0.1, 0.1);

      // fitting the histos
      f_fit->SetParameter(0, pHist[i]->GetBinContent(1)*1.1);
      if(i_inp == 0) {
	f_fit->FixParameter(2, 0);
	f_fit->FixParameter(1, -0.146);
      }
      else if(i_inp == 1) {
	f_fit->FixParameter(1, 0.3);
	f_fit->FixParameter(2, 0.5);
	//f_fit->ReleaseParameter(2);
	//f_fit->SetParameter(2, 0.1);
      }
      else if(i_inp == 2) {
	f_fit->FixParameter(1, 0.5);
	//f_fit->SetParameter(2, 0.1);
	f_fit->FixParameter(2, 0.9);
      }
      pHist[i]->Fit("f_fit", "R");
      for(int j = 0; j < 3; j++) {
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
    
      c->SaveAs(Form("plots/%s/cth_%d.pdf", lbl[i_inp].c_str(), i+1));
      c->Clear();

      // calculating pulls
      double cv[nBinsX], pv[nBinsX], dv[nBinsX];
      for(int i_cos = 0 ; i_cos < nBinsX; i_cos++) {
	double cMin = h_cth[i_inp]->GetXaxis()->GetBinLowEdge(i_cos+1);
	double cMax = h_cth[i_inp]->GetXaxis()->GetBinUpEdge(i_cos+1);
	cv[i_cos] = 0.5 * (cMax+cMin);
	
	double fitv = f_fit->Eval(cv[i_cos]);
	double datav = pHist[i]->GetBinContent(i_cos+1);
	double datau = pHist[i]->GetBinError(i_cos+1);
	if(cv[i_cos] < cR) {
	  pv[i_cos] = (datav-fitv)/datau;
	  dv[i_cos] = (datav-fitv)/fitv * 100.;
	}
	else {
	  pv[i_cos] = 0;
	  dv[i_cos] = 0;
	}
      }

      // plot the pulls
      TH1F *fp = c->DrawFrame(minX, -pull_min[i_inp], maxX, pull_min[i_inp]);
      fp->SetXTitle("|cos#theta|");
      fp->SetYTitle("pulls");
      fp->GetYaxis()->SetTitleOffset(1.3);
      fp->GetYaxis()->SetLabelOffset(0.01);
      fp->SetTitle(Form("2018 %s pulls (%.0f < p_{T} < %.0f GeV)", lbl[i_inp].c_str(), yBins[i], yBins[i+1]));

      TGraph *g_pull = new TGraph(nBinsX, cv, pv);
      g_pull->SetLineColor(kBlack);
      g_pull->SetMarkerColor(kBlack);
      g_pull->SetMarkerStyle(20);
      g_pull->Draw("p");

      TLine *clim = new TLine(cR, -pull_min[i_inp], cR, pull_min[i_inp]);
      clim->SetLineColor(kRed);
      clim->SetLineStyle(kDashed);
      clim->Draw();

      TLine *zero = new TLine(minX, 0, maxX, 0);
      zero->SetLineColor(kBlack);
      zero->SetLineStyle(kDashed);
      zero->Draw();
	
      c->SaveAs(Form("plots/%s/pulls_%d.pdf", lbl[i_inp].c_str(), i+1));
      c->Clear();

      // plot the deviations
      TH1F *fd = c->DrawFrame(minX, -dev_min[i_inp], maxX, dev_min[i_inp]);
      fd->SetXTitle("|cos#theta|");
      fd->SetYTitle("relative difference (%)");
      fd->GetYaxis()->SetTitleOffset(1.3);
      fd->GetYaxis()->SetLabelOffset(0.01);
      fd->SetTitle(Form("2018 %s rel. difference (%.0f < p_{T} < %.0f GeV)", lbl[i_inp].c_str(), yBins[i], yBins[i+1]));

      TGraph *g_dev = new TGraph(nBinsX, cv, dv);
      g_dev->SetLineColor(kBlack);
      g_dev->SetMarkerColor(kBlack);
      g_dev->SetMarkerStyle(20);
      g_dev->Draw("p");
		
      TLine *climd = new TLine(cR, -dev_min[i_inp], cR, dev_min[i_inp]);
      climd->SetLineColor(kRed);
      climd->SetLineStyle(kDashed);
      climd->Draw();

      zero->Draw();
	
      c->SaveAs(Form("plots/%s/dev_%d.pdf", lbl[i_inp].c_str(), i+1));
      c->Clear();
    }

    // making TGraphs for all the fit parameters
    TGraphErrors *graph_N  = new TGraphErrors(nBinsY, pt, par[0], ept, epar[0]);
    TGraphErrors *graph_l2 = new TGraphErrors(nBinsY, pt, par[1], ept, epar[1]);
    TGraphErrors *graph_l4 = new TGraphErrors(nBinsY, pt, par[2], ept, epar[2]);
  
    TGraph *graph_C  = new TGraph(nBinsY, pt, chi2);
    TGraph *graph_ND = new TGraph(nBinsY, pt, ndf);
    TGraph *graph_P  = new TGraph(nBinsY, pt, chiP);
    TGraph *graph_CM  = new TGraph(nBinsY, pt, cMaxVal);
      
    c->SetLogy(0);
    TH1F *fc = c->DrawFrame(yBins[0]-5, 0, yBins[nBinsY]+5, 1);
    fc->SetXTitle("p_{T} (GeV)");
    fc->SetYTitle("P(#chi^{2}, ndf)");
    fc->GetYaxis()->SetTitleOffset(1.3);
    fc->GetYaxis()->SetLabelOffset(0.01);
    fc->SetTitle(Form("2018 %s P(#chi^{2}, ndf)", lbl[i_inp].c_str()));
      
    graph_P->SetLineColor(kBlack);
    graph_P->SetMarkerColor(kBlack);
    graph_P->SetMarkerStyle(20);
    graph_P->SetMarkerSize(1.5);
    graph_P->Draw("p");
      
    c->SaveAs(Form("plots/%s/par_chiP.pdf", lbl[i_inp].c_str()));
    c->Clear();
      
    c->SetLogy();
    TH1F *fN = c->DrawFrame(yBins[0]-5, N_min[i_inp], yBins[nBinsY]+5, N_max[i_inp]);
    fN->SetXTitle("p_{T} (GeV)");
    fN->SetYTitle("N");
    fN->GetYaxis()->SetTitleOffset(1.3);
    fN->GetYaxis()->SetLabelOffset(0.01);
    fN->SetTitle(Form("2018 %s N", lbl[i_inp].c_str()));
      
    graph_N->SetLineColor(kBlack);
    graph_N->SetMarkerColor(kBlack);
    graph_N->SetMarkerStyle(20);
    graph_N->SetMarkerSize(0.75);
    graph_N->Draw("p");
      
    c->SaveAs(Form("plots/%s/par_N.pdf", lbl[i_inp].c_str()));
    c->Clear();
      
    c->SetLogy(0);
    TH1F *fl2 = c->DrawFrame(yBins[0]-5, -l2_min[i_inp], yBins[nBinsY]+5, l2_min[i_inp]);
    fl2->SetXTitle("p_{T} (GeV)");
    fl2->SetYTitle(Form("#lambda_{%s}", l2n[i_inp].c_str()));
    fl2->GetYaxis()->SetTitleOffset(1.3);
    fl2->GetYaxis()->SetLabelOffset(0.01);
    fl2->SetTitle(Form("2018 %s #lambda_{%s}", lbl[i_inp].c_str(), l2n[i_inp].c_str()));
      
    graph_l2->SetLineColor(kBlack);
    graph_l2->SetMarkerColor(kBlack);
    graph_l2->SetMarkerStyle(20);
    graph_l2->SetMarkerSize(0.75);
    graph_l2->Draw("p");
            
    TLine *zero = new TLine(yBins[0]-5, 0, yBins[nBinsY]+5, 0);
    zero->SetLineStyle(kDashed);
    zero->SetLineColor(kBlack);
    zero->Draw();

    //  TF1 *fcon = new TF1("fcon", "[0]", yBins[0], yBins[nBinsY]);
    //   graph_l2->Fit(fcon);
   
    
    c->SaveAs(Form("plots/%s/par_l2.pdf", lbl[i_inp].c_str()));
    c->Clear();
      
    TH1F *fl4 = c->DrawFrame(yBins[0]-5, -l4_min[i_inp], yBins[nBinsY]+5, l4_min[i_inp]);
    fl4->SetXTitle("p_{T} (GeV)");
    fl4->SetYTitle("#lambda_{4}");
    fl4->GetYaxis()->SetTitleOffset(1.3);
    fl4->GetYaxis()->SetLabelOffset(0.01);
    fl4->SetTitle(Form("2018 %s #lambda_{4}", lbl[i_inp].c_str()));
      
    graph_l4->SetLineColor(kBlack);
    graph_l4->SetMarkerColor(kBlack);
    graph_l4->SetMarkerStyle(20);
    graph_l4->SetMarkerSize(0.75);
    graph_l4->Draw("p");
      
    zero->Draw();

    //TF1 *fcon = new TF1("fcon", "[0]", yBins[0], yBins[nBinsY]);
    // graph_l4->Fit(fcon);
    
    c->SaveAs(Form("plots/%s/par_l4.pdf", lbl[i_inp].c_str()));
    c->Clear();
      
    // storing the fit results in tex format
    ofstream fout;
    fout.open(Form("text_output/cos_%s.tex", lbl[i_inp].c_str()));
    if(i_inp == 0) {
      fout << "\\begin{tabular}{c||c|c|c}\n";
      fout << "$\\pt$ (GeV) & $N$ & $\\lambda_{NP}$ & $\\chi^2$/ndf  \\\\\n";
    }
    else {
      fout << "\\begin{tabular}{c||c|c|c|c}\n";
      fout << "$\\pt$ (GeV) & $N$ & $\\lambda_{2}$ & $\\lambda_4$  & $\\chi^2$/ndf  \\\\\n";
    }
    fout << "\\hline\n";
    for(int i = 0; i < nBinsY; i++) {
      double pMin = h_cth[i_inp]->GetYaxis()->GetBinLowEdge(i+1);
      double pMax = h_cth[i_inp]->GetYaxis()->GetBinUpEdge(i+1);
      fout << Form("$[%.0f, %.0f]$", pMin, pMax);
      int ilim = (i_inp == 0 ? 2 : 3);
      for(int ip = 0; ip < ilim; ip++){
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
    TFile *fOut = new TFile(Form("files/%s_fitres.root", lbl[i_inp].c_str()), "recreate");
    graph_N->SetName(Form("graph_N"));
    graph_N->Write();
    graph_l2->SetName(Form("graph_l2"));
    graph_l2->Write();
    graph_l4->SetName(Form("graph_l4"));
    graph_l4->Write();
    graph_C->SetName(Form("graph_chi"));
    graph_C->Write();
    graph_ND->SetName(Form("graph_ndf"));
    graph_ND->Write();
    graph_P->SetName(Form("graph_chiP"));
    graph_P->Write();
    graph_CM->SetName(Form("graph_cosMax"));
    graph_CM->Write();
    fOut->Close();

    cout << endl << "finished fitting " << lbl[i_inp] << endl << endl;
  }
  c->Destructor();
}
