#import "../../cosMax/imp_jumpF.C"

// macro to fit the sideband/MC distributions
void fitBkgCosth()
{
  // define plotting lims and labels
  string lbl[] = {"LSB", "RSB"};
  double dev_min = 30, pull_min = 3;

  // get the SB/MC 2d maps from the repository
  TH2D **h_cth = new TH2D*[2];
  TFile *fIn = new TFile("../../PR_fit/files/bkgHist.root");
  for(int i = 0; i < 2; i++) {
    h_cth[i] = (TH2D*)fIn->Get(Form("ratioH%d_ab", i));
    h_cth[i]->SetDirectory(0);
  }
  fIn->Close();
  
  // get the binning
  int nBinsX = h_cth[0]->GetNbinsX(), nBinsY = h_cth[0]->GetNbinsY();
  const double *yBins = h_cth[0]->GetYaxis()->GetXbins()->GetArray();
  double minX = h_cth[0]->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_cth[0]->GetXaxis()->GetBinUpEdge(nBinsX);

  // get the fit range from our cosmax(pT)
  ifstream in;
  string dataS;
  in.open("../../cosMax/cosMaxFitRes.txt");
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
  
  // cycle over 2 background types
  for(int i_inp = 0; i_inp < 2; i_inp++) {

    // prepare output
    TFile *fOutR = new TFile(Form("files/%s_fitres.root", lbl[i_inp].c_str()), "recreate");
    fOutR->Close();

    // cycle over all pT bins
    for(int i = 0; i < nBinsY; i++) {
      
      pt[i] = 0.5*(yBins[i+1]+yBins[i]);
      ept[i] = 0.5*(yBins[i+1]-yBins[i]);
      
      // getting the pT bin projections of Data/MC
      pHist[i] = h_cth[i_inp]->ProjectionX(Form("%s_bin%d_1d", lbl[i_inp].c_str(), i+1), i+1, i+1);
      pHist[i]->SetTitle(Form("2017 %s/MC bin %d: [%.0f, %.0f] GeV", lbl[i_inp].c_str(), i+1, yBins[i], yBins[i+1]));

      // getting the max costh value for the fit
      double pMin = h_cth[i_inp]->GetYaxis()->GetBinLowEdge(i+1);
      double pMax = h_cth[i_inp]->GetYaxis()->GetBinUpEdge(i+1);
	
      cMaxVal[i] = jumpF(cosMax->Integral(pMin, pMax)/(pMax-pMin));

      // define the ratio fit function
      TF1 *f_fit = new TF1("f_fit", "[0]*(1+[1]*x*x+[2]*pow(x,4))", 0, cMaxVal[i]);
      f_fit->SetParNames("N", "l_2", "l_4");
      f_fit->SetParameters(pHist[i]->GetBinContent(1)*1.1, 0.1, 0.1);

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
    
      c->SaveAs(Form("plots/%s/fit_pt%d.pdf", lbl[i_inp].c_str(), i));
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
      fp->SetTitle(Form("2017 %s pulls (%.0f < p_{T} < %.0f GeV)", lbl[i_inp].c_str(), yBins[i], yBins[i+1]));

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
	
      c->SaveAs(Form("plots/%s/pulls_pt%d.pdf", lbl[i_inp].c_str(), i));
      c->Clear();

      // plot the deviations
      TH1F *fd = c->DrawFrame(minX, -dev_min, maxX, dev_min);
      fd->SetXTitle("|cos#theta|");
      fd->SetYTitle("relative difference (%)");
      fd->GetYaxis()->SetTitleOffset(1.3);
      fd->GetYaxis()->SetLabelOffset(0.01);
      fd->SetTitle(Form("2017 %s rel. difference (%.0f < p_{T} < %.0f GeV)", lbl[i_inp].c_str(), yBins[i], yBins[i+1]));

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
	
      c->SaveAs(Form("plots/%s/devs_pt%d.pdf", lbl[i_inp].c_str(), i));
      c->Clear();
    }
    
    // making TGraphs for all the fit parameters
    TGraphErrors *graph_N  = new TGraphErrors(nBinsY, pt, par[0], ept, epar[0]);
    TGraphErrors *graph_l2 = new TGraphErrors(nBinsY, pt, par[1], ept, epar[1]);
    TGraphErrors *graph_l4 = new TGraphErrors(nBinsY, pt, par[2], ept, epar[2]);
    
    TGraph *graph_C  = new TGraph(nBinsY, pt, chi2);
    TGraph *graph_ND = new TGraph(nBinsY, pt, ndf);
    TGraph *graph_P  = new TGraph(nBinsY, pt, chiP);     
    
    // storing the fit results in tex format
    ofstream fout;
    fout.open(Form("text_output/cos_%s.tex", lbl[i_inp].c_str()));
    fout << "\\begin{tabular}{c||c|c|c|c}\n";
    fout << "$\\pt$ (GeV) & $N$ & $\\lambda_{2}$ & $\\lambda_4$  & $\\chi^2$/ndf  \\\\\n";
    fout << "\\hline\n";
    for(int i = 0; i < nBinsY; i++) {
      double pMin = h_cth[i_inp]->GetYaxis()->GetBinLowEdge(i+1);
      double pMax = h_cth[i_inp]->GetYaxis()->GetBinUpEdge(i+1);
      fout << Form("$[%.0f, %.0f]$", pMin, pMax);
      for(int ip = 0; ip < 3; ip++){
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
    TFile *fOut = new TFile(Form("files/%s_fitres.root", lbl[i_inp].c_str()), "update");
    graph_N->SetName(Form("fit_N"));
    graph_N->Write();
    graph_l2->SetName(Form("fit_l2"));
    graph_l2->Write();
    graph_l4->SetName(Form("fit_l4"));
    graph_l4->Write();
    graph_C->SetName(Form("fit_chi"));
    graph_C->Write();
    graph_ND->SetName(Form("fit_ndf"));
    graph_ND->Write();
    graph_P->SetName(Form("fit_chiP"));
    graph_P->Write();
    fOut->Close();
    
    cout << endl << "finished fitting " << lbl[i_inp] << endl << endl;
  }
  c->Destructor();
  
}
