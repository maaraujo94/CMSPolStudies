// code to do the minimal |costh|<0.3 fit
// then fill the empty histo bins with extrapolated results
// and fit again with full |costh|<0.8 range
// plots histos at 0.3 and 0.8 and compares with fit functions

// main
void extCosFit()
{
  TCanvas *c = new TCanvas("", "", 700, 700);
  
  // read the coarse 2d histo in |costh|
  TFile *infile = new TFile("files/store_hist.root");
  TH2D *bHist = new TH2D();
  infile->GetObject("cHist_ab", bHist);
  bHist->SetDirectory(0);
  infile->Close();

  TH2D* hist = (TH2D*)bHist->Clone();
  
  // cout some useful info for checking purposes
  int nBinsX = hist->GetNbinsX(), nBinsY = hist->GetNbinsY();
  double minX = hist->GetXaxis()->GetBinLowEdge(1), maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);
  const double *yBins = hist->GetYaxis()->GetXbins()->GetArray();

  cout << "fitting histogram " << hist->GetName() << endl;
  cout << "X axis: " << nBinsX << " bins in [" << minX << "," << maxX << "]" << endl;
  cout << "Y axis: " << nBinsY << " bins in [" << hist->GetYaxis()->GetBinLowEdge(1) << "," << hist->GetYaxis()->GetBinUpEdge(nBinsY) << "]" << endl;
  
  // the function to be used - linear lth = m(y-y_ref)+lth_ref
  TF2 *fitMin = new TF2("fit m 2d", "[0]*(1+([1]*(y-[3])+[2])*x*x)", 0., 0.3, 12., 70.);
  fitMin->SetParameters(20, 0.1, 0.1, 20);
  fitMin->SetParNames("A", "m", "lth_ref", "pT_ref");
  fitMin->FixParameter(3, 20);
  
  bHist->Fit(fitMin, "R");

  // save fit output
  ofstream outtex;
  outtex.open("text_output/minMaxFit.tex");
  outtex << "\\begin{tabular}{c|c|c|c|c|c}\n";
  outtex << Form("Fit & max $|\\cost|$ & $A$ & $m$ & $\\lambda_\\theta(\\pt=%.0f$ GeV) & $\\chi^2$/ndf \\\\\n", fitMin->GetParameter(3));
  outtex << "\\hline\n";
  outtex << "Minimal & $0.3$ & ";
  for(int i = 0; i < 3; i++) {
    int prec = ceil(-log10(fitMin->GetParError(i)))+1;
    outtex << "$" << setprecision(prec) << fixed << fitMin->GetParameter(i) << "\\pm" << fitMin->GetParError(i) << "$ & ";
  }
  outtex << setprecision(0) << fixed <<  fitMin->GetChisquare() << "/" << fitMin->GetNDF() <<  " \\\\\n";
  outtex.close();

  // after minimal fit we want to extrapolate the histogram values
  // first read the |costh|max(pt) function
  ifstream in;
  string dataS;
  in.open("text_output/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", 0, 100);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);

  double binErr = 100.*hist->GetBinError(1, nBinsY);
  // for each bin determine the cos limit and fill with fitMin outside
  for (int i = 0; i < nBinsY; i++)
    {
      double pMin = hist->GetYaxis()->GetBinLowEdge(i+1);
      double pMax = hist->GetYaxis()->GetBinUpEdge(i+1);
      double cMaxVal = cosMax->Integral(pMin, pMax)/(pMax-pMin);

      for(int j = 0; j < nBinsX; j++) {
	if(hist->GetXaxis()->GetBinUpEdge(j+1) > cMaxVal)
	  {
	    double cMin = hist->GetXaxis()->GetBinLowEdge(j+1);
	    double cMax = hist->GetXaxis()->GetBinUpEdge(j+1);
	    double fillVal = fitMin->Integral(cMin, cMax, pMin, pMax)/((cMax-cMin)*(pMax-pMin));
	    
	    hist->SetBinContent(j+1, i+1, fillVal);
	    hist->SetBinError(j+1, i+1, binErr);
	  }    
	
      }
    }

  // save the histogram of the extrapolated info
  dataS = hist->GetName();
  hist->SetName(Form("%s_ext", dataS.c_str()));
  dataS = hist->GetTitle();
  hist->SetTitle(Form("%s extrapolated", dataS.c_str()));

  // redo the fit using the full width of the histogram
  TF2 *fitCom = new TF2("fit c 2d", "[0]*(1+([1]*(y-[3])+[2])*x*x)", 0., 0.8, 12., 70.);
  fitCom->SetParameters(20, 0.1, 0.1, 20);
  fitCom->SetParNames("A", "m", "lth_ref", "pT_ref");
  fitCom->FixParameter(3, 20);
  
  hist->Fit(fitCom, "R");
  
  // plot the histogram of the extrapolated info
  hist->SetStats(0);
  hist->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  hist->GetYaxis()->SetTitle("p_{T} (GeV)");
  hist->GetZaxis()->SetRangeUser(0,50);
  hist->Draw("hist COLZ");

  TF1 *invF = new TF1("invF", "(exp(x/[0])-[1])/[2]", 0, 0.9);
  invF->SetParameters(maxPar[0], maxPar[1], maxPar[2]);
  invF->SetLineColor(kRed);
  invF->Draw("same");
  
  c->SaveAs("plots/ratio_2d_abs_ext.pdf");
  c->Clear();

  // plot the original histogram
  bHist->SetStats(0);
  bHist->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  bHist->GetYaxis()->SetTitle("p_{T} (GeV)");
  bHist->GetZaxis()->SetRangeUser(0,50);
  bHist->Draw("hist COLZ");

  invF->Draw("same");
  
  c->SaveAs("plots/ratio_2d_abs_coarse.pdf");
  c->Clear();

  // save fit results to tex
  outtex.open("text_output/minMaxFit.tex", fstream::app);
  outtex << "Complete & $0.8$ & ";
  for(int i = 0; i < 3; i++) {
    int prec = ceil(-log10(fitCom->GetParError(i)))+1;
    outtex << "$" << setprecision(prec) << fixed << fitCom->GetParameter(i) << "\\pm" << fitCom->GetParError(i) << "$ & ";
  }
  outtex << setprecision(0) << fixed <<  fitCom->GetChisquare() << "/" << fitCom->GetNDF() << endl;
  outtex << "\\end{tabular}\n";
  outtex.close();

  // save plots to .root file
  TFile* outfile = new TFile("files/store_hist.root", "update");
  hist->Write(0, TObject::kOverwrite);
  outfile->Close();

  // now the illustrative plots
  // first the histogram plots in the given range compared to function
  TH2D *minHist = new TH2D("minHist", "Minimal fit function", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *comHist = new TH2D("comHist", "Complete fit function", nBinsX, minX, maxX, nBinsY, yBins);
  for (int i = 0; i < nBinsY; i++) {
    double pMin = hist->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = hist->GetYaxis()->GetBinUpEdge(i+1);
    for(int j = 0; j < nBinsX; j++) {
      double cMin = hist->GetXaxis()->GetBinLowEdge(j+1);
      double cMax = hist->GetXaxis()->GetBinUpEdge(j+1);

      double fillVal = fitCom->Integral(cMin, cMax, pMin, pMax)/((cMax-cMin)*(pMax-pMin));
      comHist->SetBinContent(j+1, i+1, fillVal);
	
      fillVal = fitMin->Integral(cMin, cMax, pMin, pMax)/((cMax-cMin)*(pMax-pMin));
      minHist->SetBinContent(j+1, i+1, fillVal);

    }
  }

  // plot the minimal plots
  TLine *minLine = new TLine(0.3, 12, 0.3, 70);
  minLine->SetLineColor(kRed);
  
  bHist->SetStats(0);
  bHist->GetZaxis()->SetRangeUser(0,50); 
  bHist->GetListOfFunctions()->Remove(bHist->GetFunction("fit m 2d"));
  bHist->Draw("COLZ");
  minLine->Draw("same");
  c->SaveAs("plots/fit_min_ratio.pdf");
  c->Clear();

  minHist->SetStats(0);
  minHist->GetZaxis()->SetRangeUser(0,50);
  minHist->Draw("hist COLZ");
  minLine->Draw("same");
  c->SaveAs("plots/fit_min_fit.pdf");
  c->Clear();

  // plot the complete plots
  TLine *comLine = new TLine(0.8, 12, 0.8, 70);
  comLine->SetLineColor(kRed);

  hist->SetStats(0);
  hist->GetZaxis()->SetRangeUser(0,50);
  hist->GetListOfFunctions()->Remove(hist->GetFunction("fit c 2d"));
  hist->Draw("COLZ");
  comLine->Draw("same");
  c->SaveAs("plots/fit_ext_ratio.pdf");
  c->Clear();
  
  comHist->SetStats(0);
  comHist->GetZaxis()->SetRangeUser(0,50);
  comHist->Draw("COLZ");
  comLine->Draw("same");
  c->SaveAs("plots/fit_ext_fit.pdf");
  c->Clear();

  // plot the pulls: (data_val-fit_val)/data_err
  TH2D *minPull = new TH2D("minPull", "Minimal fit pulls", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *comPull = new TH2D("comPull", "Complete fit pulls", nBinsX, minX, maxX, nBinsY, yBins);

  for (int i = 0; i < nBinsY; i++) {
    double pMin = hist->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = hist->GetYaxis()->GetBinUpEdge(i+1);
    for(int j = 0; j < nBinsX; j++) {
      double cMin = hist->GetXaxis()->GetBinLowEdge(j+1);
      double cMax = hist->GetXaxis()->GetBinUpEdge(j+1);

      if(cMin < 0.8) {
	double data_val = hist->GetBinContent(j+1, i+1);
	double data_err = hist->GetBinError(j+1, i+1);
	double fit_val = fitCom->Integral(cMin, cMax, pMin, pMax)/((cMax-cMin)*(pMax-pMin));
	comPull->SetBinContent(j+1, i+1, (data_val-fit_val)/data_err);

	if(cMin < 0.3) {
	  double data_val = bHist->GetBinContent(j+1, i+1);
	  double data_err = bHist->GetBinError(j+1, i+1);
	  double fit_val = fitMin->Integral(cMin, cMax, pMin, pMax)/((cMax-cMin)*(pMax-pMin));
	  minPull->SetBinContent(j+1, i+1, (data_val-fit_val)/data_err);
	}
      }
    }
  }

  // saving the pulls plots  
  minPull->SetStats(0);
  minPull->GetZaxis()->SetRangeUser(-5,5);
  minPull->Draw("COLZ");
  minLine->Draw("same");
  c->SaveAs("plots/fit_min_pull.pdf");
  c->Clear();

  comPull->SetStats(0);
  comPull->GetZaxis()->SetRangeUser(-5,5);
  comPull->Draw("COLZ");
  comLine->Draw("same");
  c->SaveAs("plots/fit_ext_pull.pdf");
  c->Clear();

  // plot the ratio data/func
  TH2D *minDiv = new TH2D("minDiv", "Minimal fit data/fit", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *comDiv = new TH2D("comDiv", "Complete fit data/fit", nBinsX, minX, maxX, nBinsY, yBins);

  minDiv = (TH2D*)bHist->Clone(Form("minDiv"));
  minDiv->Divide(minHist);
  minDiv->SetMaximum(2.1);
  minDiv->SetTitle("Minimal fit data/fit");
  minDiv->Draw("COLZ");
  minLine->Draw();
  c->SaveAs("plots/fit_min_quot.pdf");
  c->Clear();

  comDiv = (TH2D*)hist->Clone(Form("comDiv"));
  comDiv->Divide(comHist);
  comDiv->SetMaximum(2.1);
  comDiv->SetTitle("Complete fit data/fit");
  comDiv->Draw("COLZ");
  comLine->Draw();
  c->SaveAs("plots/fit_ext_quot.pdf");
  c->Clear();

  // plot the 1D coarse bins with fit function superimposed (integ over pt bin)
  // the projections before and after extrapolation
  TH1D *pMinHist[nBinsY], *pComHist[nBinsY];
  TH1D *fpMinHist[nBinsY], *fpComHist[nBinsY];
  double pjMax[] = {20, 20, 25, 30, 35, 40, 50, 55};
  for(int i = 1; i <= nBinsY; i++) {
    pMinHist[i-1] = bHist->ProjectionX(Form("coarse_bin%d_1d_min", i), i, i);
    pMinHist[i-1]->SetTitle(Form("p_{T} bin %d: [%.0f, %.0f] GeV", i, yBins[i-1], yBins[i]));

    pComHist[i-1] = hist->ProjectionX(Form("coarse_bin%d_1d_com", i), i, i);
    pComHist[i-1]->SetTitle(Form("p_{T} bin %d: [%.0f, %.0f] GeV", i, yBins[i-1], yBins[i]));

    fpMinHist[i-1] = minHist->ProjectionX(Form("coarse_bin%d_1d_min_f", i), i, i);
    fpMinHist[i-1]->SetTitle(Form("Minimal fit c bin %d: [%.0f, %.0f] GeV", i, yBins[i-1], yBins[i]));

    fpComHist[i-1] = comHist->ProjectionX(Form("coarse_bin%d_1d_com_f", i), i, i);
    fpComHist[i-1]->SetTitle(Form("Complete fit c bin %d: [%.0f, %.0f] GeV", i, yBins[i-1], yBins[i]));

    pComHist[i-1]->SetStats(0);
    pComHist[i-1]->SetLineColor(kBlack);
    pMinHist[i-1]->SetLineColor(kBlue);
    pMinHist[i-1]->SetStats(0);
    fpComHist[i-1]->SetLineColor(kRed);
    fpMinHist[i-1]->SetLineColor(kViolet);
    
    pComHist[i-1]->GetYaxis()->SetRangeUser(0, pjMax[i-1]);
    pComHist[i-1]->Draw();
    pMinHist[i-1]->Draw("same");
    fpComHist[i-1]->Draw("same");
    fpMinHist[i-1]->Draw("same");
    
    c->SaveAs(Form("plots/proj_pt%d.pdf", i));

    // draw limits of each fit
    TLine *min1D = new TLine(0.3, c->GetUymin(), 0.3, c->GetUymax());
    min1D->SetLineStyle(kDashed);
    min1D->SetLineColor(kBlue);
    min1D->Draw();
    TLine *com1D = new TLine(0.8, c->GetUymin(), 0.8, c->GetUymax());
    com1D->SetLineStyle(kDashed);
    com1D->SetLineColor(kBlack);
    com1D->Draw();

    c->SaveAs(Form("plots/proj_pt%d.pdf", i));
    c->Clear();
  }
  
  c->Destructor();
}

