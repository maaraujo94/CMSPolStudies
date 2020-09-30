// code to do the minimal |costh|<0.3 fit
// then fill the empty histo bins with extrapolated results
// and fit again with full |costh|<0.8 range
// plots histos at 0.3 and 0.8 and compares with fit functions

TH2D* hist;

// the fit function with all the free normalization parameters
double fit_func(double *xx, double *par)
{
  // attribute x vars
  double cos = xx[0], pt = xx[1];

  // attribute parameters
  int NA = hist->GetNbinsY();
  double A[NA];
  for(int i = 0; i < NA; i++) A[i] = par[i];
  double m = par[NA], lth_ref = par[NA+1], pt_ref = par[NA+2];

  int i;
  double N = 0;
  for(i = 1; i <= NA; i++) 
    if(pt >= hist->GetYaxis()->GetBinLowEdge(i) && pt < hist->GetYaxis()->GetBinUpEdge(i)) {
      N = A[i-1];
      break;
    }
  if( pt == hist->GetYaxis()->GetBinUpEdge(NA)) N = A[NA-1];
  
  if(N==0) cout << N;
  
  return N*(1+(m*(pt-pt_ref)+lth_ref)*cos*cos);
}

// main
void extCosFit()
{
  double M_q = 1.;//3.097; // using the J/psi mass
  string fit_type = "linear";

  TCanvas *c = new TCanvas("", "", 700, 700);
  
  // read the coarse 2d histo in |costh|
  TFile *infile = new TFile("files/ratioHist.root");
  TH2D *bHist = new TH2D();
  infile->GetObject("ratioHist_ab", bHist);
  bHist->SetDirectory(0);
  infile->Close();

  hist = (TH2D*)bHist->Clone();
  
  // save some useful info
  int nBinsX = hist->GetNbinsX(), nBinsY = hist->GetNbinsY();
  double minX = hist->GetXaxis()->GetBinLowEdge(1), maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);
  const double *yBins = hist->GetYaxis()->GetXbins()->GetArray();

  double pT_ref = 30.;
  double histoMax = 7.7;
  double ratioMax = 3.2;
  double pullMax = 4.0;
  double mult = 1.4;
  //  double pjMax[] = {40, 40, 40, 40, 40, 40, 40, 40};
  
  // the function to be used - linear lth = m(y-y_ref)+lth_ref
  TF2 *fitMin = new TF2("fit m 2d", "fit_func", 0., 0.3, 25./M_q, 200./M_q, nBinsY+3, 2);
  for(int i = 0; i < nBinsY; i++){
    fitMin->SetParameter(i, 15);
    fitMin->SetParName(i, Form("A_%d", i+1));
  }
  fitMin->SetParameter(nBinsY, 0.1);
  fitMin->SetParameter(nBinsY+1, 0.1);
  fitMin->SetParameter(nBinsY+2, pT_ref);
  fitMin->FixParameter(nBinsY+2, pT_ref);
  if(fit_type == "constant") fitMin->FixParameter(nBinsY, 0);
  fitMin->SetParName(nBinsY, "m");
  fitMin->SetParName(nBinsY+1, "lth_ref");
  fitMin->SetParName(nBinsY+2, "xi_ref");
    
  bHist->Fit(fitMin, "R");
  
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
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", 0, 215);
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

  // redo the fit using the full width of the histogram
  TF2 *fitCom = new TF2("fit c 2d", "fit_func", 0., 0.8, 25./M_q, 200./M_q, nBinsY+3, 2);
  for(int i = 0; i < nBinsY; i++){
    fitCom->SetParameter(i, 15);
    fitCom->SetParName(i, Form("A_%d", i+1));
  }
  fitCom->SetParameter(nBinsY, 0.1);
  fitCom->SetParameter(nBinsY+1, 0.1);
  fitCom->SetParameter(nBinsY+2, pT_ref);
  fitCom->FixParameter(nBinsY+2, pT_ref);
  if(fit_type == "constant") fitCom->FixParameter(nBinsY, 0);
  fitCom->SetParName(nBinsY, "m");
  fitCom->SetParName(nBinsY+1, "lth_ref");
  fitCom->SetParName(nBinsY+2, "xi_ref");
  
  TFitResultPtr fitres = hist->Fit(fitCom, "RS");
  
  // plot the histogram of the extrapolated info
  dataS = hist->GetName();
  hist->SetName(Form("%s_ext", dataS.c_str()));
  dataS = hist->GetTitle();
  hist->SetTitle(Form("%s extrapolated (%s)", dataS.c_str(), fit_type.c_str()));

  cout << "histo variation:" << endl;
  cout <<  hist->GetMaximum()  << " / " << bHist->GetMaximum() << endl;
  
  hist->SetStats(0);
  hist->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  hist->GetYaxis()->SetTitle("p_{T}/M");
  hist->GetZaxis()->SetRangeUser(0, histoMax);
  hist->Draw("hist COLZ");

  TF1 *invF = new TF1("invF", "(exp(x/[0])-[1])/[2]", 0, 0.9);
  invF->SetParameters(maxPar[0], maxPar[1], maxPar[2]);
  invF->SetLineColor(kRed);
  invF->Draw("same");
  
  c->SaveAs(Form("plots/fit_%s/ratio_2d_abs_ext.pdf", fit_type.c_str()));
  c->Clear();

  // plot the original histogram
  bHist->SetStats(0);
  bHist->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  bHist->GetYaxis()->SetTitle("p_{T}/M");
  bHist->GetZaxis()->SetRangeUser(0,histoMax);
  bHist->Draw("hist COLZ");

  invF->Draw("same");
  
  c->SaveAs(Form("plots/fit_%s/ratio_2d_abs_coarse.pdf", fit_type.c_str()));
  c->Clear();

  // save fit results to tex
  ofstream outtex;
  outtex.open(Form("text_output/fit_%s/minMaxFit.tex", fit_type.c_str()));
  outtex << "\\begin{tabular}{c|c|c}\n";
  outtex << "$|\\cost|_{max}$ & 0.3 & 0.8\\\\\n";
  outtex << "\\hline\n";
  
  for(int i=0; i<nBinsY; i++) {
    int prec_min = ceil(-log10(abs(fitMin->GetParError(i))))+1;
    prec_min = max(prec_min, 0);
    int prec_com = ceil(-log10(abs(fitCom->GetParError(i))))+1;
    prec_com = max(prec_com, 0);    
    outtex << Form("$A_%d$ & $", i+1) << setprecision(prec_min) << fixed << fitMin->GetParameter(i) << "\\pm" << fitMin->GetParError(i) << "$ & $" << setprecision(prec_com) << fixed << fitCom->GetParameter(i) << "\\pm" << fitCom->GetParError(i) << "$ \\\\\n";
  }
  
  int prec_min = ceil(-log10(abs(fitMin->GetParError(nBinsY))))+1;
  prec_min = max(prec_min, 0);
  int prec_com = ceil(-log10(abs(fitCom->GetParError(nBinsY))))+1;
  prec_com = max(prec_com, 0);    
  outtex << "$m$ & $" << setprecision(prec_min) << fixed << fitMin->GetParameter(nBinsY) << "\\pm" << fitMin->GetParError(nBinsY) << "$ & $" << setprecision(prec_com) << fixed << fitCom->GetParameter(nBinsY) << "\\pm" << fitCom->GetParError(nBinsY) << "$ \\\\\n";

  prec_min = ceil(-log10(abs(fitMin->GetParError(nBinsY+1))))+1;
  prec_min = max(prec_min, 0);
  prec_com = ceil(-log10(abs(fitCom->GetParError(nBinsY+1))))+1;
  prec_com = max(prec_com, 0);
  outtex << Form("$\\lambda_\\theta(\\pt/M = %.1f$) & $", pT_ref) << setprecision(prec_min) << fixed << fitMin->GetParameter(nBinsY+1) << "\\pm" << fitMin->GetParError(nBinsY+1) << "$ & $" << setprecision(prec_com) << fixed << fitCom->GetParameter(nBinsY+1) << "\\pm" << fitCom->GetParError(nBinsY+1) << "$ \\\\\n";

  outtex << "\\hline\n";
  outtex << "$\\chi^2$/ndf & " << setprecision(0) << fixed <<  fitMin->GetChisquare() << "/" << fitMin->GetNDF() << " & " << fitCom->GetChisquare() << "/" << fitCom->GetNDF() << endl;
  outtex << "\\end{tabular}\n";
  
  outtex.close();
  
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
  TLine *minLine = new TLine(0.3, 25./M_q, 0.3, 200./M_q);
  minLine->SetLineColor(kRed);

  cout << "histo variation:" << endl;
  cout <<  minHist->GetMaximum()  << " / " << comHist->GetMaximum() << endl;
  
  bHist->SetStats(0);
  bHist->GetZaxis()->SetRangeUser(0, histoMax); 
  bHist->GetListOfFunctions()->Remove(bHist->GetFunction("fit m 2d"));
  bHist->SetTitle(Form("%s (%s)", bHist->GetTitle(), fit_type.c_str()));
  bHist->Draw("COLZ");
  minLine->Draw("same");
  c->SaveAs(Form("plots/fit_%s/fit_min_ratio.pdf", fit_type.c_str()));
  c->Clear();
  
  minHist->SetStats(0);
  minHist->GetZaxis()->SetRangeUser(0, histoMax);
  minHist->SetTitle(Form("%s (%s)", minHist->GetTitle(), fit_type.c_str()));
  minHist->Draw("hist COLZ");
  minLine->Draw("same");
  c->SaveAs(Form("plots/fit_%s/fit_min_fit.pdf", fit_type.c_str()));
  c->Clear();

  // plot the complete plots
  TLine *comLine = new TLine(0.8, 25./M_q, 0.8, 200./M_q);
  comLine->SetLineColor(kRed);

  hist->SetStats(0);
  hist->GetZaxis()->SetRangeUser(0, histoMax);
  hist->GetListOfFunctions()->Remove(hist->GetFunction("fit c 2d"));
  hist->Draw("COLZ");
  comLine->Draw("same");
  c->SaveAs(Form("plots/fit_%s/fit_ext_ratio.pdf", fit_type.c_str()));
  c->Clear();
  
  comHist->SetStats(0);
  comHist->GetZaxis()->SetRangeUser(0, histoMax);
  comHist->SetTitle(Form("%s (%s)", comHist->GetTitle(), fit_type.c_str()));
  comHist->Draw("COLZ");
  comLine->Draw("same");
  c->SaveAs(Form("plots/fit_%s/fit_ext_fit.pdf", fit_type.c_str()));
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

  cout << "pulls variation:" << endl;
  cout <<  minPull->GetMinimum() << " - " << minPull->GetMaximum() << endl;
  cout <<  comPull->GetMinimum() << " - " << comPull->GetMaximum() << endl;
  
  // saving the pulls plots  
  minPull->SetStats(0);
  minPull->GetZaxis()->SetRangeUser(-1.*pullMax, pullMax);
  minPull->Draw("COLZ");
  minPull->SetTitle(Form("%s (%s)", minPull->GetTitle(), fit_type.c_str()));
  minLine->Draw("same");
  c->SaveAs(Form("plots/fit_%s/fit_min_pull.pdf", fit_type.c_str()));
  c->Clear();

  comPull->SetStats(0);
  comPull->GetZaxis()->SetRangeUser(-1.*pullMax, pullMax);
  comPull->SetTitle(Form("%s (%s)", comPull->GetTitle(), fit_type.c_str()));
  comPull->Draw("COLZ");
  comLine->Draw("same");
  c->SaveAs(Form("plots/fit_%s/fit_ext_pull.pdf", fit_type.c_str()));
  c->Clear();

  // plot the ratio data/func
  TH2D *minDiv = new TH2D("minDiv", "Minimal fit data/fit", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *comDiv = new TH2D("comDiv", "Complete fit data/fit", nBinsX, minX, maxX, nBinsY, yBins);

  minDiv = (TH2D*)bHist->Clone(Form("minDiv"));
  minDiv->Divide(minHist);
  minDiv->SetMaximum();
  comDiv = (TH2D*)hist->Clone(Form("comDiv"));
  comDiv->Divide(comHist);
  comDiv->SetMaximum();
  
  cout << "ratio variation:" << endl;
  cout <<  minDiv->GetMaximum()  << " / " << comDiv->GetMaximum() << endl;
  
  minDiv->SetTitle("Minimal fit data/fit");
  minDiv->SetMaximum(ratioMax);
  minDiv->SetTitle(Form("%s (%s)", minDiv->GetTitle(), fit_type.c_str()));
  minDiv->Draw("COLZ");
  minLine->Draw();
  c->SaveAs(Form("plots/fit_%s/fit_min_quot.pdf", fit_type.c_str()));
  c->Clear();
  
  comDiv->SetTitle("Complete fit data/fit");
  comDiv->SetMaximum(ratioMax);
  comDiv->SetTitle(Form("%s (%s)", comDiv->GetTitle(), fit_type.c_str()));
  comDiv->Draw("COLZ");
  comLine->Draw();
  c->SaveAs(Form("plots/fit_%s/fit_ext_quot.pdf", fit_type.c_str()));
  c->Clear();

  // plot the 1D coarse bins with fit function superimposed (integ over pt bin)
  // the projections before and after extrapolation
  TH1D *pMinHist[nBinsY], *pComHist[nBinsY];
  TH1D *fpMinHist[nBinsY], *fpComHist[nBinsY];
  for(int i = 1; i <= nBinsY; i++) {
    pMinHist[i-1] = bHist->ProjectionX(Form("coarse_bin%d_1d_min", i), i, i);
    pMinHist[i-1]->SetTitle(Form("p_{T} bin %d: [%.0f, %.0f] GeV", i, yBins[i-1]*M_q, yBins[i]*M_q));

    pComHist[i-1] = hist->ProjectionX(Form("coarse_bin%d_1d_com", i), i, i);
    pComHist[i-1]->SetTitle(Form("p_{T} bin %d: [%.0f, %.0f] GeV", i, yBins[i-1]*M_q, yBins[i]*M_q));

    fpMinHist[i-1] = minHist->ProjectionX(Form("coarse_bin%d_1d_min_f", i), i, i);
    fpMinHist[i-1]->SetTitle(Form("Minimal fit c bin %d: [%.0f, %.0f] GeV", i, yBins[i-1]*M_q, yBins[i]*M_q));

    fpComHist[i-1] = comHist->ProjectionX(Form("coarse_bin%d_1d_com_f", i), i, i);
    fpComHist[i-1]->SetTitle(Form("Complete fit c bin %d: [%.0f, %.0f] GeV", i, yBins[i-1]*M_q, yBins[i]*M_q));

    pComHist[i-1]->SetStats(0);
    pComHist[i-1]->SetLineColor(kBlack);
    pMinHist[i-1]->SetLineColor(kBlue);
    pMinHist[i-1]->SetStats(0);
    fpComHist[i-1]->SetLineColor(kRed);
    fpMinHist[i-1]->SetLineColor(kViolet);
    
    pComHist[i-1]->GetYaxis()->SetRangeUser(0, mult*fitCom->GetParameter(i-1));
    pComHist[i-1]->Draw();
    pMinHist[i-1]->Draw("same");
    fpComHist[i-1]->Draw("same");
    //fpMinHist[i-1]->Draw("same");
    
    c->SaveAs(Form("plots/fit_%s/proj_pt%d.pdf", fit_type.c_str(), i));

    // draw limits of each fit
    TLine *min1D = new TLine(0.3, c->GetUymin(), 0.3, c->GetUymax());
    min1D->SetLineStyle(kDashed);
    min1D->SetLineColor(kBlue);
    min1D->Draw();
    TLine *com1D = new TLine(0.8, c->GetUymin(), 0.8, c->GetUymax());
    com1D->SetLineStyle(kDashed);
    com1D->SetLineColor(kBlack);
    com1D->Draw();

    c->SaveAs(Form("plots/fit_%s/proj_pt%d.pdf", fit_type.c_str(), i));
    c->Clear();
  }

  // plots of the fit quantities A_i and lambda_theta(pT)
  TH1D *histA = new TH1D("histA", "A_{i} (p_{T}/M)", nBinsY, yBins);
  for(int i = 0; i < nBinsY; i++) {
    histA->SetBinContent(i+1, fitCom->GetParameter(i));
    histA->SetBinError(i+1, fitCom->GetParError(i));
  }
  histA->SetStats(0);
  histA->GetXaxis()->SetTitle("p_{T}/M");
  histA->GetYaxis()->SetRangeUser(15, 210);
  histA->SetTitle(Form("%s (%s)", histA->GetTitle(), fit_type.c_str()));
  histA->Draw("error");
  histA->Draw("hist same");
  c->SaveAs(Form("plots/fit_%s/fit_A_vals.pdf", fit_type.c_str()));
  c->Clear();

  // lambda_theta(xi) plotted as uncertainty band
  TF1 *funcL = new TF1("funcL", "[0]*(x-[2])+[1]", 0, 215);
  funcL->SetParameters(fitCom->GetParameter(nBinsY), fitCom->GetParameter(nBinsY+1), fitCom->GetParameter(nBinsY+2));

  const int sizeunc = 201-25;
  float l_band[4][sizeunc];
  double dfm, dfl, ln = 1e4;
  double m = fitCom->GetParameter(nBinsY), errm = fitCom->GetParError(nBinsY);
  double l = fitCom->GetParameter(nBinsY+1), errl = fitCom->GetParError(nBinsY+1);
  double errml = fitres->GetCovarianceMatrix()(nBinsY, nBinsY+1);

  for(int i=0; i<sizeunc; i++) {
    // x- and y-vals
    l_band[0][i] = (25.+i)/M_q; 
    l_band[1][i] = funcL->Eval(l_band[0][i]);

    // x uncertainty
    l_band[2][i]=0.5/M_q;

    // uncert contribution from m
    if(m != 0) {
      funcL->SetParameter(0, m + errm/ln);
      dfm = (funcL->Eval(l_band[0][i])-l_band[1][i]) / (errm/ln);
      funcL->SetParameter(0, m);
    
      // uncert contrib from lth_ref
      funcL->SetParameter(1, l + errl/ln);
      dfl = (funcL->Eval(l_band[0][i])-l_band[1][i]) / (errl/ln);
      funcL->SetParameter(1, l);
      
      // y uncertainty
      l_band[3][i] = sqrt(dfm*dfm*errm*errm + dfl*dfl*errl*errl + 2*dfm*dfl*errml);
    }
    else {
      l_band[3][i] = errl;
    }
  }

  TH1F *func = c->DrawFrame(0, -1, 215, 1);
  func->SetXTitle("p_{T}/M");
  func->SetYTitle("#lambda_{#theta}");
  func->GetYaxis()->SetTitleOffset(1.3);
  func->GetYaxis()->SetLabelOffset(0.01);
  func->SetTitle(Form("#lambda_{#theta} (p_{T}/M) (%s)", fit_type.c_str()));
  
  TGraphErrors* l_unc = new TGraphErrors(sizeunc, l_band[0], l_band[1], l_band[2], l_band[3]);
  l_unc->SetLineColor(kBlack);
  l_unc->SetFillColorAlpha(kBlue, 0.5);
  l_unc->SetTitle(Form("#lambda_{#theta} (p_{T}/M) (%s)", fit_type.c_str()));
  l_unc->Draw("ce3");
  
  TF1 *zero = new TF1("zero", "0", 0, 215);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw("lsame");

  c->SaveAs(Form("plots/fit_%s/fit_lth_pt.pdf", fit_type.c_str()));
  c->Clear();

  c->Destructor();

  // save lambda_theta band to root file
  TFile* outfile = new TFile("files/fit_res_2d.root", "update");
  hist->SetName(Form("%s_%s", hist->GetName(), fit_type.c_str()));
  hist->Write(0, TObject::kOverwrite);
  l_unc->SetName(Form("lth_%s", fit_type.c_str()));
  l_unc->Write(0, TObject::kOverwrite);
  histA->SetName(Form("A_%s", fit_type.c_str()));
  histA->Write(0, TObject::kOverwrite);
  outfile->Close();

  
}


