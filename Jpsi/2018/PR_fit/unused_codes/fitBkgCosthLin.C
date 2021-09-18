#import "../cosMax/imp_jumpF.C"

//pt bins defined globally for access from functions
const int nPtBins = 7;
double ptBins[nPtBins+1];

TF1 *cosMax;

// tf2 function parser
// parameters: N (per bin), lambda2 (linear), lambda4 (linear)
double sum_func(double *xx, double *pp)
{
  // get costh, pt and corresp pt bin
  double cos = xx[0], pt = xx[1];
  int pt_bin;
  for(int i = 0; i < nPtBins; i++)
    if(ptBins[i] <= pt && ptBins[i+1] >= pt)
      pt_bin = i;

  
  double pMin = ptBins[pt_bin], pMax = ptBins[pt_bin+1];
  double pAv = 0.5 * (pMax + pMin);
  
  // get costh region where fit is done
  double cMaxVal = jumpF(cosMax->Integral(pMin, pMax)/(pMax-pMin));
  if(cos > cMaxVal)
    {
      TF2::RejectPoint();
      return 0;
    }

  // get fit parameters from input
  double N = pp[pt_bin];
  double ld2 = pp[nPtBins] * pAv + pp[nPtBins+1];
  double ld4 = pp[nPtBins+2] * pAv + pp[nPtBins+3];

  double func = N * (1 + ld2 * pow(cos,2) + ld4 * pow(cos,4));
  return func;
}

// macro to fit the sideband/MC distributions
void fitBkgCosthLin()
{
  string lbl[] = {"LSB", "RSB"};
  double dev_min = 30, pull_min = 3;

  // get the SB/MC 2d maps from the repository
  TH2D **h_cth = new TH2D*[2];
  TFile *fIn = new TFile("files/bkgHist.root");
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
  for(int i = 0; i <= nPtBins; i++) ptBins[i] = yBins[i];

  // get the fit range from our cosmax(pT)
  ifstream in;
  string dataS;
  in.open("../cosMax/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", yBins[0]-10, yBins[nBinsY]+10);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);
  
  // cycle over 2 background types
  TCanvas *c = new TCanvas("", "", 700, 700);
  for(int i_inp = 0; i_inp < 2; i_inp++) {
    
    // get the 1D bins
    TH1D *pHist[nBinsY];
    // getting the pT bin projections of Data/MC
    for(int i = 0; i < nPtBins; i++) {
      pHist[i] = h_cth[i_inp]->ProjectionX(Form("%s_bin%d_1d", lbl[i_inp].c_str(), i+1), i+1, i+1);
      pHist[i]->SetTitle(Form("2018 %s/MC bin %d: [%.0f, %.0f] GeV", lbl[i_inp].c_str(), i+1, yBins[i], yBins[i+1]));
    }

    // define the ratio fit function
    TF2 *f_fit = new TF2("f_fit", sum_func, minX, maxX, ptBins[0], ptBins[nPtBins], nPtBins+4, 2);
    // define free parameters - N
    for(int i = 0; i < nPtBins; i++) {
      f_fit->SetParName(i, Form("N_%d", i));
      f_fit->SetParameter(i, pHist[i]->GetBinContent(1)*1.1);
    }
    // define linear parameters - lambda_2, lambda_4
    f_fit->SetParName(nPtBins, "m_ld_2");
    f_fit->SetParameter(nPtBins, 0.01);
    f_fit->SetParName(nPtBins+1, "b_ld_2");
    f_fit->SetParameter(nPtBins+1, 0.1);
    f_fit->SetParName(nPtBins+2, "m_ld_4");
    f_fit->SetParameter(nPtBins+2, 0.01);
    f_fit->SetParName(nPtBins+3, "b_ld_4");
    f_fit->SetParameter(nPtBins+3, 0.1);
    // fit the 2d function to the costh:pT map
    h_cth[i_inp]->Fit("f_fit");

    // tf1 for plotting in the 1D bins
    TF1 *f_1d = new TF1("f_1d", "[0]*(1+[1]*x*x+[2]*pow(x,4))", minX, maxX);
    f_1d->SetParNames("N", "l_2", "l_4");
    
    double par[3][nBinsY], epar[3][nBinsY];
    double cMaxVal[nBinsY], pt[nBinsY], ept[nBinsY];

    // cycle over all pT bins
    for(int i = 0; i < nBinsY; i++) {

      pt[i] = 0.5*(yBins[i+1]+yBins[i]);
      ept[i] = 0.5*(yBins[i+1]-yBins[i]);

      // storing parameters
      par[0][i] = f_fit->GetParameter(i);
      epar[0][i] = f_fit->GetParError(i);
      par[1][i] = f_fit->GetParameter(nPtBins) * pt[i] + f_fit->GetParameter(nPtBins+1);
      epar[1][i] = sqrt(pow(f_fit->GetParError(nPtBins) * pt[i], 2) + pow(f_fit->GetParError(nPtBins+1), 2));
      par[2][i] = f_fit->GetParameter(nPtBins+2) * pt[i] + f_fit->GetParameter(nPtBins+3);
      epar[2][i] = sqrt(pow(f_fit->GetParError(nPtBins+2) * pt[i], 2) + pow(f_fit->GetParError(nPtBins+3), 2));
      
      // getting the max costh value for the fit
      double pMin = h_cth[i_inp]->GetYaxis()->GetBinLowEdge(i+1);
      double pMax = h_cth[i_inp]->GetYaxis()->GetBinUpEdge(i+1);
      
      cMaxVal[i] = jumpF(cosMax->Integral(pMin, pMax)/(pMax-pMin));

      // initializing f_1d and plotting
      f_1d->SetRange(0, cMaxVal[i]);
      f_1d->SetParameters(par[0][i], par[1][i], par[2][i]);
      
      // plot the fit
      pHist[i]->SetStats(0);
      pHist[i]->SetMinimum(0);
      pHist[i]->SetMaximum(par[0][i]*2.5);
      pHist[i]->SetMarkerColor(kBlack);
      pHist[i]->SetLineColor(kBlack);
      pHist[i]->GetXaxis()->SetTitle("|cos#theta|");
      pHist[i]->Draw();
      f_1d->SetLineColor(kBlue);
      f_1d->Draw("same");
    
      c->SaveAs(Form("plots/%slin/fit_pt%d.pdf", lbl[i_inp].c_str(), i));
      c->Clear();

      // calculating pulls
      double cv[nBinsX], pv[nBinsX], dv[nBinsX];
      for(int i_cos = 0 ; i_cos < nBinsX; i_cos++) {
	double cMin = h_cth[i_inp]->GetXaxis()->GetBinLowEdge(i_cos+1);
	double cMax = h_cth[i_inp]->GetXaxis()->GetBinUpEdge(i_cos+1);
	cv[i_cos] = 0.5 * (cMax+cMin);
	
	double fitv = f_1d->Eval(cv[i_cos]);
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
      fp->SetTitle(Form("2018 %s pulls (%.0f < p_{T} < %.0f GeV)", lbl[i_inp].c_str(), yBins[i], yBins[i+1]));

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
	
      c->SaveAs(Form("plots/%slin/pulls_pt%d.pdf", lbl[i_inp].c_str(), i));
      c->Clear();

      // plot the deviations
      TH1F *fd = c->DrawFrame(minX, -dev_min, maxX, dev_min);
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
		
      TLine *climd = new TLine(cMaxVal[i], -dev_min, cMaxVal[i], dev_min);
      climd->SetLineColor(kRed);
      climd->SetLineStyle(kDashed);
      climd->Draw();

      zero->Draw();
	
      c->SaveAs(Form("plots/%slin/devs_pt%d.pdf", lbl[i_inp].c_str(), i));
      c->Clear();
    }

    // storing the free parameters
    TFile *fOut = new TFile(Form("files/%slin_fitres.root", lbl[i_inp].c_str()), "recreate");
    string parlab[] = {"N", "l2", "l4"};
    
    for(int i_p = 0; i_p < 3; i_p++) {
      TGraphErrors *g_par = new TGraphErrors(nPtBins, pt, par[i_p], ept, epar[i_p]);
      g_par->Write(Form("fit_%s", parlab[i_p].c_str()));
    }
    
    TLine *l_chi = new TLine(yBins[0], TMath::Prob(f_fit->GetChisquare(), f_fit->GetNDF()), yBins[nBinsY], TMath::Prob(f_fit->GetChisquare(), f_fit->GetNDF()));
    l_chi->Write(Form("fit_chiP"));

    // parameters of linear function for lambdas
    double lv2[2] = {f_fit->GetParameter(nPtBins), f_fit->GetParameter(nPtBins+1)};
    double lv4[2] = {f_fit->GetParameter(nPtBins+2), f_fit->GetParameter(nPtBins+3)};
    TGraphErrors *l_l2 = new TGraphErrors(2, lv2, lv2, lv2, lv2);
    TGraphErrors *l_l4 = new TGraphErrors(2, lv4, lv4, lv4, lv4);
    l_l2->Write(Form("ld2_lin"));
    l_l4->Write(Form("ld4_lin"));
    
    fOut->Close();
 
    // storing the fit results in tex format
    /* ofstream fout;
    fout.open(Form("text_output/cos_%s2d.tex", lbl[i_inp].c_str()));
    fout << "\\begin{tabular}{c||c}\n";
    fout << "$\\pt$ (GeV) & $N$  \\\\\n";
    fout << "\\hline\n";
    for(int i = 0; i < nBinsY; i++) {
      double pMin = h_cth[i_inp]->GetYaxis()->GetBinLowEdge(i+1);
      double pMax = h_cth[i_inp]->GetYaxis()->GetBinUpEdge(i+1);
      fout << Form("$[%.0f, %.0f]$", pMin, pMax);
      if(epar[0][i] > 0) {
	int p_norm = 1.;
	if(epar[0][i] < 1 ) 
	  p_norm = ceil(-log10(epar[0][i]))+1;
	fout << " & " <<  setprecision(p_norm) << fixed << par[0][i] << " $\\pm$ " << epar[0][i];
	  }
	else
	  fout << " & " <<  setprecision(2) << fixed << par[0][i];
      fout  << "\\\\\n";
    }
    fout << "\\end{tabular}\n";
    fout.close();
 
    ofstream fout2;
    fout2.open(Form("text_output/cosA_%s2d.tex", lbl[i_inp].c_str()));
    fout2 << "\\begin{tabular}{c|c||c}\n";
    fout2 << "$\\lambda_{2}$ & $\\lambda_4$ & $\\chi^{2}$/ndf \\\\\n";
    fout2 << "\\hline\n";
    // lambda_2
    double val = f_fit->GetParameter(nPtBins);
    double unc = f_fit->GetParError(nPtBins);
    int p_norm = 1.;
    if(unc < 1) p_norm = ceil(-log10(unc))+1;	
    fout2 <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
    // lambda_4
    val = f_fit->GetParameter(nPtBins+1);
    unc = f_fit->GetParError(nPtBins+1);
    p_norm = 1.;
    if(unc < 1) p_norm = ceil(-log10(unc))+1;	
    fout2 <<  setprecision(p_norm) << fixed << val << " $\\pm$ " << unc << " & ";
    // chi^2
    fout2 << setprecision(0) << f_fit->GetChisquare() << "/" << f_fit->GetNDF() << "\\\\\n";
    fout2 << "\\end{tabular}\n";
    fout2.close(); */

    cout<< endl << "finished fitting " << lbl[i_inp] << endl << endl;
  }
  c->Destructor();

}
