// code to do the fit (2d histos for each pT bin)

#import "../ptbins.C"
TH2D **h_NPc = new TH2D*[nPtBins];
TF2 **f_w = new TF2*[nPtBins];

const double gPI = TMath::Pi();

double w_dist(double costh, double ph_deg, double lth, double lph, double ltp)
{
  double phi = ph_deg*gPI/180.;
  ltp = 0;
  
  double th_term = costh*costh;
  double ph_term = sqrt(1-costh*costh)*cos(2.*phi);
  double tp_term = 2*sqrt(1-costh*costh)*costh*costh*cos(phi);

  double w = 3./(4.*gPI*(3.+lth))*(1+lth*th_term+lph*ph_term+ltp*tp_term);
  return w;
}

double myFunction(double lth, double lph, double ltp, double A, double id)
{
  double wval, hval, herr;
  double xmin, xmax, ymin, ymax;
  double chisq = 0, factor;

  int i = (int)id;
  
  int nBinsX = h_NPc[0]->GetNbinsX(), nBinsY = h_NPc[0]->GetNbinsY();
  double minX = h_NPc[0]->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_NPc[0]->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/(nBinsX);
  double minY = h_NPc[0]->GetYaxis()->GetBinLowEdge(1);
  double maxY = h_NPc[0]->GetYaxis()->GetBinUpEdge(nBinsY);
  double dY = (maxY-minY)/(nBinsY);

  f_w[i]->SetParameters(lth, lph, ltp, A);
  
  for(int iX = 0; iX < nBinsX; iX++) {
    xmin = minX+dX*iX;
    xmax = minX+dX*(iX+1);
    for(int iY = 0; iY < nBinsY; iY++) {
      ymin = minY+dY*iY;
      ymax = minY+dY*(iY+1);

      hval = h_NPc[i]->GetBinContent(iX+1, iY+1);
      if(hval==0) continue;

      if(iX+2 < nBinsX)
	if(h_NPc[i]->GetBinContent(iX+2, iY+1) == 0)
	  continue;
      if(iY+2 < nBinsX)
	if(h_NPc[i]->GetBinContent(iX+1, iY+2) == 0)
	  continue;
      if(iX > 0)
	if(h_NPc[i]->GetBinContent(iX, iY+1) == 0)
	  continue;
      if(iY > 0)
	if(h_NPc[i]->GetBinContent(iX+1, iY) == 0)
	  continue;
      
      wval = f_w[i]->Integral(xmin, xmax, ymin, ymax); 
      herr = h_NPc[i]->GetBinError(iX+1, iY+1);

      factor = pow((wval-hval)/herr, 2);
      chisq += factor;
    }
  }

  return chisq;
}

//Minuit-specific "wrapping"
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg)
{
  result=myFunction(par[0], par[1], par[2], par[3], par[4]);
}


// main
void doFit_ltpz()
{
  TH2D **h_fit = new TH2D*[nPtBins];
  TH2D **h_pull= new TH2D*[nPtBins];
  
  // read the histos from subtraction
  TFile *infile = new TFile("files/bkgSubRes.root");
  for(int i = 0; i < nPtBins; i++) {
    infile->GetObject(Form("h_NPcs_%d", i), h_NPc[i]);
    h_NPc[i]->SetDirectory(0);
  }
  infile->Close();

  int nBinsX = h_NPc[0]->GetNbinsX(), nBinsY = h_NPc[0]->GetNbinsY();
  double minX = h_NPc[0]->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_NPc[0]->GetXaxis()->GetBinUpEdge(nBinsX);
  double minY = h_NPc[0]->GetYaxis()->GetBinLowEdge(1);
  double maxY = h_NPc[0]->GetYaxis()->GetBinUpEdge(nBinsY);
  double dX = (maxX-minX)/(nBinsX);
  double dY = (maxY-minY)/(nBinsY);

  // the fit function to be used
  for(int i = 0; i < nPtBins; i++) {
    f_w[i] = new TF2(Form("f_w_%d", i), "[3]*w_dist(x,y,[0],[1],[2])", 0, 1, 0, 90);
    f_w[i]->SetParNames("l_th", "l_ph", "l_tp", "A");

    //also define the fit ratio histo
    h_fit[i] = new TH2D(Form("fit_ratio_%d", i), "Data/Fit Ratio", nBinsX, minX, maxX, nBinsY, minY, maxY);
    h_pull[i] = new TH2D(Form("fit_pull_%d", i), "Fit Pulls", nBinsX, minX, maxX, nBinsY, minY, maxY);
  }

  // the cycle to fit each bin and store fit results
  TCanvas *c = new TCanvas("", "", 700, 700);    
  TFile *outfile = new TFile("files/finalFitRes_ltpz.root", "recreate");

  double parA[nPtBins], eparA[nPtBins];
  double parL[3][nPtBins], eparL[3][nPtBins];
  double chi2[nPtBins], ndf[nPtBins], chiP[nPtBins];
  double pt[nPtBins], ept[nPtBins];

  TFitter* fit = new TFitter(5);
  fit->SetFCN(minuitFunction);
  
  for(int i = 0; i < nPtBins; i++) {
    //  for(int i = 0; i < 1; i++) {
    double in_a = h_NPc[i]->GetBinContent(1, 1);
    
    fit->SetParameter(0, "l_th", 0.1, 0.01, -1, 1);
    fit->SetParameter(1, "l_ph", -0.1, 0.01, -1, 1);
    fit->SetParameter(2, "l_tp", 0.0, 0.01, -1, 1);
    fit->FixParameter(2);
    fit->SetParameter(3, "A", in_a, in_a/10., 0, 0);
    fit->SetParameter(4, "i_pt", i, i, 0, 0);
    fit->FixParameter(4);
    
    // get pt vars
    double pMin = ptBins[i];
    double pMax = ptBins[i+1];
    pt[i] = (pMax+pMin)/2.;
    ept[i] = (pMax-pMin)/2.;

    fit->ExecuteCommand("MIGRAD",0,0);
    fit->ExecuteCommand("MIGRAD",0,0);

    for(int j = 0; j< 3; j++) {
      parL[j][i] = fit->GetParameter(j);
      eparL[j][i] = fit->GetParError(j);
    }
    parA[i] = fit->GetParameter(3);
    eparA[i] = fit->GetParError(3);
    chi2[i] = myFunction(parL[0][i], parL[1][i], parL[2][i], parA[i], i);
    f_w[i]->SetParameters(parL[0][i], parL[1][i], parL[2][i], parA[i]);

    // ndf will require the bin removal first
    ndf[i] = 0;
    double xmin, xmax, ymin, ymax;
    for(int iX = 0; iX < nBinsX; iX++) {
      xmin = minX+dX*iX;
      xmax = minX+dX*(iX+1);
      for(int iY = 0; iY < nBinsY; iY++) {
	ymin = minY+dY*iY;
	ymax = minY+dY*(iY+1);
	if(h_NPc[i]->GetBinContent(iX+1, iY+1)==0) continue;
	if(iX+2 < nBinsX)
	  if(h_NPc[i]->GetBinContent(iX+2, iY+1) == 0)
	    continue;
	if(iY+2 < nBinsX)
	  if(h_NPc[i]->GetBinContent(iX+1, iY+2) == 0)
	    continue;
	if(iX > 0)
	  if(h_NPc[i]->GetBinContent(iX, iY+1) == 0)
	    continue;
	if(iY > 0)
	  if(h_NPc[i]->GetBinContent(iX+1, iY) == 0)
	    continue;
	ndf[i]+=1;

	//get histo of data/fit ratio
	double hval = h_NPc[i]->GetBinContent(iX+1, iY+1);
	double  herr = h_NPc[i]->GetBinError(iX+1, iY+1);
	double wval = f_w[i]->Integral(xmin, xmax, ymin, ymax);      

	h_fit[i]->SetBinContent(iX+1, iY+1, hval/wval);
	h_fit[i]->SetBinError(iX+1, iY+1, herr/wval);
	h_pull[i]->SetBinContent(iX+1, iY+1, (hval-wval)/herr);
	h_pull[i]->SetBinError(iX+1, iY+1, 0);
      }
    }
    ndf[i] -= fit->GetNumberFreeParameters();
    chiP[i] = TMath::Prob(chi2[i], ndf[i]);

    // plotting everything - will need to adapt fit reuslts since its 2d
    f_w[i]->Write();
    h_NPc[i]->Write();
    h_fit[i]->Write();
    h_pull[i]->Write();

    cout << endl << endl;
  }
  
  // make and save the TGraph with the fit results
     TGraphErrors *graphA = new TGraphErrors(nPtBins, pt, parA, ept, eparA);
  TGraphErrors **graphL = new TGraphErrors*[3];
  for(int j = 0; j < 3; j++) {
    graphL[j] = new TGraphErrors(nPtBins, pt, parL[j], ept, eparL[j]);
  }
  TGraph *graphC = new TGraph(nPtBins, pt, chi2);
  TGraph *graphN = new TGraph(nPtBins, pt, ndf);
  TGraph *graphP = new TGraph(nPtBins, pt, chiP);

  string lbl[] = {"theta", "phi", "thph"};
  
  graphA->SetName(Form("graph_A"));
  for(int j = 0; j < 3; j++) 
    graphL[j]->SetName(Form("graph_lambda_%s", lbl[j].c_str()));
  graphC->SetName(Form("graph_chisquare"));
  graphN->SetName(Form("graph_NDF"));
  graphP->SetName(Form("graph_chiP"));
 
  graphA->Write();
  for(int j = 0; j < 3; j++) 
    graphL[j]->Write();
  graphC->Write();
  graphN->Write();
  graphP->Write();
    
  outfile->Close();

  c->Destructor();
}
