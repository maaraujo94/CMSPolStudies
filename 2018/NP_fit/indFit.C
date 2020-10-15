// code to do the individual fit
// as this fits each slice on its own, can use the right fit range already

// function to parse a string into components separated by "deli"
vector< string > parseString( string line, string deli) {
  vector< string > out;
  string aux = line;
  size_t pos = 0;
  
  while( (pos = aux.find(deli)) != string::npos)
    {
      out.push_back( aux.substr( 0, pos ) );
      aux.erase( 0, pos + deli.length() );	    
    }
  out.push_back( aux );
  
  return out;
}

// main
void indFit()
{

  TFile *outf = new TFile("files/fit_res_1d.root", "recreate");
  outf->Close();

  int bin = 0;

  TCanvas *c = new TCanvas("","",700,700);
  
  string fileN[3] = {"", "_hpt", "_vhpt"};
  for(int i_pt = 0; i_pt < 3; i_pt++) {
    // read the coarse histos in |costh|
    TFile *infile = new TFile("files/ratioHist.root");
    TH2D *hist = new TH2D();
    infile->GetObject(Form("ratioHist_ab%s", fileN[i_pt].c_str()), hist);
    hist->SetDirectory(0);

    int nBinsX = hist->GetNbinsX(), nBinsY = hist->GetNbinsY();
    const double *yBins = hist->GetYaxis()->GetXbins()->GetArray();

    infile->Close();
  
    // get the 1d plots
    TH1D *pHist[nBinsY];
    for(int i = 1; i <= nBinsY; i++) {
      bin++;
      pHist[i-1] = hist->ProjectionX(Form("fine_bin%d_1d_min", bin), i, i);
      pHist[i-1]->SetTitle(Form("p_{T} bin %d: [%.0f, %.0f] GeV", bin, yBins[i-1], yBins[i]));
    }
  
    // cout some useful info for checking purposes
    cout << "fitting histogram " << hist->GetName() << endl;
    cout << "X axis: " << nBinsX << " bins in [" << hist->GetXaxis()->GetBinLowEdge(1) << "," << hist->GetXaxis()->GetBinUpEdge(nBinsX) << "]" << endl;
    cout << "Y axis: " << nBinsY << " bins in [" << hist->GetYaxis()->GetBinLowEdge(1) << "," << hist->GetYaxis()->GetBinUpEdge(nBinsY) << "]" << endl;

    // the fit function to be used
    TF1 *fit1d = new TF1("fit f 1d", "[0]*(1+[1]*x*x)", 0, 1);
    fit1d->SetParNames("A", "l_th");
    fit1d->SetParameters(1., 0.1);

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
 
    // the cycle to fit each bin and store fit results
    TFile *outfile = new TFile("files/fit_res_1d.root", "update");

    double parA[2][nBinsY], parL[2][nBinsY], chi2[nBinsY], ndf[nBinsY], chiP[nBinsY];
    double bins[nBinsY], ebins[nBinsY], cMaxVal[nBinsY];
 
  for(int i = 0; i < nBinsY; i++) {
    double pMin = hist->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = hist->GetYaxis()->GetBinUpEdge(i+1);
    
    cMaxVal[i] = cosMax->Integral(pMin, pMax)/(pMax-pMin);
    double cR = floor(cMaxVal[i]*10.)/10.;
    if(cMaxVal[i]-cR>0.05) cR += 0.05;
    fit1d->SetRange(0, cR);
    cMaxVal[i] = cR;
    
    pHist[i]->Fit(fit1d, "R");
    parA[0][i] = fit1d->GetParameter(0);
    parA[1][i] = fit1d->GetParError(0);
    parL[0][i] = fit1d->GetParameter(1);
    parL[1][i] = fit1d->GetParError(1);
    chi2[i] = fit1d->GetChisquare();
    ndf[i] = fit1d->GetNDF();
    chiP[i] = TMath::Prob(chi2[i], ndf[i]);
    bins[i] = (pMax+pMin)/2.;
    ebins[i] = (pMax-pMin)/2.;    
    pHist[i]->Write();
  }

  // make and save the TGraph with the fit results and max costh used
  TGraphErrors *graphA = new TGraphErrors(nBinsY, bins, parA[0], ebins, parA[1]);
  TGraphErrors *graphL = new TGraphErrors(nBinsY, bins, parL[0], ebins, parL[1]);
  TGraph *graphC = new TGraph(nBinsY, bins, chi2);
  TGraph *graphN = new TGraph(nBinsY, bins, ndf);
  TGraph *graphP = new TGraph(nBinsY, bins, chiP);
  TGraph *graphCM = new TGraph(nBinsY, bins, cMaxVal);

  graphA->SetName(Form("graph_A%s", fileN[i_pt].c_str()));
  graphL->SetName(Form("graph_lambda%s", fileN[i_pt].c_str()));
  graphC->SetName(Form("graph_chisquare%s", fileN[i_pt].c_str()));
  graphN->SetName(Form("graph_NDF%s", fileN[i_pt].c_str()));
  graphP->SetName(Form("graph_chiP%s", fileN[i_pt].c_str()));
  graphCM->SetName(Form("graph_cosMax%s", fileN[i_pt].c_str()));
 
  graphA->Write();
  graphL->Write();
  graphC->Write();
  graphN->Write();
  graphP->Write();
  graphCM->Write();
  outfile->Close();
  }
  c->Destructor();
}
