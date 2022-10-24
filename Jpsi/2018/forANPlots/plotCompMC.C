#import "../cosMax/imp_jumpF.C"

void plotCompMC()
{
  // get the ratio and base MCs
  TFile *fin = new TFile("files/MCcompStore.root");
  int n_c = 3;
  int n_b = 6;
  TH2D **hr_MC = new TH2D*[n_c];
  TH2D **h_MC = new TH2D*[n_b];
  for(int i = 0; i < n_c; i++) {
    hr_MC[i] = (TH2D*)fin->Get(Form("ratio_%d", i));
    hr_MC[i]->SetDirectory(0);
  }
  for(int i = 0; i < n_b; i++) {
    h_MC[i] = (TH2D*)fin->Get(Form("mcH%d", i));
    h_MC[i]->SetDirectory(0);
  }
  fin->Close();

  //set the colors for the histos
  int colpt[] = {kRed+1, kViolet+2, kViolet+2, kGreen+3, kGreen+3, kBlue};
  
  // get the fit range from our cosmax(pT)
  ifstream in;
  string dataS;
  in.open("../cosMax/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", 20, 130);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);

  // constant function for check
  TF1 *fc = new TF1("fc", "[0]", 0, 1);

  // plotting one by one for each interval
  TCanvas *c = new TCanvas("", "", 700, 700);    

  for(int i = 0; i < n_c; i++) {
    // get the binning
    int nBinsX = hr_MC[i]->GetNbinsX(), nBinsY = hr_MC[i]->GetNbinsY();
    const double *yBins = hr_MC[i]->GetYaxis()->GetXbins()->GetArray();
    
    for(int i_y = 0; i_y < nBinsY; i_y++)
      {
	// first the ratios
	// get the 1d projection
	TH1D *pHist = hr_MC[i]->ProjectionX(Form("bin%d_%d", i, i_y), i_y+1, i_y+1);
	pHist->SetTitle(Form("MC ratio %.0f < p_{T} < %.0f GeV", yBins[i_y], yBins[i_y+1]));
	pHist->SetStats(0);
	pHist->SetLineColor(kBlack);
	pHist->SetMarkerColor(kBlack);
	pHist->SetMinimum(0);
	pHist->SetMaximum(pHist->GetBinContent(1)*1.5);
	pHist->GetXaxis()->SetTitle("|cos#theta_{HX}|");
	pHist->Draw("error");

	// get max costheta
	double cMaxVal = jumpF(cosMax->Integral(yBins[i_y], yBins[i_y+1])/(yBins[i_y+1]-yBins[i_y]));
	
	// fit to constant
	fc->SetRange(0, cMaxVal);
	fc->SetParameter(0, pHist->GetBinContent(1));
	pHist->Fit("fc", "R");

	TLine *c_lim = new TLine(cMaxVal, 0, cMaxVal, pHist->GetMaximum());
	c_lim->SetLineStyle(kDashed);
	c_lim->SetLineColor(kBlack);
	c_lim->Draw();

	c->SaveAs(Form("plots/compMC/bin%d_%d.pdf", i, i_y));
	c->Clear();

	// then the base dists
	// get the 1d projections
	TH1D *pHistL = h_MC[i*2]->ProjectionX(Form("base1_%d_%d", i, i_y), i_y+1, i_y+1);
	pHistL->SetTitle(Form("MC distributions %.0f < p_{T} < %.0f GeV", yBins[i_y], yBins[i_y+1]));
	pHistL->SetStats(0);
	pHistL->SetLineColor(colpt[i*2]);
	pHistL->SetMarkerColor(colpt[i*2]);
	pHistL->Scale(1./pHistL->Integral());
	pHistL->SetMinimum(0);
	pHistL->SetMaximum(pHistL->GetBinContent(1)*1.5);
	pHistL->GetXaxis()->SetTitle("|cos#theta_{HX}|");
	pHistL->Draw("error");

	TH1D *pHistH = h_MC[i*2+1]->ProjectionX(Form("base2_%d_%d", i, i_y), i_y+1, i_y+1);
	pHistH->SetStats(0);
	pHistH->SetLineColor(colpt[i*2+1]);
	pHistH->SetMarkerColor(colpt[i*2+1]);
	pHistH->Scale(1./pHistH->Integral());
	pHistH->Draw("error same");

	TLine *c_limB = new TLine(cMaxVal, 0, cMaxVal, pHistL->GetMaximum());
	c_limB->SetLineStyle(kDashed);
	c_limB->SetLineColor(kBlack);
	c_limB->Draw();

	c->SaveAs(Form("plots/compMC/binB%d_%d.pdf", i, i_y));
	c->Clear();

      }
  }

  c->Destructor();
  
}
