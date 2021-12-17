// macro to get the maximum |costh| in each pT bin that we can fit

TF1 *fitf;

double jumpF(double pt, double pt_z)
{
  if(pt<pt_z) return 0;
  
  double b_val = fitf->Eval(pt);
  
  double b_round = floor(b_val*10.)/10.;
  if(b_val-b_round > 0.025) {
    b_round += 0.05;
    if(b_val - b_round > 0.025)
      b_round += 0.05;
  }

  return b_round;
}

void getCosMin()
{
  double aux = 0;

  ofstream fout;
  fout.open("cos_max.txt");
  fout.close();  
  
  // read the coarse histos in |costh|
  TFile *infile = new TFile("histoStore.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("ratioH_ab"), hist);
  hist->SetDirectory(0);
  infile->Close();

  int nBinsY = hist->GetNbinsY();
  const double* yBins = hist->GetYaxis()->GetXbins()->GetArray();  
  
  TH1D *pHist[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    pHist[i] = hist->ProjectionX(Form("fine_bin%d_1d_ab", i+1), i+1, i+1);
  }
    
  double cosMax[nBinsY];
  for(int i = nBinsY-1; i >= 0; i--) {
    int nBinsX = pHist[i]->GetNbinsX();
    for(int j = 1; j < nBinsX; j++) {
      if(pHist[i]->GetBinError(j) > 2.*pHist[i]->GetBinError(j+1) && i == nBinsY-1) {
	cosMax[i] = pHist[i]->GetBinLowEdge(j+1);
	cout << i+1 << " " << yBins[i] << " " << yBins[i+1] << " " << cosMax[i] << endl;
	break;
      }
      else if(pHist[i]->GetBinError(j) > 2.*pHist[i]->GetBinError(j+1) && i < nBinsY-1 && pHist[i]->GetBinLowEdge(j+1) <= cosMax[i+1]) {
	cosMax[i] = pHist[i]->GetBinLowEdge(j+1);
	cout << i+1 << " " << yBins[i] << " " << yBins[i+1] << " " << cosMax[i] << endl;
	break;
      }
      else if(i<nBinsY-1 && cosMax[i+1] == 0.05) {
	cosMax[i] = 0;
      }
      else if(i<nBinsY-1) {
	cosMax[i] = cosMax[i+1];
      }
    }
  }

  fout.open("cos_min.txt", std::ofstream::app);
  for(int i = 0; i < nBinsY; i++) {
    fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    if(cosMax[i] == 0) aux = yBins[i+1];
  }
  fout.close();
  
  // plot the results in a histogram and fit it
  TCanvas *can = new TCanvas("", "", 700, 700);
  
  // the histogram
  TH1F* costh = new TH1F("name", "cos#theta limits", nBinsY, yBins);
  costh->GetXaxis()->SetTitle("p_{T} (GeV)");
  costh->GetYaxis()->SetTitle("|cos#theta_{HX}|_{max}");
  for(int i = 0; i < nBinsY; i++) {
    costh->SetBinContent(i+1, cosMax[i]);
    costh->SetBinError(i+1, 0.025);
  }
  costh->SetLineColor(kBlack);
  costh->SetStats(0);
  costh->Draw("");
  
  // the fit function: a logarithm
  // fitf = new TF1("fitf", "[0]*log([1]+[2]*x)", aux, yBins[nBinsY]);
  // fitf = new TF1("fitf", "[0]*sqrt(abs([1]+[2]*x))", aux, yBins[nBinsY]);
  fitf = new TF1("fitf", "[0]*(1-exp([1]+[2]*x))", aux, yBins[nBinsY]);
  fitf->SetParameters(10, -1, -6e-2);
  fitf->SetLineColor(kBlue);
  costh->Fit("fitf", "R");
  costh->Draw();
  cout << "chi^2/ndf = " << fitf->GetChisquare() << "/" << fitf->GetNDF() << endl;

  TF1 *fit_j = new TF1("fit_j", "jumpF(x,[0])", yBins[0], yBins[nBinsY]);
  fit_j->SetParameter(0, aux);
  fit_j->SetLineColor(kGreen);
  fit_j->Draw("same");

  
  can->SaveAs("costh_min.pdf");
  can->Clear();
  
  // save the fit results to a txt file
  ofstream outfile;
  outfile.open("cosMinFitRes.txt");
  outfile << "[a]*(1-exp([b]+[c]*pT), pT > [d]" << endl;
  outfile << "a\t e_a\t b\t e_b\t c\t e_c\t d\t chi2\t ndf" << endl;
  for(int i = 0; i < 3; i++)
    outfile << fitf->GetParameter(i) << "\t" << fitf->GetParError(i) << "\t";
  outfile << aux << "\t";
  outfile << fitf->GetChisquare() << "\t" << fitf->GetNDF() << endl;
  outfile.close();

  outfile.open("cosMinFitRes.tex");
  outfile << "\\begin{tabular}{c|c|c|c|c}\n";
  outfile << "$a$ & $b$ & $c$ & $d$ & $\\chi^2$/ndf \\\\\n";
  outfile << "\\hline\n";
  for(int i = 0; i < 3; i++) {
    int prec = ceil(-log10(fitf->GetParError(i)))+1;
    outfile << "$" << setprecision(prec) << fixed << fitf->GetParameter(i) << "\\pm" << fitf->GetParError(i) << "$ & ";
  }
  outfile << setprecision(0) << fixed << aux << " & ";
  outfile << setprecision(0) << fixed <<  fitf->GetChisquare() << "/" << fitf->GetNDF() << endl;
  outfile << "\\end{tabular}\n";
  outfile.close();
}
