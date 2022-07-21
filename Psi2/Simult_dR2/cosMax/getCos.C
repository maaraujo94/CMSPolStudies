// macro to get the maximum |costh| in each pT bin that we can fit

TF1 *fitf;

double jumpF(double pt)
{
  double b_val = fitf->Eval(pt);
  
  double b_round = floor(b_val*10.)/10.;
  if(b_val-b_round > 0.025) {
    b_round += 0.05;
    if(b_val - b_round > 0.025)
      b_round += 0.05;
  }

  return b_round;
}

void getCos()
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
  for(int i = 0; i < nBinsY; i++) {
    int nBinsX = pHist[i]->GetNbinsX();
    for(int j = 1; j < nBinsX; j++) {
      if(pHist[i]->GetBinError(j+2) > 2.*pHist[i]->GetBinError(j+1) && i == 0) {
	cosMax[i] = pHist[i]->GetBinLowEdge(j+2);
	cout << i << " " << yBins[i] << " " << yBins[i+1] << " " << cosMax[i] << endl;
	break;
      }
      else if(pHist[i]->GetBinError(j+2) > 2.*pHist[i]->GetBinError(j+1) && i > 0 && pHist[i]->GetBinLowEdge(j+2) >= cosMax[i-1]) {
	cosMax[i] = pHist[i]->GetBinLowEdge(j+2);
	cout << i << " " << yBins[i] << " " << yBins[i+1] << " " << cosMax[i] << endl;
	break;
      }
      else if(i>0) {
	cosMax[i] = cosMax[i-1];
      }
    }
  }
  
  fout.open("cos_max.txt", std::ofstream::app);
  for(int i = 0; i < nBinsY; i++) {
    if(i == 0 && cosMax[i] >= aux)
      fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    else if(i == 0 && cosMax[i] < aux) {
      cosMax[i] = aux;
      fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    }
    else if(cosMax[i] >= cosMax[i-1])
      fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    else {
      cosMax[i] = cosMax[i-1];
      fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    }
    aux = cosMax[i];
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
  
  // the fit function: a logarithm
  fitf = new TF1("fitf", "[0]*log([1]+[2]*x)", yBins[0], yBins[nBinsY]);
  fitf->SetParameters(0.1, -10, 8);
  fitf->SetLineColor(kBlue);
  costh->Fit("fitf");
  costh->Draw("");
  cout << "chi^2/ndf = " << fitf->GetChisquare() << "/" << fitf->GetNDF() << endl;

  TF1 *fit_j = new TF1("fit_j", "jumpF(x)", yBins[0], yBins[nBinsY]);
  fit_j->SetLineColor(kGreen);
  fit_j->Draw("same");

  
  can->SaveAs("costh_lim.pdf");
  can->Clear();

  // save the fit results to a txt file
  ofstream outfile;
  outfile.open("cosMaxFitRes.txt");
  outfile << "[a]*log([b]+[c]*pT)" << endl;
  outfile << "a\t e_a\t b\t e_b\t c\t e_c\t chi2\t ndf" << endl;
  for(int i = 0; i < 3; i++)
    outfile << fitf->GetParameter(i) << "\t" << fitf->GetParError(i) << "\t";
  outfile << fitf->GetChisquare() << "\t" << fitf->GetNDF() << endl;
  outfile.close();

  outfile.open("cosMaxFitRes.tex");
  outfile << "\\begin{tabular}{c|c|c|c}\n";
  outfile << "$a$ & $b$ & $c$ & $\\chi^2$/ndf \\\\\n";
  outfile << "\\hline\n";
  for(int i = 0; i < 3; i++) {
    int prec = ceil(-log10(fitf->GetParError(i)))+1;
    outfile << "$" << setprecision(prec) << fixed << fitf->GetParameter(i) << "\\pm" << fitf->GetParError(i) << "$ & ";
  }
  outfile << setprecision(0) << fixed <<  fitf->GetChisquare() << "/" << fitf->GetNDF() << endl;
  outfile << "\\end{tabular}\n";
  outfile.close();
}
