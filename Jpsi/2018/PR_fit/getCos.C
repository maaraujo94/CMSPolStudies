// macro to get the maximum |costh| in each pT bin that we can fit

void getCos()
{
  double aux = 0;

  ofstream fout;
  fout.open("text_output/cos_max.txt");
  fout.close();  
  
  // read the coarse histos in |costh|
  TFile *infile = new TFile("files/ratioHist.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("ratioHist_ab"), hist);
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
      if(pHist[i]->GetBinError(j+2) > 2.*pHist[i]->GetBinError(j+1)) {
	cosMax[i] = pHist[i]->GetBinLowEdge(j+2);
	cout << i << " " << cosMax[i] << endl;
	break;
      }
    }
  }
  
  fout.open("text_output/cos_max.txt", std::ofstream::app);
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
  TF1 *fitf = new TF1("fitf", "[0]*log([1]+[2]*x)", yBins[0], yBins[nBinsY]);
  fitf->SetParameters(0.1, -10, 8);
  fitf->SetLineColor(kBlue);
  costh->Fit("fitf");
  costh->Draw("");
  cout << "chi^2/ndf = " << fitf->GetChisquare() << "/" << fitf->GetNDF() << endl;
  
  can->SaveAs("plots/costh_lim.pdf");
  can->Clear();

  // save the fit results to a txt file
  ofstream outfile;
  outfile.open("text_output/cosMaxFitRes.txt");
  outfile << "[a]*log([b]+[c]*pT/M)" << endl;
  outfile << "a\t e_a\t b\t e_b\t c\t e_c\t chi2\t ndf" << endl;
  for(int i = 0; i < 3; i++)
    outfile << fitf->GetParameter(i) << "\t" << fitf->GetParError(i) << "\t";
  outfile << fitf->GetChisquare() << "\t" << fitf->GetNDF() << endl;
  outfile.close();

  outfile.open("text_output/cosMaxFitRes.tex");
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
