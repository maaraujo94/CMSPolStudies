// fitting the function that will then give the costh_max
// currently commented out: plotting the obtained limits in the 1D histos

void fitCosMax()
{
  double M_q = 3.097;
  
  // get .txt where I stored each costh_max
  ifstream cosmax;
  string data;
  cosmax.open("text_output/cos_max.txt");

  const int nbins = 29;
  double costhLim[nbins], ptLims[nbins+1], aux;
  for(int i = 0; i < nbins; i++) {
    cosmax >> ptLims[i] >> aux >> costhLim[i];
    ptLims[i] /= M_q;
  }
  cosmax.close();
  ptLims[nbins] = aux/M_q;

  // plot the results in a histogram and fit it
  TCanvas *can = new TCanvas("", "", 700, 700);
  
  // the histogram
  TH1F* costh = new TH1F("name", "cos#theta limits", nbins, ptLims);
  costh->GetXaxis()->SetTitle("p_{T}/M");
  costh->GetYaxis()->SetTitle("|cos#theta_{HX}|_{max}");
  for(int i = 0;i < nbins; i++) {
    costh->SetBinContent(i+1, costhLim[i]);
    costh->SetBinError(i+1, 0.025);
  }
  costh->GetXaxis()->SetRangeUser(0, 25);
  costh->SetLineColor(kBlack);
  costh->SetStats(0);
  
  // the fit function: a logarithm
  TF1 *fitf = new TF1("fitf", "[0]*log([1]+[2]*x)", 12/M_q, 70/M_q);
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
