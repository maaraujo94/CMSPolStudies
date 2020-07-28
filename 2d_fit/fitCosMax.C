// fitting the function that will then give the costh_max
// currently commented out: plotting the obtained limits in the 1D histos

void fitCosMax()
{
  // get .txt where I stored each costh_max
  ifstream cosmax;
  string data;
  cosmax.open("cos_max.txt");

  const int nbins = 29;
  double costhLim[nbins], ptLims[nbins+1], aux;
  getline(cosmax, data);
  getline(cosmax, data);
  for(int i = 0; i < nbins; i++) {
    cosmax >> ptLims[i] >> aux >> costhLim[i];
  }
  cosmax.close();
  ptLims[nbins] = aux;

  // plot the results in a histogram and fit it
  TCanvas *can = new TCanvas("", "", 700, 700);
  
  // the histogram
  TH1F* costh = new TH1F("name", "cos#theta limits", nbins, ptLims);
  costh->GetXaxis()->SetTitle("p_{T} (GeV)");
  costh->GetYaxis()->SetTitle("|cos#theta_{HX}|_{max}");
  for(int i = 0;i < nbins; i++) {
    costh->SetBinContent(i+1, costhLim[i]);
    costh->SetBinError(i+1, 0.025);
  }
  costh->GetXaxis()->SetRangeUser(10, 70);
  costh->SetLineColor(kBlack);
  costh->SetStats(0);
  
  // the fit function: a logarithm
  TF1 *fitf = new TF1("fitf", "[0]*log([1]+[2]*x)", 12, 70);
  fitf->SetParameters(0.1, -10, 1);
  fitf->SetLineColor(kBlue);
  costh->Fit("fitf");
  costh->Draw("");
  cout << "chi^2/ndf = " << fitf->GetChisquare() << "/" << fitf->GetNDF() << endl;

  can->SaveAs("plots/costh_lim.pdf");
  can->Clear();

  // read fine binned 2D histo and get the 1D histos w the costh cut
  /*  TFile *infile = new TFile("ratioHist.root");
  TH2D *hist = new TH2D();
  infile->GetObject("ratioHist", hist);
  hist->SetDirectory(0);
  infile->Close();

  // Get the 1D projection in the pT bins (with costh lims marked)
  // coarse binned
  int ptBins[9]={0, 4, 8, 10, 13, 18, 23, 27, 29};
  int binEdges[9] = {12, 16, 20, 24, 30, 40, 50, 60, 70};
  TH1D *rBins[8];
  for(int i=0; i<8; i++) {
    rBins[i] = hist->ProjectionX(Form("bin%d",i+1), ptBins[i]+1, ptBins[i+1]);
    rBins[i]->SetTitle(Form("bin %d: [%d, %d] GeV", i+1, binEdges[i], binEdges[i+1]));
    rBins[i]->SetStats(0);
    rBins[i]->Draw();
    can->SaveAs(Form("plots/ratio_bin%d.pdf", i+1));
    
    double costh_cut = fitf->Integral(binEdges[i], binEdges[i+1])/(binEdges[i+1]-binEdges[i]);
    TLine *cut_m = new TLine(-costh_cut, can->GetUymin(), -costh_cut, can->GetUymax());
    cout << i << " " << can->GetUymax() << endl;
    cut_m->SetLineColor(kRed);
    cut_m->Draw();
    TLine *cut_p = new TLine(costh_cut, can->GetUymin(), costh_cut, can->GetUymax());
    cut_p->SetLineColor(kRed);
    cut_p->Draw();
    
    can->SaveAs(Form("plots/ratio_bin%d.pdf", i+1));
    can->Clear();
  }

  // fine binned
  int nPtBins = hist->GetNbinsY();
  TH1D *fBins[nPtBins];
  for(int i=1; i<=nPtBins; i++) {
    fBins[i-1] = hist->ProjectionX(Form("fine_bin%d",i), i, i+1);
    fBins[i-1]->Draw();
    fBins[i-1]->SetTitle(Form("bin %d: [%.1f, %.1f] GeV", i, hist->GetYaxis()->GetBinLowEdge(i), hist->GetYaxis()->GetBinUpEdge(i)));
    fBins[i-1]->SetStats(0);
    can->SaveAs(Form("plots/ratio_fine_bin%d.pdf", i));

    double costh_cut = fitf->Integral(hist->GetYaxis()->GetBinLowEdge(i), hist->GetYaxis()->GetBinUpEdge(i))/(hist->GetYaxis()->GetBinUpEdge(i) - hist->GetYaxis()->GetBinLowEdge(i));
    TLine *cut_m = new TLine(-costh_cut, can->GetUymin(), -costh_cut, can->GetUymax());
    cut_m->SetLineColor(kRed);
    cut_m->Draw();
    TLine *cut_p = new TLine(costh_cut, can->GetUymin(), costh_cut, can->GetUymax());
    cut_p->SetLineColor(kRed);
    cut_p->Draw();

    costh_cut = costhLim[i-1];
    TLine *cut_mH = new TLine(-costh_cut, can->GetUymin(), -costh_cut, can->GetUymax());
    cut_mH->SetLineColor(kGreen);
    cut_mH->Draw();
    TLine *cut_pH = new TLine(costh_cut, can->GetUymin(), costh_cut, can->GetUymax());
    cut_pH->SetLineColor(kGreen);
    cut_pH->Draw();
      
    can->SaveAs(Form("plots/ratio_fine_bin%d.pdf", i));
    can->Clear();
  }
  */
  
}
