void plotCos()
{
  // read the coarse histos in |costh|
  TFile *infile = new TFile("histoStore.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("ratioH_ab"), hist);
  hist->SetDirectory(0);
  infile->Close();

  int nBinsX = hist->GetNbinsX(), nBinsY = hist->GetNbinsY();
  const double *yBins = hist->GetYaxis()->GetXbins()->GetArray();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);

  TH2D *h_inv = new TH2D("h_inv", "data/MC", nBinsY, yBins, nBinsX, minX, maxX);
  for(int i = 0; i < nBinsX; i++) {
    for(int j = 0; j < nBinsY; j++) {
      h_inv->SetBinContent(j+1, i+1, hist->GetBinContent(i+1, j+1));
      h_inv->SetBinError(j+1, i+1, hist->GetBinError(i+1, j+1));
    }
  }
  
  // read the two functions
  ifstream in;
  string dataS;
  in.open("cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", yBins[0]-10, yBins[nBinsY]+10);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);

  ifstream in2;
  in.open("cosMinFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double minPar[3];
  in >> minPar[0] >> aux >> minPar[1] >> aux >> minPar[2] >> aux >> aux;
  in.close();
  
  TF1 *cosMin = new TF1("cosMin", "[0]*(1-exp([1]+[2]*x))", aux, yBins[nBinsY]);
  cosMin->SetParameters(minPar[0], minPar[1], minPar[2]);

  TCanvas *c = new TCanvas("", "", 900, 900);

  double aux1 = h_inv->GetMaximum();
  double aux2 = h_inv->GetMaximum(aux1);
  h_inv->SetMaximum(aux2);
  h_inv->SetStats(0);
  h_inv->GetXaxis()->SetTitle("p_{T} (GeV)");
  h_inv->GetYaxis()->SetTitle("|cos#theta|");
  h_inv->Draw("colz");
  cosMin->Draw("lsame");
  cosMax->Draw("lsame");

  c->SaveAs("plotCos.pdf");
  c->Clear();
  c->Destructor();
}
