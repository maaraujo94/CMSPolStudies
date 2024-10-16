// macro to get the maximum |costh| in each pT bin that we can fit

#import "../ptbins.C"

void getCos()
{
  double aux = 0;

  ofstream fout;
  fout.open("cos_max.txt");
  fout.close();  
  
  // read the 2d histos
  TFile *infile = new TFile("../../Simult2d_CS/PR_fit/files/histoStore.root");
  TH2D **h2d = new TH2D*[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    infile->GetObject(Form("rPRH_%d", i), h2d[i]);
    h2d[i]->SetDirectory(0);
  }
  infile->Close();

  double cosMax[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    int nBinsX = h2d[i]->GetNbinsX();
    for(int j = 1; j < nBinsX; j++) {
      if(h2d[i]->GetBinContent(j,1)!=0) {
	cosMax[i] = h2d[i]->GetXaxis()->GetBinLowEdge(j);
	cout << i << " " << ptBins[i] << " " << ptBins[i+1] << " " << cosMax[i] << endl;
	break;
      }
    }
  }
  
  fout.open("cos_max.txt", std::ofstream::app);
  for(int i = 0; i < nPtBins; i++) {
    fout << ptBins[i] << "\t" << ptBins[i+1] << "\t" << cosMax[i] << endl;
  }
  fout.close();

  // plot the results in a histogram and fit it
  TCanvas *can = new TCanvas("", "", 700, 700);
  
  // the histogram
  TH1F* costh = new TH1F("name", "cos#theta limits", nPtBins, ptBins);
  costh->GetXaxis()->SetTitle("p_{T} (GeV)");
  costh->GetYaxis()->SetTitle("|cos#theta_{CS}|_{max}");
  for(int i = 0; i < nPtBins; i++) {
    costh->SetBinContent(i+1, cosMax[i]);
    costh->SetBinError(i+1, 0.025);
  }
  costh->SetLineColor(kBlack);
  costh->SetStats(0);
  costh->Draw("");
  
  can->SaveAs("costh_lim.pdf");
  can->Clear();
  can->Destructor();  
}
