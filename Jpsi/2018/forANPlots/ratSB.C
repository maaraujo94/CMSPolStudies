// macro to plot the ratios between SB costh in individual plots and pt integrated dist

void ratSB()
{
  // get the dists in each pt bin
  TFile *fin = new TFile("../PR_fit/files/bkgHist.root");
  TH2D *LSBHist2 = (TH2D*)fin->Get("dataH0_ab");
  LSBHist2->SetDirectory(0);
  TH2D *RSBHist2 = (TH2D*)fin->Get("dataH1_ab");
  RSBHist2->SetDirectory(0);
  fin->Close();

  // get the binning
  int nBinsX = LSBHist2->GetNbinsX(), nBinsY = LSBHist2->GetNbinsY();
  const double *yBins = LSBHist2->GetYaxis()->GetXbins()->GetArray();
  double minX = LSBHist2->GetXaxis()->GetBinLowEdge(1);
  double maxX = LSBHist2->GetXaxis()->GetBinUpEdge(nBinsX);

  // define the 1d histos
  TH1D **LSBHist = new TH1D*[nBinsY+1];
  TH1D **RSBHist = new TH1D*[nBinsY+1];
  for(int i = 0; i < nBinsY; i++) {
    LSBHist[i] = LSBHist2->ProjectionX(Form("LSB_%d", i), i+1, i+1);
    LSBHist[i]->SetTitle(Form("LSB Data (%.0f < p_{T} < %.0f GeV)", yBins[i], yBins[i+1]));

    RSBHist[i] = RSBHist2->ProjectionX(Form("RSB_%d", i), i+1, i+1);
    RSBHist[i]->SetTitle(Form("RSB Data (%.0f < p_{T} < %.0f GeV)", yBins[i], yBins[i+1]));
  }
  LSBHist[nBinsY] = new TH1D(Form("LSB_full"), Form("LSB Data (%.0f < p_{T} < %.0f GeV)", yBins[0], yBins[nBinsY]), nBinsX, minX, maxX);
  RSBHist[nBinsY] = new TH1D(Form("RSB_full"), Form("RSB Data (%.0f < p_{T} < %.0f GeV)", yBins[0], yBins[nBinsY]), nBinsX, minX, maxX);

  // open files and read TTrees - fill integrated dist
  TFile *finD = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/data18_cos.root");
  TTree *treeD = (TTree*)finD->Get("data_cos");
  
  int dEvt = treeD->GetEntries();

  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  // cycle over data, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(abs(data_lt) < 0.01 && data_pt > yBins[0] && data_pt < yBins[nBinsY]) { // lt, pT cut
	// LSB - mass cut
	if(data_m > 2.92 && data_m < 2.95) {
	  LSBHist[nBinsY]->Fill(abs(cos(data_th)));
	}
	// RSB - mass cut
	if(data_m > 3.21 && data_m < 3.28) {
	  RSBHist[nBinsY]->Fill(abs(cos(data_th)));
	}
      }
    }
  
  
  // plot costh histograms  
  TH1D **LSBRat = new TH1D*[nBinsY];
  TH1D **RSBRat = new TH1D*[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    LSBRat[i] = (TH1D*)LSBHist[i]->Clone(Form("LSBR_%d", i));
    LSBRat[i]->Sumw2();
    double scF = LSBHist[nBinsY]->Integral(1, nBinsX/2)/LSBHist[i]->Integral(1, nBinsX/2);
    LSBRat[i]->Scale(scF);
    LSBRat[i]->Divide(LSBHist[nBinsY]);
    
    RSBRat[i] = (TH1D*)RSBHist[i]->Clone(Form("RSBR_%d", i));
    scF = RSBHist[nBinsY]->Integral(1, nBinsX/2)/RSBHist[i]->Integral(1, nBinsX/2);
    RSBRat[i]->Scale(scF);
    RSBRat[i]->Sumw2();
    RSBRat[i]->Divide(RSBHist[nBinsY]);    
  }    

  TCanvas *c = new TCanvas("", "", 700, 700);

  // plot the ratios
  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);

  for(int i = 0; i < nBinsY; i+=2) {
    LSBRat[i]->SetStats(0);
    LSBRat[i]->SetLineColor(i/2+1);
    LSBRat[i]->SetMarkerStyle(20);
    LSBRat[i]->SetMarkerColor(i/2+1);
    LSBRat[i]->SetMinimum(0);
    LSBRat[i]->SetTitle("LSB data ratios");
    LSBRat[i]->SetMaximum(2.);
    if(i==0) LSBRat[i]->Draw("error");
    else LSBRat[i]->Draw("error same");

    leg->AddEntry(LSBRat[i], Form("Bin %d: [%.0f, %.0f]", i+1, yBins[i], yBins[i+1]), "pl");
  }
  leg->Draw();

  TLine *zero = new TLine(0, 1, 1, 1);
  zero->SetLineColor(kViolet);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  
  c->SaveAs("plots/LSB_rat.pdf");
  c->Clear();

  for(int i = 0; i < nBinsY; i+=2) {
    RSBRat[i]->SetStats(0);
    RSBRat[i]->SetLineColor(i/2+1);
    RSBRat[i]->SetMarkerStyle(20);
    RSBRat[i]->SetMarkerColor(i/2+1);
    RSBRat[i]->SetMinimum(0);
    RSBRat[i]->SetTitle("RSB data ratios");
    RSBRat[i]->SetMaximum(2.);
    if(i==0) RSBRat[i]->Draw("error");
    else RSBRat[i]->Draw("error same");

  }
  leg->Draw();
  zero->Draw();

  c->SaveAs("plots/RSB_rat.pdf");
  c->Clear();

  // plot the histos
  for(int i = 0; i < nBinsY; i+=2) {
    // scale histo to match full
    LSBHist[i]->Scale(LSBHist[nBinsY]->Integral(1, nBinsX/2)/LSBHist[i]->Integral(1, nBinsX/2));
    LSBHist[i]->SetStats(0);
    LSBHist[i]->SetLineColor(i/2+1);
    LSBHist[i]->SetMarkerStyle(20);
    LSBHist[i]->SetMarkerColor(i/2+1);
    LSBHist[i]->SetMinimum(0);
    LSBHist[i]->SetTitle("LSB data histos");
    LSBHist[i]->SetMaximum(LSBHist[i]->GetMaximum()*2.);
    if(i==0) LSBHist[i]->Draw("error");
    else LSBHist[i]->Draw("error same");

  }
  LSBHist[nBinsY]->SetStats(0);
  LSBHist[nBinsY]->SetLineColor(6);
  LSBHist[nBinsY]->SetMarkerStyle(20);
  LSBHist[nBinsY]->SetMarkerColor(6);
  LSBHist[nBinsY]->Draw("error same");
  leg->AddEntry(LSBHist[nBinsY], Form("p_{T}-integrated"), "pl");

  leg->Draw();

  c->SaveAs("plots/LSB_hist.pdf");
  c->Clear();

  for(int i = 0; i < nBinsY; i+=2) {
    // scale histo to match full
    RSBHist[i]->Scale(RSBHist[nBinsY]->Integral(1, nBinsX/2)/RSBHist[i]->Integral(1, nBinsX/2));
    RSBHist[i]->SetStats(0);
    RSBHist[i]->SetLineColor(i/2+1);
    RSBHist[i]->SetMarkerColor(i/2+1);
    RSBHist[i]->SetMarkerStyle(20);
    RSBHist[i]->SetMinimum(0);
    RSBHist[i]->SetTitle("RSB data histos");
    RSBHist[i]->SetMaximum(RSBHist[i]->GetMaximum()*2.);
    if(i==0) RSBHist[i]->Draw("error");
    else RSBHist[i]->Draw("error same");

  }
  RSBHist[nBinsY]->SetStats(0);
  RSBHist[nBinsY]->SetLineColor(6);
  RSBHist[nBinsY]->SetMarkerStyle(20);
  RSBHist[nBinsY]->SetMarkerColor(6);
  RSBHist[nBinsY]->Draw("error same");

  leg->Draw();

  c->SaveAs("plots/RSB_hist.pdf");
  c->Clear();

}
