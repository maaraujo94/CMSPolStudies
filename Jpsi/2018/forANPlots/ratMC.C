// macro to plot the ratios between MC costh in individual plots and pt integrated dist

int DO_FILL = 0;

void ratMC()
{
  // prepare binning and histograms for plots
  const int nBinsY = 7;
  double yBins[nBinsY+1];
  for(int i=0; i<3; i++) yBins[i] = 7.*i+25.;
  for(int i=0; i<4; i++) yBins[i+3] = 46.+10.*i;
  yBins[7] = 120;
  for(int i=0; i<nBinsY+1; i++) cout << yBins[i] << ",";
  cout << endl;

  double minX = 0, maxX = 1, nBinsX = 20;

  // define the 1d histos
  TH1D **MCHist = new TH1D*[nBinsY+1];
  for(int i = 0; i < nBinsY; i++) {
    MCHist[i] = new TH1D(Form("MC_%d", i), Form("MC (%.0f < p_{T} < %.0f GeV)", yBins[i], yBins[i+1]), nBinsX, minX, maxX);
  }
  MCHist[nBinsY] = new TH1D(Form("MC_full"), Form("MC (%.0f < p_{T} < %.0f GeV)", yBins[0], yBins[nBinsY]), nBinsX, minX, maxX);

  if (DO_FILL == 1) {
    // open files and read TTrees - fill all dists
    TFile *fin2 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MC18_cos.root");
    TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
    TFile *fin3 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCh18_cos.root");
    TTree *treeM2 = (TTree*)fin3->Get("MC_cos");
    TFile *fin4 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MCvh18_cos.root");
    TTree *treeM3 = (TTree*)fin4->Get("MC_cos");
    
    int m1Evt = treeM1->GetEntries();
    int m2Evt = treeM2->GetEntries();
    int m3Evt = treeM3->GetEntries();
    
    // definitions to store data and MC events
    Double_t mc_th, mc_pt, mc_lt, mc_m;
    
    treeM1->SetBranchAddress("theta", &mc_th);
    treeM1->SetBranchAddress("dimPt", &mc_pt);
    treeM1->SetBranchAddress("Mass", &mc_m);
    treeM1->SetBranchAddress("lt", &mc_lt);
  
    // cycle over mc, fill the costh histogram acc to binning
    for(int i = 0; i < m1Evt; i++)
      {
	treeM1->GetEntry(i);
	if(abs(mc_lt) < 0.01 && mc_pt > yBins[0] && mc_pt < 46 && mc_m > 3.0 && mc_m < 3.2) { // lt, pT, mass cut
	  MCHist[nBinsY]->Fill(abs(cos(mc_th))); // fill full histo
	  // fill binned histos
	  for(int i_p = 0; i_p < nBinsY; i_p++)
	    if(mc_pt > yBins[i_p] && mc_pt < yBins[i_p+1]) {
	      MCHist[i_p]->Fill(abs(cos(mc_th)));
	    }
	}
      }

    for(int i = 0; i < nBinsY+1; i++) {
      cout << MCHist[i]->Integral() << " ";
    }
    cout << endl;
  
    treeM2->SetBranchAddress("theta", &mc_th);
    treeM2->SetBranchAddress("dimPt", &mc_pt);
    treeM2->SetBranchAddress("Mass", &mc_m);
    treeM2->SetBranchAddress("lt", &mc_lt);
  
    // cycle over mc, fill the costh histogram acc to binning
    for(int i = 0; i < m2Evt; i++)
      {
	treeM2->GetEntry(i);
	if(abs(mc_lt) < 0.01 && mc_pt > 46 && mc_pt < 66 && mc_m > 3.0 && mc_m < 3.2) { // lt, pT, mass cut
	  MCHist[nBinsY]->Fill(abs(cos(mc_th))); // fill full histo
	  // fill binned histos
	  for(int i_p = 0; i_p < nBinsY; i_p++)
	    if(mc_pt > yBins[i_p] && mc_pt < yBins[i_p+1]) {
	      MCHist[i_p]->Fill(abs(cos(mc_th)));
	    }
	}
      }

    for(int i = 0; i < nBinsY+1; i++) {
      cout << MCHist[i]->Integral() << " ";
    }
    cout << endl;
  
    treeM3->SetBranchAddress("theta", &mc_th);
    treeM3->SetBranchAddress("dimPt", &mc_pt);
    treeM3->SetBranchAddress("Mass", &mc_m);
    treeM3->SetBranchAddress("lt", &mc_lt);
  
    // cycle over mc, fill the costh histogram acc to binning
    for(int i = 0; i < m3Evt; i++)
      {
	treeM3->GetEntry(i);
	if(abs(mc_lt) < 0.01 && mc_pt > 66 && mc_pt < yBins[nBinsY] && mc_m > 3.0 && mc_m < 3.2) { // lt, pT, mass cut
	  MCHist[nBinsY]->Fill(abs(cos(mc_th))); // fill full histo
	  // fill binned histos
	  for(int i_p = 0; i_p < nBinsY; i_p++)
	    if(mc_pt > yBins[i_p] && mc_pt < yBins[i_p+1]) {
	      MCHist[i_p]->Fill(abs(cos(mc_th)));
	    }
	}
      }

    for(int i = 0; i < nBinsY+1; i++) {
      cout << MCHist[i]->Integral() << " ";
    }
    cout << endl;

    // store MC histos
    TFile *fout = new TFile("store_SBcomp.root", "recreate");
    for(int i = 0; i < nBinsY+1; i++)
      MCHist[i]->Write();
    fout->Close();
  }
  else {
    // read MC histos
    TFile *fin = new TFile("store_SBcomp.root");
    for(int i = 0; i < nBinsY; i++) {
      fin->GetObject(Form("MC_%d", i), MCHist[i]);
      MCHist[i]->SetDirectory(0);
    }
    fin->GetObject(Form("MC_full"), MCHist[nBinsY]);
    MCHist[nBinsY]->SetDirectory(0);
    fin->Close();
  }
  
  // plot costh histograms  
  TH1D **MCRat = new TH1D*[nBinsY];
  TH1D **RSBRat = new TH1D*[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    MCRat[i] = (TH1D*)MCHist[i]->Clone(Form("MCR_%d", i));
    MCRat[i]->Sumw2();
    double scF = MCHist[nBinsY]->Integral(1, nBinsX/2)/MCHist[i]->Integral(1, nBinsX/2);
    MCRat[i]->Scale(scF);
    MCRat[i]->Divide(MCHist[nBinsY]);
  }    

  TCanvas *c = new TCanvas("", "", 700, 700);

  // plot the ratios
  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);

  for(int i = 0; i < nBinsY; i+=2) {
    MCRat[i]->SetStats(0);
    MCRat[i]->SetLineColor(i/2+1);
    MCRat[i]->SetMarkerStyle(20);
    MCRat[i]->SetMarkerColor(i/2+1);
    MCRat[i]->SetMinimum(0.75);
    MCRat[i]->SetTitle("MC data ratios");
    MCRat[i]->SetMaximum(1.25);
    if(i==0) MCRat[i]->Draw("error");
    else MCRat[i]->Draw("error same");

    leg->AddEntry(MCRat[i], Form("Bin %d: [%.0f, %.0f]", i+1, yBins[i], yBins[i+1]), "pl");
  }
  leg->Draw();

  TLine *zero = new TLine(0, 1, 1, 1);
  zero->SetLineColor(kViolet);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  
  c->SaveAs("plots/MC_rat.pdf");
  c->Clear();

  // plot the histos
  for(int i = 0; i < nBinsY; i+=2) {
    // scale histo to match full
    MCHist[i]->Scale(MCHist[nBinsY]->Integral(1, nBinsX/2)/MCHist[i]->Integral(1, nBinsX/2));
    MCHist[i]->SetStats(0);
    MCHist[i]->SetLineColor(i/2+1);
    MCHist[i]->SetMarkerStyle(20);
    MCHist[i]->SetMarkerColor(i/2+1);
    MCHist[i]->SetMinimum(0);
    MCHist[i]->SetTitle("MC data histos");
    MCHist[i]->SetMaximum(MCHist[i]->GetMaximum()*2.);
    if(i==0) MCHist[i]->Draw("error");
    else MCHist[i]->Draw("error same");

  }
  MCHist[nBinsY]->SetStats(0);
  MCHist[nBinsY]->SetLineColor(6);
  MCHist[nBinsY]->SetMarkerStyle(20);
  MCHist[nBinsY]->SetMarkerColor(6);
  MCHist[nBinsY]->Draw("error same");
  leg->AddEntry(MCHist[nBinsY], Form("p_{T}-integrated"), "pl");

  leg->Draw();

  c->SaveAs("plots/MC_hist.pdf");
  c->Clear();


}
