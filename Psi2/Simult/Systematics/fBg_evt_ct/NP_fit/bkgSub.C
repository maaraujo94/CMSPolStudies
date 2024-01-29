// macro to subtract background for pol correction

void bkgSub()
{
  // PART 1 - all the inputs
  // PR SR data distribution over all pT bins
  TH2D *h_NP2d = new TH2D(); // base NP SR 2d map
  TH2D *h_MC2d = new TH2D(); // MC 2d map
  TFile *inHist = new TFile("../../../PR_fit/files/histoStore.root");
  inHist->GetObject("NPH", h_NP2d);
  h_NP2d->SetDirectory(0);
  inHist->GetObject("MCH", h_MC2d);
  h_MC2d->SetDirectory(0);
  inHist->Close();
  
  // get the binning
  int nBinsX = h_NP2d->GetNbinsX(), nBinsY = h_NP2d->GetNbinsY();
  const double *yBins = h_NP2d->GetYaxis()->GetXbins()->GetArray();
  double minX = h_NP2d->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_NP2d->GetXaxis()->GetBinUpEdge(nBinsX);

  // get the bkg dists
  TH1D **h_SB = new TH1D*[nBinsY]; // SB background 1d histos
  TFile *inBkg = new TFile("files/bkgCosModel.root");
  for(int i = 0; i < nBinsY; i++) {
    inBkg->GetObject(Form("h_SB_%d", i), h_SB[i]);
    h_SB[i]->SetDirectory(0);
  }
  inBkg->Close();
  
  // bkg fraction in NP SR
  TH2D *h_fb2d = new TH2D();
  TFile *inFracSB = new TFile("../bkgFits/files/bkgFrac_NP.root");
  inFracSB->GetObject("h_fbkg", h_fb2d);
  h_fb2d->SetDirectory(0);
  inFracSB->Close();

  TFile *fout = new TFile("files/bkgSubRes.root", "recreate");
  TCanvas *c =  new TCanvas("", "", 900, 900);
  
  // histos with all the 2d dists to be stored
  TH2D *h_NPs = new TH2D("h_NP", "NP/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_SBs = new TH2D("h_SB", "SB/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_NPcs = new TH2D("h_NPc", "Pure NP/MC", nBinsX, minX, maxX, nBinsY, yBins);

  // this part is done for every pT bin
  for(int i = 0; i < nBinsY; i++) {
    double pt_min = yBins[i], pt_max = yBins[i+1];
    double pt_avg = 0.5 * (pt_max + pt_min);
    
    // get number of data events in PR SR
    double N_sig = h_NP2d->Integral(1, nBinsX, i+1, i+1);
    
    // get the base data and MC 1d projections
    TH1D *h_NP = h_NP2d->ProjectionX(Form("h_NP_%d", i), i+1, i+1);
    TH1D *h_MC = h_MC2d->ProjectionX(Form("h_MC_%d", i), i+1, i+1);
    // get the fbkg 1d projections - easier to propagate unc
    TH1D *h_fbkg = h_fb2d->ProjectionX(Form("h_fbkg_%d", i), i+1, i+1);

    // scale SB dist to unity integral;
    h_SB[i]->Scale(1. / h_SB[i]->Integral());

    // get the proper scaling factor out of the f_bkg
    h_SB[i]->Multiply(h_fbkg); // propagating unc
    h_SB[i]->Scale(N_sig);

    // PART 4 - signal extraction
    // define the pure PR histo
    TH1D *h_NPc = new TH1D(Form("h_NPc_%d", i), "pure NP cos#theta", nBinsX, minX, maxX);
    // subtract the background dist from the data dist
    h_NPc->Sumw2();
    h_NPc->Add(h_NP, h_SB[i], 1, -1); // NP part
 
    // PART 5 - storing    
    // fill the 2d histos
    for(int iX = 0; iX < nBinsX; iX++) {
      h_NPs->SetBinContent(iX+1, i+1, h_NP->GetBinContent(iX+1));
      h_NPs->SetBinError(iX+1, i+1, h_NP->GetBinError(iX+1));
      
      h_SBs->SetBinContent(iX+1, i+1, h_SB[i]->GetBinContent(iX+1));
      h_SBs->SetBinError(iX+1, i+1, h_SB[i]->GetBinError(iX+1));
      
      h_NPcs->SetBinContent(iX+1, i+1, h_NPc->GetBinContent(iX+1));
      h_NPcs->SetBinError(iX+1, i+1, h_NPc->GetBinError(iX+1));
    }
  
  }
  c->Destructor();
  
  // PART 6 - output
  TH2D *h_NPB = new TH2D();
  h_NPB = (TH2D*)h_NPs->Clone("h_NPB");
  h_NPB->Write();
  h_NPs->Sumw2();
  h_NPs->Divide(h_MC2d);
  h_NPs->Write();

  TH2D *h_SBB = new TH2D();
  h_SBB = (TH2D*)h_SBs->Clone("h_SBB");
  h_SBB->Write();
  h_SBs->Sumw2();
  h_SBs->Divide(h_MC2d);
  h_SBs->Write();

  TH2D *h_NPcB = new TH2D();
  h_NPcB = (TH2D*)h_NPcs->Clone("h_NPcB");
  h_NPcB->Write();
  h_NPcs->Sumw2();
  h_NPcs->Divide(h_MC2d);
  h_NPcs->Write();

  fout->Close();

}
