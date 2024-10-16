// macro to subtract background for pol correction

#import "../ptbins.C"

void bkgSub()
{
  // PART 1 - all the inputs
  // PR SR data and MC distribution over all pT bins
  TH2D **h_NP2d = new TH2D*[nPtBins]; // base NP SR 2d map
  TH2D **h_MC2d = new TH2D*[nPtBins]; // MC 2d map
  TFile *inHist = new TFile("../PR_fit/files/histoStore.root");
  for(int i = 0; i < nPtBins; i++) {
    inHist->GetObject(Form("NPH_%d", i), h_NP2d[i]);
    h_NP2d[i]->SetDirectory(0);
    inHist->GetObject(Form("MCH_%d", i), h_MC2d[i]);
    h_MC2d[i]->SetDirectory(0);
  }
  inHist->Close();

  // get the binning
  int nBinsX = h_NP2d[0]->GetNbinsX(), nBinsY = h_NP2d[0]->GetNbinsY();
  double minX = h_NP2d[0]->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_NP2d[0]->GetXaxis()->GetBinUpEdge(nBinsX);
  double minY = h_NP2d[0]->GetYaxis()->GetBinLowEdge(1);
  double maxY = h_NP2d[0]->GetYaxis()->GetBinUpEdge(nBinsY);

  // get the bkg distributions
  TH2D **h_SB = new TH2D*[nPtBins]; // SB background 1d histos
  TFile *inBkg = new TFile("files/bkgCosModel.root");
  for(int i = 0; i < nPtBins; i++) {
    inBkg->GetObject(Form("h_SB_%d", i), h_SB[i]);
    h_SB[i]->SetDirectory(0);
  }
  inBkg->Close();
    
  // bkg fraction in PR SR -> 2d cos vs pt object
  // must lose unc and just store value vs pT
  TH2D *h_fb2d = new TH2D();
  TFile *inFracSB = new TFile("../../Simult/bkgFits/files/bkgFrac_NP.root");
  inFracSB->GetObject("h_fbkg", h_fb2d);
  h_fb2d->SetDirectory(0);
  inFracSB->Close();

  // prepare output
  TFile *fout = new TFile("files/bkgSubRes.root", "recreate");
  TCanvas *c =  new TCanvas("", "", 900, 900);
  
  // histos with all the 2d dists to be stored
  TH2D **h_NPs = new TH2D*[nPtBins];
  TH2D **h_SBs = new TH2D*[nPtBins];
  TH2D **h_NPcs = new TH2D*[nPtBins];

  TH2D **h_NPB = new TH2D*[nPtBins];
  TH2D **h_SBB = new TH2D*[nPtBins];
  TH2D **h_NPcB = new TH2D*[nPtBins];
  
  // this part is done for every pT bin
  for(int i = 0; i < nPtBins; i++) {
    h_NPs[i] = new TH2D(Form("h_NPs_%d", i), "NP/MC", nBinsX, minX, maxX, nBinsY, minY, maxY);
    h_SBs[i] = new TH2D(Form("h_SBs_%d", i), "SB/MC", nBinsX, minX, maxX, nBinsY, minY, maxY);
    h_NPcs[i] = new TH2D(Form("h_NPcs_%d", i), "Pure NP/MC", nBinsX, minX, maxX, nBinsY, minY, maxY);
    
    double pt_min = ptBins[i], pt_max = ptBins[i+1];
    double pt_avg = 0.5 * (pt_max + pt_min);
    
    // get number of data events in PR SR
    double N_sig = h_NP2d[i]->Integral(1, nBinsX, 1, nBinsY);
        
    // scale background dists to unity integral;
    cout << h_SB[i]->Integral() << endl;
    h_SB[i]->Scale(1. / h_SB[i]->Integral());

    // PART 3 - scaling background
    double fbg_pt = h_fb2d->GetBinContent(1, i+1);
    cout << fbg_pt << endl;
    h_SB[i]->Scale(fbg_pt); // no propagating unc
    h_SB[i]->Scale(N_sig);

    // PART 4 - signal extraction
    // define the pure NP histo
    TH2D *h_NPc = new TH2D(Form("h_NPc_%d", i), "pure NP cos#theta", nBinsX, minX, maxX, nBinsY, minY, maxY);
    // subtract the background dist from the data dist
    h_NPc->Sumw2();
    h_NPc->Add(h_NP2d[i], h_SB[i], 1, -1); // NP part

    // PART 5 - storing
    // fill the 2d histos
    for(int iX = 1; iX < nBinsX+1; iX++) {
      for(int iY = 1; iY < nBinsY+1; iY++) {
	h_NPs[i]->SetBinContent(iX, iY, h_NP2d[i]->GetBinContent(iX, iY));
	h_NPs[i]->SetBinError(iX, iY, h_NP2d[i]->GetBinError(iX, iY));
	
	h_SBs[i]->SetBinContent(iX, iY, h_SB[i]->GetBinContent(iX, iY));
	h_SBs[i]->SetBinError(iX, iY, h_SB[i]->GetBinError(iX, iY));
	
	h_NPcs[i]->SetBinContent(iX, iY, h_NPc->GetBinContent(iX, iY));
	h_NPcs[i]->SetBinError(iX, iY, h_NPc->GetBinError(iX, iY));
      }
    }

    // PART 6 - output
    h_NPB[i] = (TH2D*)h_NPs[i]->Clone(Form("h_NPB_%d", i));
    h_NPB[i]->Write();
    h_NPs[i]->Sumw2();
    h_NPs[i]->Divide(h_MC2d[i]);
    h_NPs[i]->Write();

    h_SBB[i] = (TH2D*)h_SBs[i]->Clone(Form("h_SBB_%d", i));
    h_SBB[i]->Write();
    h_SBs[i]->Sumw2();
    h_SBs[i]->Divide(h_MC2d[i]);
    h_SBs[i]->Write();

    h_NPcB[i] = (TH2D*)h_NPcs[i]->Clone(Form("h_NPcB_%d", i));
    h_NPcB[i]->Write();
    h_NPcs[i]->Sumw2();
    h_NPcs[i]->Divide(h_MC2d[i]);
    h_NPcs[i]->Write();

  }
  c->Destructor();
  
  
  fout->Close();

}
