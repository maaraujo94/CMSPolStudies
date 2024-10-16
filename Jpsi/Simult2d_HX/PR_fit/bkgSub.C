// macro to subtract background for pol correction

#import "../ptbins.C"

void bkgSub()
{
  // PART 1 - all the inputs
  // PR SR data and MC distribution over all pT bins
  TH2D **h_PR2d = new TH2D*[nPtBins]; // base PR SR 2d map
  TH2D **h_MC2d = new TH2D*[nPtBins]; // MC 2d map
  TFile *inHist = new TFile("files/histoStore.root");
  for(int i = 0; i < nPtBins; i++) {
    inHist->GetObject(Form("PRH_%d", i), h_PR2d[i]);
    h_PR2d[i]->SetDirectory(0);
    inHist->GetObject(Form("MCH_%d", i), h_MC2d[i]);
    h_MC2d[i]->SetDirectory(0);
  }
  inHist->Close();

  // NP distribution corrected for mass bkg contamination
  TH2D **h_NP2d = new TH2D*[nPtBins]; // base NP SR 2d map
  TFile *inNP = new TFile("../NP_fit/files/bkgSubRes.root");
  for(int i = 0; i < nPtBins; i++) {
    inNP->GetObject(Form("h_NPcB_%d", i), h_NP2d[i]);
    h_NP2d[i]->SetDirectory(0);
  }
  inNP->Close();
  
  // get the binning
  int nBinsX = h_PR2d[0]->GetNbinsX(), nBinsY = h_PR2d[0]->GetNbinsY();
  double minX = h_PR2d[0]->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_PR2d[0]->GetXaxis()->GetBinUpEdge(nBinsX);
  double minY = h_PR2d[0]->GetYaxis()->GetBinLowEdge(1);
  double maxY = h_PR2d[0]->GetYaxis()->GetBinUpEdge(nBinsY);

  // get the bkg distributions
  TH2D **h_SB = new TH2D*[nBinsY]; // SB background 1d histos
  TFile *inBkg = new TFile("files/bkgCosModel.root");
  for(int i = 0; i < nPtBins; i++) {
    inBkg->GetObject(Form("h_SB_%d", i), h_SB[i]);
    h_SB[i]->SetDirectory(0);
  }
  inBkg->Close();
    
  // bkg fraction in PR SR
  // NP fraction in NP SR - corrected for mass bkg contamination
  TH2D *h_fb2d = new TH2D();
  TH2D *h_fnp2d = new TH2D();
  TFile *inFracSB = new TFile("../../Simult/bkgFits/files/bkgFrac.root");
  inFracSB->GetObject("h_fbkg", h_fb2d);
  h_fb2d->SetDirectory(0);
  inFracSB->Close();
  TFile *inFracNP = new TFile("../../Simult/PR_fit/files/NPFrac.root");
  inFracNP->GetObject("h_fnp", h_fnp2d);
  h_fnp2d->SetDirectory(0);
  inFracNP->Close();

  // prepare output
  TFile *fout = new TFile("files/bkgSubRes.root", "recreate");
  TCanvas *c =  new TCanvas("", "", 900, 900);
  
  // histos with all the 2d dists to be stored
  TH2D **h_Datas = new TH2D*[nPtBins];
  TH2D **h_NPs = new TH2D*[nPtBins];
  TH2D **h_SBs = new TH2D*[nPtBins];
  TH2D **h_PRs = new TH2D*[nPtBins];
  TH2D **h_Js = new TH2D*[nPtBins];

  TH2D **h_DataB = new TH2D*[nPtBins];
  TH2D **h_NPB = new TH2D*[nPtBins];
  TH2D **h_SBB = new TH2D*[nPtBins];
  TH2D **h_PRB = new TH2D*[nPtBins];
  TH2D **h_JB = new TH2D*[nPtBins];
 
  // this part is done for every pT bin
  for(int i = 0; i < nPtBins; i++) {
    h_Datas[i] = new TH2D(Form("h_Data_%d", i), "Peak/MC", nBinsX, minX, maxX, nBinsY, minY, maxY);
    h_NPs[i] = new TH2D(Form("h_NP_%d", i), "NP/MC", nBinsX, minX, maxX, nBinsY, minY, maxY);
    h_SBs[i] = new TH2D(Form("h_SB_%d", i), "SB/MC", nBinsX, minX, maxX, nBinsY, minY, maxY);
    h_PRs[i] = new TH2D(Form("h_PR_%d", i), "Prompt/MC", nBinsX, minX, maxX, nBinsY, minY, maxY);
    h_Js[i] = new TH2D(Form("h_J_%d", i), "Prompt J/#psi/MC", nBinsX, minX, maxX, nBinsY, minY, maxY);

    double pt_min = ptBins[i], pt_max = ptBins[i+1];
    double pt_avg = 0.5 * (pt_max + pt_min);
    
    // get number of data events in PR SR
    double N_sig = h_PR2d[i]->Integral(1, nBinsX, 1, nBinsY);
        
    // scale background dists to unity integral;
    h_NP2d[i]->Scale(1. / h_NP2d[i]->Integral());
    h_SB[i]->Scale(1. / h_SB[i]->Integral());

    // PART 3 - scaling background
    double fnp_pt = h_fnp2d->GetBinContent(1, i+1);
    double fbg_pt = h_fb2d->GetBinContent(1, i+1);
    h_NP2d[i]->Scale(fnp_pt); // propagating unc
    h_NP2d[i]->Scale(N_sig);
    h_SB[i]->Scale(fbg_pt); // propagating unc
    h_SB[i]->Scale(N_sig);

    // PART 4 - signal extraction
    // define the pure PR histo
    TH2D *h_justPR = new TH2D(Form("h_justPR_%d", i), "prompt cos#theta", nBinsX, minX, maxX, nBinsY, minY, maxY);
    // subtract the background dist from the data dist
    h_justPR->Sumw2();
    h_justPR->Add(h_PR2d[i], h_NP2d[i], 1, -1); // NP part
    // define the pure PR histo
    TH2D *h_purePR = new TH2D(Form("h_purePR_%d", i), "prompt J/#psi cos#theta", nBinsX, minX, maxX, nBinsY, minY, maxY);
    // subtract the background dist from the data dist
    h_purePR->Sumw2();
    h_purePR->Add(h_justPR, h_SB[i], 1, -1); // sideband part

    // PART 5 - storing

    h_PR2d[i]->Sumw2();
    
    // fill the 2d histos
    for(int iX = 1; iX < nBinsX+1; iX++) {
      for(int iY = 1; iY < nBinsY+1; iY++) {
	h_Datas[i]->SetBinContent(iX, iY, h_PR2d[i]->GetBinContent(iX, iY));
	h_Datas[i]->SetBinError(iX, iY, h_PR2d[i]->GetBinError(iX, iY));
      
	h_NPs[i]->SetBinContent(iX, iY, h_NP2d[i]->GetBinContent(iX, iY));
	h_NPs[i]->SetBinError(iX, iY, h_NP2d[i]->GetBinError(iX, iY));
      
	h_SBs[i]->SetBinContent(iX, iY, h_SB[i]->GetBinContent(iX, iY));
	h_SBs[i]->SetBinError(iX, iY, h_SB[i]->GetBinError(iX, iY));
      
	h_PRs[i]->SetBinContent(iX, iY, h_justPR->GetBinContent(iX, iY));
	h_PRs[i]->SetBinError(iX, iY, h_justPR->GetBinError(iX, iY));
      
	h_Js[i]->SetBinContent(iX, iY, h_purePR->GetBinContent(iX, iY));
	h_Js[i]->SetBinError(iX, iY, h_purePR->GetBinError(iX, iY));
      }
    }
    // PART 6 - output
    h_DataB[i] = (TH2D*)h_Datas[i]->Clone(Form("h_DataB_%d", i));
    h_DataB[i]->Write();
    h_Datas[i]->Sumw2();
    h_Datas[i]->Divide(h_MC2d[i]);
    h_Datas[i]->Write();

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

    h_PRB[i] = (TH2D*)h_PRs[i]->Clone(Form("h_PRB_%d", i));
    h_PRB[i]->Write();
    h_PRs[i]->Sumw2();
    h_PRs[i]->Divide(h_MC2d[i]);
    h_PRs[i]->Write();

    h_JB[i] = (TH2D*)h_Js[i]->Clone(Form("h_JB_%d", i));
    h_JB[i]->Write();
    h_Js[i]->Sumw2();
    h_Js[i]->Divide(h_MC2d[i]);
    h_Js[i]->Write();
  
  }
  c->Destructor();
  
  
  fout->Close();

}
