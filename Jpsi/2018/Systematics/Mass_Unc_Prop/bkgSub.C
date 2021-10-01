// macro to subtract background for pol correction

#import "../../cosMax/imp_jumpF.C"

void bkgSub()
{
  // PART 1 - all the inputs
  // PR SR data distribution over all pT bins
  TH2D *h_PR2d = new TH2D(); // base PR SR 2d map
  TH2D *h_NP2d = new TH2D(); // base NP SR 2d map
  TH2D *h_MC2d = new TH2D(); // MC 2d map
  TFile *inHist = new TFile("../../PR_fit/files/histoStore.root");
  inHist->GetObject("dataH_ab", h_PR2d);
  h_PR2d->SetDirectory(0);
  inHist->GetObject("NPH_ab", h_NP2d);
  h_NP2d->SetDirectory(0);
  inHist->GetObject("mcH_ab", h_MC2d);
  h_MC2d->SetDirectory(0);
  inHist->Close();
  
  // get the binning
  int nBinsX = h_PR2d->GetNbinsX(), nBinsY = h_PR2d->GetNbinsY();
  const double *yBins = h_PR2d->GetYaxis()->GetXbins()->GetArray();
  double minX = h_PR2d->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_PR2d->GetXaxis()->GetBinUpEdge(nBinsX);

  // get the bkg model functions
  TH2D *h_SB2dr = new TH2D();
  TFile *inBkg = new TFile("files/bkgCosModel.root");
  inBkg->GetObject("h_SB", h_SB2dr);
  h_SB2dr->SetDirectory(0);
  inBkg->Close();
  // get the bkg/MC * MC version right away
  TH2D *h_SB2d = new TH2D();
  h_SB2d = (TH2D*)h_SB2dr->Clone("h_bkg");
  h_SB2d->Sumw2();
  h_SB2d->Multiply(h_MC2d);
  
  // get the fit range from our cosmax(pT)
  ifstream in;
  string dataS;
  in.open("../../cosMax/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", yBins[0]-10, yBins[nBinsY]+10);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);
  
  // bkg fraction in PR SR
  TH2D *h_fb2d = new TH2D();
  TFile *inFracS = new TFile("files/bkgFrac.root");
  inFracS->GetObject("h_fbkg", h_fb2d);
  h_fb2d->SetDirectory(0);
  inFracS->Close();

  // NP fraction - still need to add uncertainties
  TGraphErrors *g_fNP = new TGraphErrors();
  TFile *inFracN = new TFile("../../PR_fit/files/ltfit.root");
  inFracN->GetObject("fit_b_fNP", g_fNP);
  inFracN->Close();
  
  TFile *fout = new TFile("files/bkgSubRes.root", "recreate");
  TCanvas *c =  new TCanvas("", "", 900, 900);
  
  // histos with all the 2d dists to be stored
  TH2D *h_Datas = new TH2D("h_Data", "Peak/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_NPs = new TH2D("h_NP", "NP/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_SBs = new TH2D("h_SB", "SB/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_PRs = new TH2D("h_PR", "Prompt/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_Js = new TH2D("h_J", "Prompt J/#psi/MC", nBinsX, minX, maxX, nBinsY, yBins);

  // this part is done for every pT bin
  for(int i = 0; i < nBinsY; i++) {
    double pt_min = yBins[i], pt_max = yBins[i+1];
    double pt_avg = 0.5 * (pt_max + pt_min);
    
    // get number of data events in PR SR
    double N_sig = h_PR2d->Integral(1, nBinsX, i+1, i+1);
    
    // getting the max costh value for the fit, cR
    double cMaxVal = jumpF(cosMax->Integral(pt_min, pt_max)/(pt_max-pt_min));
    
    // get the data and MC 1d projections
    TH1D *h_PR = h_PR2d->ProjectionX(Form("h_PRSR_%d", i), i+1, i+1);
    TH1D *h_NP = h_NP2d->ProjectionX(Form("h_NP_%d", i), i+1, i+1);
    TH1D *h_MC = h_MC2d->ProjectionX(Form("h_MC_%d", i), i+1, i+1);
    TH1D *h_SB = h_SB2d->ProjectionX(Form("h_SB_%d", i), i+1, i+1);
    // get the fbkg 1d projections - easier to propagate unc
    TH1D *h_fbkg = h_fb2d->ProjectionX(Form("h_fbkg_%d", i), i+1, i+1);
    
    // scale background dists to unity integral;
    h_NP->Scale(1. / h_NP->Integral());
    h_SB->Scale(1. / h_SB->Integral());

    // PART 3 - scaling background
    // get the proper scaling factor out of the f_NP
    double f_NP = g_fNP->GetY()[i]/100.;
    double scFac = f_NP * N_sig;
    h_NP->Scale(scFac);

    // get the proper scaling factor out of the f_bkg
    h_SB->Multiply(h_fbkg); // propagating unc
    h_SB->Scale(N_sig);

    // PART 4 - signal extraction
    // define the pure PR histo
    TH1D *h_justPR = new TH1D(Form("h_justPR_%d", i), "prompt cos#theta", nBinsX, minX, maxX);
    // subtract the background dist from the data dist
    h_justPR->Sumw2();
    h_justPR->Add(h_PR, h_NP, 1, -1); // NP part
    // define the pure PR histo
    TH1D *h_purePR = new TH1D(Form("h_purePR_%d", i), "prompt J/#psi cos#theta", nBinsX, minX, maxX);
    // subtract the background dist from the data dist
    h_purePR->Sumw2();
    h_purePR->Add(h_justPR, h_SB, 1, -1); // sideband part

    // PART 5 - storing

    h_PR->Sumw2();
    
    // fill the 2d histos
    for(int iX = 0; iX < nBinsX; iX++) {
      h_Datas->SetBinContent(iX+1, i+1, h_PR->GetBinContent(iX+1));
      h_Datas->SetBinError(iX+1, i+1, h_PR->GetBinError(iX+1));
      
      h_NPs->SetBinContent(iX+1, i+1, h_NP->GetBinContent(iX+1));
      h_NPs->SetBinError(iX+1, i+1, h_NP->GetBinError(iX+1));
      
      h_SBs->SetBinContent(iX+1, i+1, h_SB->GetBinContent(iX+1));
      h_SBs->SetBinError(iX+1, i+1, h_SB->GetBinError(iX+1));
      
      h_PRs->SetBinContent(iX+1, i+1, h_justPR->GetBinContent(iX+1));
      h_PRs->SetBinError(iX+1, i+1, h_justPR->GetBinError(iX+1));
      
      h_Js->SetBinContent(iX+1, i+1, h_purePR->GetBinContent(iX+1));
      h_Js->SetBinError(iX+1, i+1, h_purePR->GetBinError(iX+1));
    }
  
  }
  c->Destructor();
  
  // PART 6 - output
  TH2D *h_DataB = new TH2D();
  h_DataB = (TH2D*)h_Datas->Clone("h_DataB");
  h_DataB->Write();
  h_Datas->Sumw2();
  h_Datas->Divide(h_MC2d);
  h_Datas->Write();

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

  TH2D *h_PRB = new TH2D();
  h_PRB = (TH2D*)h_PRs->Clone("h_PRB");
  h_PRB->Write();
  h_PRs->Sumw2();
  h_PRs->Divide(h_MC2d);
  h_PRs->Write();

  TH2D *h_JB = new TH2D();
  h_JB = (TH2D*)h_Js->Clone("h_JB");
  h_JB->Write();
  h_Js->Sumw2();
  h_Js->Divide(h_MC2d);
  h_Js->Write();
  
  fout->Close();

}
