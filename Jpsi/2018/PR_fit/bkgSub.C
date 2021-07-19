#import "../cosMax/imp_jumpF.C"

// macro to subtract background for pol correction

// macro for rounding to integers
int do_round(double val)
{
  int valR = (int)val;
  if (val-valR > 0.5) return valR+1;
  else return valR;
}

void bkgSub()
{
  // PART 1 - all the inputs
  // PR SR data distribution over all pT bins
  TH2D *h_PR2d = new TH2D(); // base PR SR 2d map
  TH2D *h_NP2d = new TH2D(); // base NP SR 2d map
  TH2D *h_MC2d = new TH2D(); // MC 2d map
  TFile *inHist = new TFile("files/histoStore.root");
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
  TH1D **h_SBr = new TH1D*[nBinsY]; // SB background 1d histos
  TFile *inBkg = new TFile("files/bkgCosModel.root");
  for(int i = 0; i < nBinsY; i++) {
    inBkg->GetObject(Form("h_SB_%d", i), h_SBr[i]);
    h_SBr[i]->SetDirectory(0);
  }
  inBkg->Close();
  
  // get the fit range from our cosmax(pT)
  ifstream in;
  string dataS;
  in.open("../cosMax/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", yBins[0]-10, yBins[nBinsY]+10);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);
  
  // bkg fraction in PR SR
  TF1 *f_fSB = new TF1();
  TFile *inFracS = new TFile("files/bkgFrac.root");
  inFracS->GetObject("fit_SB", f_fSB);
  inFracS->Close();

  TGraphErrors *g_fNP = new TGraphErrors();
  TFile *inFracN = new TFile("files/ltfit.root");
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
    
    // get the base data and MC 1d projections
    TH1D *h_PR = h_PR2d->ProjectionX(Form("h_PRSR_%d", i), i+1, i+1);
    TH1D *h_NP = h_NP2d->ProjectionX(Form("h_NP_%d", i), i+1, i+1);
    TH1D *h_MC = h_MC2d->ProjectionX(Form("h_MC_%d", i), i+1, i+1);

    // scale NP dist to unity integral;
    h_NP->Scale(1. / h_NP->Integral());
    
    // multiply analytical function by MC - SB
    TH1D *h_SB = new TH1D();
    h_SB = (TH1D*)h_SBr[i]->Clone(Form("h_PRSB_%d", i));
    h_SB->Sumw2();
    h_SB->Multiply(h_MC);
    // scale SB dist to unity integral;
    h_SB->Scale(1. / h_SB->Integral());

    // PART 3 - scaling background
    // get the proper scaling factor out of the f_bkg
    double f_NP = g_fNP->GetY()[i]/100.;
    double scFac = f_NP * N_sig;
    // scale the background dist
    h_NP->Scale(scFac);

    // get the proper scaling factor out of the f_bkg
    double f_SB = f_fSB->Eval(pt_avg);
    scFac = f_SB * N_sig;
    // scale the background dist
    h_SB->Scale(scFac);

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
  h_Datas->Sumw2();
  h_Datas->Divide(h_MC2d);
  h_Datas->Write();

  h_NPs->Sumw2();
  h_NPs->Divide(h_MC2d);
  h_NPs->Write();

  h_SBs->Sumw2();
  h_SBs->Divide(h_MC2d);
  h_SBs->Write();

  h_PRs->Sumw2();
  h_PRs->Divide(h_MC2d);
  h_PRs->Write();

  h_Js->Sumw2();
  h_Js->Divide(h_MC2d);
  h_Js->Write();
  
  fout->Close();

}
