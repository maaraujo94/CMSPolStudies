// macro to get the fraction of LSB contributing to SR (fL)

int do_round(double val)
{
  int valR = (int)val;
  if (val-valR > 0.5) return valR+1;
  else return valR;
}

void getfL()
{
  // get the binning
  TH2D *h_pr = new TH2D();
  TFile *inPR = new TFile("/home/mariana/Documents/2021_PhD_work/CERN/0618_plotMass/02_fitCosth/bkgHist.root");
  inPR->GetObject("ratioH0_ab", h_pr);
  h_pr->SetDirectory(0);

  int nBinsX = h_pr->GetNbinsX(), nBinsY = h_pr->GetNbinsY();
  const double *yBins = h_pr->GetYaxis()->GetXbins()->GetArray();
  double minX = h_pr->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_pr->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/(double)nBinsX;

  inPR->Close();

  // now get the mass background fit function
  TF1 *fMass = new TF1("fMass", "exp(-x/[0])", 2.9, 3.3);
  // define same function as above but *m
  TF1 *mMass = new TF1("mMass", "exp(-x/[0])*x", 2.9, 3.3);
  // get fMass parameters
  TFile *inFMass = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/2018/bkgFits/files/mfit_data.root");
  TGraphErrors *m_ld = (TGraphErrors*)inFMass->Get("fit_lambda");
  inFMass->Close();
  
  // this part is done for every pT bin
  double avg_LSB[nBinsY], avg_RSB[nBinsY], fL[nBinsY];
  double pt_v[nBinsY], pt_e[nBinsY], zeros[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    pt_v[i] = 0.5*(yBins[i+1]+yBins[i]);
    pt_e[i] = 0.5*(yBins[i+1]-yBins[i]);
    zeros[i] = 1e-4;
    
    // set the parameters of the binned functions
    fMass->SetParameter(0, m_ld->GetY()[i]*1e-3);
    mMass->SetParameter(0, m_ld->GetY()[i]*1e-3);
    
    // get avg mass in sidebands, get f_LSB
    avg_LSB[i] = mMass->Integral(2.92, 2.95)/fMass->Integral(2.92, 2.95);
    avg_RSB[i] = mMass->Integral(3.21, 3.28)/fMass->Integral(3.21, 3.28);
    fL[i] = (3.1 - avg_LSB[i]) / (avg_RSB[i] - avg_LSB[i]);
    
    fL[i] *= 100;
  }

  
  TGraphErrors *g_mL = new TGraphErrors(nBinsY, pt_v, avg_LSB, pt_e, zeros); 
  TGraphErrors *g_mR = new TGraphErrors(nBinsY, pt_v, avg_RSB, pt_e, zeros); 
  TGraphErrors *g_fL = new TGraphErrors(nBinsY, pt_v, fL, pt_e, zeros);

  
  TCanvas *c2 =  new TCanvas("", "", 600, 800);

  
  c2->Divide(1, 3);
  c2->cd(1);
  
  TH1F *fmL = gPad->DrawFrame(20, 2.92, 105, 2.95);
  fmL->SetXTitle("p_{T} (GeV)");
  fmL->SetYTitle("<m_{LSB}> (GeV)");
  fmL->GetYaxis()->SetTitleOffset(1.3);
  fmL->GetYaxis()->SetLabelOffset(0.01);
  fmL->SetTitle("<m_{LSB}> (p_{T})");

  g_mL->SetLineColor(kBlack);
  g_mL->SetMarkerStyle(20);
  g_mL->SetMarkerSize(.75);
  g_mL->SetMarkerColor(kBlack);
  g_mL->Draw("psame");

  c2->cd(2);
  
  TH1F *fmR = gPad->DrawFrame(20, 3.21, 105, 3.28);
  fmR->SetXTitle("p_{T} (GeV)");
  fmR->SetYTitle("<m_{RSB}> (GeV)");
  fmR->GetYaxis()->SetTitleOffset(1.3);
  fmR->GetYaxis()->SetLabelOffset(0.01);
  fmR->SetTitle("<m_{RSB}> (p_{T})");

  g_mR->SetLineColor(kBlue);
  g_mR->SetMarkerStyle(20);
  g_mR->SetMarkerSize(.75);
  g_mR->SetMarkerColor(kBlue);
  g_mR->Draw("psame");

  c2->cd(3);
  
  TH1F *ffL = gPad->DrawFrame(20, 50, 105, 55);
  ffL->SetXTitle("p_{T} (GeV)");
  ffL->SetYTitle("f_{L} (%)");
  ffL->GetYaxis()->SetTitleOffset(1.3);
  ffL->GetYaxis()->SetLabelOffset(0.01);
  ffL->SetTitle("f_{L} (p_{T})");

  g_fL->SetLineColor(kRed);
  g_fL->SetMarkerStyle(20);
  g_fL->SetMarkerSize(.75);
  g_fL->SetMarkerColor(kRed);
  g_fL->Draw("psame");

  TF1 *fcons = new TF1("fcons", "[0]", yBins[0]-5, yBins[nBinsY]+5);
  fcons->SetParameter(0,53.3);
  fcons->SetLineColor(kBlack);
  fcons->SetLineStyle(kDashed);
  //fcons->Draw("same");
  
  c2->SaveAs(Form("plots/fL.pdf"));
  c2->Clear();
  c2->Destructor();
}
