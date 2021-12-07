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
  TFile *inPR = new TFile("../../PR_fit/files/bkgHist.root");
  inPR->GetObject("ratioH0_ab", h_pr);
  h_pr->SetDirectory(0);

  int nBinsX = h_pr->GetNbinsX(), nBinsY = h_pr->GetNbinsY();
  const double *yBins = h_pr->GetYaxis()->GetXbins()->GetArray();
  double minX = h_pr->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_pr->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/(double)nBinsX;

  inPR->Close();

  // define aux vals for plotting
  double m_min[] = {3.4, 3.57, 3.82};
  double m_max[] = {3.52, 3.81, 4.0};

  // now get the mass background fit function
  TF1 *fMass = new TF1("fMass", "exp(-x/[0])", 3.35, 4.0);
  // define same function as above but *m
  TF1 *mMass = new TF1("mMass", "exp(-x/[0])*x", 3.35, 4.0);
  // get fMass parameters
  TFile *inFMass = new TFile("files/mfit.root");
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
    fMass->SetParameter(0, m_ld->GetY()[i]);
    mMass->SetParameter(0, m_ld->GetY()[i]);
    
    // get avg mass in sidebands, get f_LSB
    avg_LSB[i] = mMass->Integral(m_min[0], m_max[0])/fMass->Integral(m_min[0], m_max[0]);
    avg_RSB[i] = mMass->Integral(m_min[2], m_max[2])/fMass->Integral(m_min[2], m_max[2]);
    double avg_sig = (m_max[1]+m_min[1])/2.;

    fL[i] = (avg_sig - avg_LSB[i]) / (avg_RSB[i] - avg_LSB[i]);
    
    fL[i] *= 100;

    cout << i << " " << fL[i] << endl;
  }

  TGraphErrors *g_mL = new TGraphErrors(nBinsY, pt_v, avg_LSB, pt_e, zeros); 
  TGraphErrors *g_mR = new TGraphErrors(nBinsY, pt_v, avg_RSB, pt_e, zeros); 
  TGraphErrors *g_fL = new TGraphErrors(nBinsY, pt_v, fL, pt_e, zeros);
  
  TCanvas *c2 =  new TCanvas("", "", 600, 800);
  
  c2->Divide(1, 3);
  c2->cd(1);
  
  TH1F *fmL = gPad->DrawFrame(20, m_min[0], 125, m_max[0]);
  fmL->SetXTitle("p_{T} (GeV)");
  fmL->SetYTitle("<m_{LSB}> (GeV)");
  fmL->GetYaxis()->SetTitleOffset(1.3);
  fmL->GetYaxis()->SetLabelOffset(0.01);
  fmL->SetTitle("2018 <m_{LSB}> (p_{T})");

  g_mL->SetLineColor(kBlack);
  g_mL->SetMarkerStyle(20);
  g_mL->SetMarkerSize(.75);
  g_mL->SetMarkerColor(kBlack);
  g_mL->Draw("psame");

  c2->cd(2);
  
  TH1F *fmR = gPad->DrawFrame(20, m_min[2], 125, m_max[2]);
  fmR->SetXTitle("p_{T} (GeV)");
  fmR->SetYTitle("<m_{RSB}> (GeV)");
  fmR->GetYaxis()->SetTitleOffset(1.3);
  fmR->GetYaxis()->SetLabelOffset(0.01);
  fmR->SetTitle("2018 <m_{RSB}> (p_{T})");

  g_mR->SetLineColor(kBlue);
  g_mR->SetMarkerStyle(20);
  g_mR->SetMarkerSize(.75);
  g_mR->SetMarkerColor(kBlue);
  g_mR->Draw("psame");

  c2->cd(3);
  
  TH1F *ffL = gPad->DrawFrame(20, 50, 125, 55);
  ffL->SetXTitle("p_{T} (GeV)");
  ffL->SetYTitle("f_{L} (%)");
  ffL->GetYaxis()->SetTitleOffset(1.3);
  ffL->GetYaxis()->SetLabelOffset(0.01);
  ffL->SetTitle("2018 f_{L} (p_{T})");

  g_fL->SetLineColor(kRed);
  g_fL->SetMarkerStyle(20);
  g_fL->SetMarkerSize(.75);
  g_fL->SetMarkerColor(kRed);
  g_fL->Draw("psame");
  
  c2->SaveAs(Form("plots/fL.pdf"));
  c2->Clear();
  c2->Destructor();

  TFile *fout = new TFile("files/store_fL.root", "recreate");

  TF1 *fc = new TF1("fc", "[0]", yBins[0], yBins[nBinsY]);
  fc->SetParameter(0,50);
  g_fL->Fit(fc);
  double f_avg = fc->GetParameter(0)/100.;
  int nm = ceil(-log10(f_avg))+2;	
  f_avg = do_round(f_avg*pow(10, nm))/pow(10, nm);

  double fl_fix[nBinsY];
  for(int i = 0; i < nBinsY; i++) fl_fix[i] = f_avg;  
  TGraphErrors *g_flF = new TGraphErrors(nBinsY, pt_v, fl_fix, pt_e, zeros);
  g_flF->Write("g_fL");
  fout->Close();
}
