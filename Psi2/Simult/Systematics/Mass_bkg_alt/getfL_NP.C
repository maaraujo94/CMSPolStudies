// macro to get the fraction of LSB contributing to SR (fL)

void getfL_NP()
{
  // get the mass background fit function
  TF1 *fMass = new TF1("fMass", "1.-x*[0]", 3.35, 4.);
  // define same function as above but *m
  TF1 *mMass = new TF1("mMass", "(1.-x*[0])*x", 3.35, 4.);
  // get fMass parameters
  TFile *inFMass = new TFile("files/mfit_NP.root");
  TGraphErrors *m_ld = (TGraphErrors*)inFMass->Get("fit_m_bkg");
  inFMass->Close();

  // get the binning
  int nBinsY = m_ld->GetN();
  const double *yBins = m_ld->GetX();

  double m_min[] = {3.4, 3.57, 3.82};
  double m_max[] = {3.52, 3.81, 4.0};

  // this part is done for every pT bin
  double avg_LSB[nBinsY], avg_RSB[nBinsY], fL[nBinsY];
  double pt_v[nBinsY], pt_e[nBinsY], zeros[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    pt_v[i] = yBins[i];
    pt_e[i] = m_ld->GetEX()[i];
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
  
  TH1F *fmL = gPad->DrawFrame(15, m_min[0], 105, m_max[0]);
  fmL->SetXTitle("p_{T} (GeV)");
  fmL->SetYTitle("<m_{LSB}> (GeV)");
  fmL->GetYaxis()->SetTitleOffset(1.3);
  fmL->GetYaxis()->SetLabelOffset(0.01);
  fmL->SetTitle("Run 2 <m_{LSB}> (p_{T})");

  g_mL->SetLineColor(kBlack);
  g_mL->SetMarkerStyle(20);
  g_mL->SetMarkerSize(.75);
  g_mL->SetMarkerColor(kBlack);
  g_mL->Draw("psame");

  c2->cd(2);
  
  TH1F *fmR = gPad->DrawFrame(15, m_min[2], 105, m_max[2]);
  fmR->SetXTitle("p_{T} (GeV)");
  fmR->SetYTitle("<m_{RSB}> (GeV)");
  fmR->GetYaxis()->SetTitleOffset(1.3);
  fmR->GetYaxis()->SetLabelOffset(0.01);
  fmR->SetTitle("Run 2 <m_{RSB}> (p_{T})");

  g_mR->SetLineColor(kBlue);
  g_mR->SetMarkerStyle(20);
  g_mR->SetMarkerSize(.75);
  g_mR->SetMarkerColor(kBlue);
  g_mR->Draw("psame");

  c2->cd(3);
  
  TH1F *ffL = gPad->DrawFrame(15, 50, 105, 55);
  ffL->SetXTitle("p_{T} (GeV)");
  ffL->SetYTitle("f_{L} (%)");
  ffL->GetYaxis()->SetTitleOffset(1.3);
  ffL->GetYaxis()->SetLabelOffset(0.01);
  ffL->SetTitle("Run 2 f_{L} (p_{T})");

  g_fL->SetLineColor(kRed);
  g_fL->SetMarkerStyle(20);
  g_fL->SetMarkerSize(.75);
  g_fL->SetMarkerColor(kRed);
  g_fL->Draw("psame");
  
  c2->SaveAs(Form("plots/fL_NP.pdf"));
  c2->Clear();
  c2->Destructor();

  TFile *fout = new TFile("files/store_fL_NP.root", "recreate");

  for(int i = 0; i < nBinsY; i++) g_fL->GetY()[i]/=100.;
  g_fL->Write("g_fL");
  fout->Close();
}
