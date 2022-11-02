// macro to get the fraction of LSB contributing to SR (fL)

int do_round(double val)
{
  int valR = (int)val;
  if (val-valR > 0.5) return valR+1;
  else return valR;
}

void getfL()
{
  // get the mass background fit function
  TF1 *fMass = new TF1("fMass", "exp(-x/[0])", 2.9, 3.3);
  // define same function as above but *m
  TF1 *mMass = new TF1("mMass", "exp(-x/[0])*x", 2.9, 3.3);
  // get fMass parameters
  TFile *inFMass = new TFile("../bkgFits/files/mfit_NP.root");
  TGraphErrors *m_ld = (TGraphErrors*)inFMass->Get("fit_lambda");
  inFMass->Close();

  // get the binning
  int nBinsY = m_ld->GetN();
  const double *yBins = m_ld->GetX();
  
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
    avg_LSB[i] = mMass->Integral(2.92, 2.95)/fMass->Integral(2.92, 2.95);
    avg_RSB[i] = mMass->Integral(3.21, 3.28)/fMass->Integral(3.21, 3.28);
    fL[i] = (3.1 - avg_LSB[i]) / (avg_RSB[i] - avg_LSB[i]);
    
    fL[i] *= 100;

    cout << i << " " << fL[i] << endl;
  }

  TGraphErrors *g_mL = new TGraphErrors(nBinsY, pt_v, avg_LSB, pt_e, zeros); 
  TGraphErrors *g_mR = new TGraphErrors(nBinsY, pt_v, avg_RSB, pt_e, zeros); 
  TGraphErrors *g_fL = new TGraphErrors(nBinsY, pt_v, fL, pt_e, zeros);
  
  TCanvas *c2 =  new TCanvas("", "", 600, 800);
  
  c2->Divide(1, 3);
  c2->cd(1);
  
  TH1F *fmL = gPad->DrawFrame(20, 2.92, 125, 2.95);
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
  
  TH1F *fmR = gPad->DrawFrame(20, 3.21, 125, 3.28);
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

  for(int i = 0; i < nBinsY; i++) g_fL->GetY()[i]/=100.;
  g_fL->Write("g_fL");
  fout->Close();
}
