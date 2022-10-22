#import "../cosMax/imp_jumpF.C"

// macro to plot fitted costh dists for slides
void plotCosFit()
{
  // read the histos from subtraction - normalized by f_NP/f_bkg
  TFile *infile = new TFile("../PR_fit/files/bkgSubRes.root");
  TH2D **h_rat = new TH2D*[5];
  string lbl[] = {"Data", "NP", "PR", "J", "SB"};
  for(int i = 0; i < 5; i++) {
    infile->GetObject(Form("h_%s", lbl[i].c_str()), h_rat[i]);
    h_rat[i]->SetDirectory(0);
  }
  infile->Close();

  // get the binning
  int nBinsX = h_rat[0]->GetNbinsX(), nBinsY = h_rat[0]->GetNbinsY();
  const double *yBins = h_rat[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *h_rat1d[5][nBinsY];
  for(int i_t = 0; i_t < 5; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      h_rat1d[i_t][i-1] = h_rat[i_t]->ProjectionX(Form("bin%d_%d_r", i, i_t+1), i, i);
    }
  }

  // the fit function to be used
  TF1 **fit1d = new TF1*[4];
  for(int i = 0; i < 4; i++) {
    fit1d[i] = new TF1(Form("fit_%d", i), "[0]*(1+[1]*x*x)", 0, 1);
    fit1d[i]->SetParNames("A", "l_th");
  }
  
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

 
  // the cycle to plot each bin
  TCanvas *c = new TCanvas("", "", 700, 700);    
  int cols[] = {kViolet-1, kRed, kBlack, kBlue, kGreen};
    
  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_rat[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_rat[0]->GetYaxis()->GetBinUpEdge(i+1);

    // get max costheta
    double cMaxVal = jumpF(cosMax->Integral(pMin, pMax)/(pMax-pMin));

    // fit the 4 functions
    for(int i_t = 0; i_t < 4; i_t++) {
      fit1d[i_t]->SetRange(0, cMaxVal);
      fit1d[i_t]->SetParameters(h_rat1d[i_t][i]->GetBinContent(1)*1.1, 0.1);

      h_rat1d[i_t][i]->Fit(fit1d[i_t], "R0");
    }

    // just peak/MC costh
    h_rat1d[0][i]->SetTitle("");
    h_rat1d[0][i]->SetStats(0);
    h_rat1d[0][i]->SetLineColor(cols[0]);
    h_rat1d[0][i]->SetMarkerColor(cols[0]);
    h_rat1d[0][i]->SetMinimum(0);
    h_rat1d[0][i]->SetMaximum(h_rat1d[0][i]->GetBinContent(1)*1.5);
    h_rat1d[0][i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    h_rat1d[0][i]->Draw("error");
    fit1d[0]->SetLineColor(cols[0]);
    fit1d[0]->SetLineStyle(kDashed);
    fit1d[0]->Draw("same");
    
    TLatex lcr1;
    lcr1.SetTextSize(0.04);
    lcr1.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.9, "2018");
    lcr1.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.85, Form("%.1f-%.1f GeV", pMin, pMax));
    lcr1.SetTextColor(cols[0]);
    lcr1.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.8, "Peak/MC");
    lcr1.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.75, Form("#lambda_{#theta} = %.3f #pm %.3f", fit1d[0]->GetParameter(1), fit1d[0]->GetParError(1)));
    
    c->SaveAs(Form("plots/ratioFinal/fits/bin1F_%d.pdf", i));

    // peak/MC, NP/MC, PR/MC costh
    // NP
    h_rat1d[1][i]->SetStats(0);
    h_rat1d[1][i]->SetLineColor(cols[1]);
    h_rat1d[1][i]->SetMarkerColor(cols[1]);
    h_rat1d[1][i]->Draw("error same");
 
    // PR
    h_rat1d[2][i]->SetStats(0);
    h_rat1d[2][i]->SetLineColor(cols[2]);
    h_rat1d[2][i]->SetMarkerColor(cols[2]);
    h_rat1d[2][i]->Draw("error same");
    fit1d[2]->SetLineColor(cols[2]);
    fit1d[2]->SetLineStyle(kDashed);
    fit1d[2]->Draw("same");

    lcr1.SetTextColor(cols[1]);
    lcr1.DrawLatex(0.15, h_rat1d[1][i]->GetMaximum()*1.1, "NP/MC");

    lcr1.SetTextColor(cols[2]);
    lcr1.DrawLatex(0.15, h_rat1d[2][i]->GetMaximum()*0.7, "PR/MC");
    lcr1.DrawLatex(0.15, h_rat1d[2][i]->GetMaximum()*0.6, Form("#lambda_{#theta} = %.3f #pm %.3f", fit1d[2]->GetParameter(1), fit1d[2]->GetParError(1)));

        
    c->SaveAs(Form("plots/ratioFinal/fits/bin2F_%d.pdf", i));
    c->Clear();


    // peak/MC, NP/MC, PR/MC, SB/MC,  costh
    // peak
    h_rat1d[0][i]->Draw("error");
    fit1d[0]->Draw("same");

    // NP
    h_rat1d[1][i]->Draw("error same");
 
    // PR
    h_rat1d[2][i]->Draw("error same");
    fit1d[2]->Draw("same");

    // Jpsi
    h_rat1d[3][i]->SetStats(0);
    h_rat1d[3][i]->SetLineColor(cols[3]);
    h_rat1d[3][i]->SetMarkerColor(cols[3]);
    h_rat1d[3][i]->Draw("error same");
    fit1d[3]->SetLineColor(cols[3]);
    fit1d[3]->SetLineStyle(kDashed);
    fit1d[3]->Draw("same");

    // bkg
    h_rat1d[4][i]->SetStats(0);
    h_rat1d[4][i]->SetLineColor(cols[4]);
    h_rat1d[4][i]->SetMarkerColor(cols[4]);
    h_rat1d[4][i]->Draw("error same");
 

    TLatex lcr3;
    lcr3.SetTextSize(0.04);
    lcr3.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.9, "2018");
    lcr3.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.85, Form("%.1f-%.1f GeV", pMin, pMax));

    lcr3.SetTextColor(cols[0]);
    lcr3.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.8, "Peak/MC");
    lcr3.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.75, Form("#lambda_{#theta} = %.3f #pm %.3f", fit1d[0]->GetParameter(1), fit1d[0]->GetParError(1)));
    
    lcr3.SetTextColor(cols[1]);
    lcr3.DrawLatex(0.15, h_rat1d[1][i]->GetMaximum()*1.1, "NP/MC");
    
    lcr3.SetTextColor(cols[2]);
    lcr3.DrawLatex(0.15, h_rat1d[2][i]->GetMaximum()*1.0, "PR/MC");
    lcr3.DrawLatex(0.15, h_rat1d[2][i]->GetMaximum()*0.93, Form("#lambda_{#theta} = %.3f #pm %.3f", fit1d[2]->GetParameter(1), fit1d[2]->GetParError(1)));

    lcr3.SetTextColor(cols[3]);
    lcr3.DrawLatex(0.15, h_rat1d[3][i]->GetMaximum()*0.75, "J/#psi/MC");
    lcr3.DrawLatex(0.15, h_rat1d[3][i]->GetMaximum()*0.65, Form("#lambda_{#theta} = %.3f #pm %.3f", fit1d[3]->GetParameter(1), fit1d[3]->GetParError(1)));

    lcr3.SetTextColor(cols[4]);
    lcr3.DrawLatex(0.15, h_rat1d[4][i]->GetMaximum()*1.1, "bkg/MC");

    c->SaveAs(Form("plots/ratioFinal/fits/bin3F_%d.pdf", i));
    c->Clear();

    cout << endl;
  }
  
  c->Destructor();
}
