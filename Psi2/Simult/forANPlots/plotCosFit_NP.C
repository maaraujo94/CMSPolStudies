#import "../cosMax/imp_jumpF.C"

// macro to plot fitted costh dists for slides
void plotCosFit_NP()
{
  // read the histos from subtraction - normalized by f_NP/f_bkg
  TFile *infile = new TFile("../NP_fit/files/bkgSubRes.root");
  const int n_inp = 3;
  TH2D **h_rat = new TH2D*[n_inp];
  string lbl[] = {"NP", "SB", "NPc"};
  for(int i = 0; i < n_inp; i++) {
    infile->GetObject(Form("h_%s", lbl[i].c_str()), h_rat[i]);
    h_rat[i]->SetDirectory(0);
  }
  infile->Close();

  // get the binning
  int nBinsX = h_rat[0]->GetNbinsX(), nBinsY = h_rat[0]->GetNbinsY();
  const double *yBins = h_rat[0]->GetYaxis()->GetXbins()->GetArray();

  // get the 1d plots
  TH1D *h_rat1d[n_inp][nBinsY];
  for(int i_t = 0; i_t < n_inp; i_t++) {
    for(int i = 1; i <= nBinsY; i++) {
      h_rat1d[i_t][i-1] = h_rat[i_t]->ProjectionX(Form("bin%d_%d_r", i, i_t+1), i, i);
    }
  }

  // the fit function to be used
  TF1 **fit1d = new TF1*[n_inp];
  for(int i = 0; i < n_inp; i++) {
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
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  int cols[] = {kRed+3, kGreen+2, kRed};
    
  for(int i = 0; i < nBinsY; i++) {
    // get pt vars
    double pMin = h_rat[0]->GetYaxis()->GetBinLowEdge(i+1);
    double pMax = h_rat[0]->GetYaxis()->GetBinUpEdge(i+1);

    // get max costheta
    double cMaxVal = jumpF(cosMax->Integral(pMin, pMax)/(pMax-pMin))-0.05;

    // fit the 2 functions
    for(int i_t = 0; i_t < n_inp; i_t++) {
      fit1d[i_t]->SetRange(0, cMaxVal);
      fit1d[i_t]->SetParameters(h_rat1d[i_t][i]->GetBinContent(1)*1.1, 0.1);

      h_rat1d[i_t][i]->Fit(fit1d[i_t], "R0");
    }

    // draw all 3 fit results
    h_rat1d[0][i]->SetTitle("");
    h_rat1d[0][i]->SetStats(0);
    h_rat1d[0][i]->SetLineColor(cols[0]);
    h_rat1d[0][i]->SetMarkerColor(cols[0]);
    h_rat1d[0][i]->SetMinimum(0);
    h_rat1d[0][i]->SetMaximum(h_rat1d[0][i]->GetBinContent(1)*1.3);
    h_rat1d[0][i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    h_rat1d[0][i]->Draw("error");
    fit1d[0]->SetLineColor(cols[0]);
    fit1d[0]->SetLineStyle(kDashed);
    fit1d[0]->Draw("same");

    /* h_rat1d[1][i]->SetLineColor(cols[1]);
    h_rat1d[1][i]->SetMarkerColor(cols[1]);
    h_rat1d[1][i]->Draw("error same");*/
    
     h_rat1d[2][i]->SetLineColor(cols[2]);
    h_rat1d[2][i]->SetMarkerColor(cols[2]);
    h_rat1d[2][i]->Draw("error same");
    fit1d[2]->SetLineColor(cols[2]);
    fit1d[2]->SetLineStyle(kDashed);
    fit1d[2]->Draw("same");

    TLatex lcr1;
    lcr1.SetTextSize(0.04);
    lcr1.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.9, "Run 2");
    lcr1.DrawLatex(0.7, h_rat1d[0][i]->GetMaximum()*0.85, Form("%.0f-%.0f GeV", pMin, pMax));
    lcr1.SetTextColor(cols[0]);
    lcr1.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.87, "NPS/MC");
    lcr1.DrawLatex(0.15, h_rat1d[0][i]->GetMaximum()*0.82, Form("#lambda_{#theta} = %.3f #pm %.3f", fit1d[0]->GetParameter(1), fit1d[0]->GetParError(1)));
    //lcr1.SetTextColor(cols[1]);
    //lcr1.DrawLatex(0.15, h_rat1d[1][i]->GetMaximum()*1.5, "bkg^{NP}/MC");
    lcr1.SetTextColor(cols[2]);
    lcr1.DrawLatex(0.15, h_rat1d[2][i]->GetMaximum()*0.85, "#psi(2S)_{NP}/MC");
    lcr1.DrawLatex(0.15, h_rat1d[2][i]->GetMaximum()*0.72, Form("#lambda_{#theta} = %.3f #pm %.3f", fit1d[2]->GetParameter(1), fit1d[2]->GetParError(1)));
    
    c->SaveAs(Form("plots/ratioFinal_NP/fits/bin1F_%d.pdf", i));
    c->Clear();
  }
  
  c->Destructor();
}
