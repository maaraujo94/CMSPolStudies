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
  TH2D *h_MC2d = new TH2D(); // MC 2d map
  TFile *inHist = new TFile("files/histoStore.root");
  inHist->GetObject("dataH_ab", h_PR2d);
  h_PR2d->SetDirectory(0);
  inHist->GetObject("mcH_ab", h_MC2d);
  h_MC2d->SetDirectory(0);
  inHist->Close();
  
  // get the binning
  int nBinsX = h_PR2d->GetNbinsX(), nBinsY = h_PR2d->GetNbinsY();
  const double *yBins = h_PR2d->GetYaxis()->GetXbins()->GetArray();
  double minX = h_PR2d->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_PR2d->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/(double)nBinsX;

  // get the bkg model functions
  TH1D **h_NPr = new TH1D*[nBinsY]; // NP background 1d histos
  TH1D **h_SBr = new TH1D*[nBinsY]; // SB background 1d histos
  TFile *inBkg = new TFile("files/bkgCosModel.root");
  for(int i = 0; i < nBinsY; i++) {
    inBkg->GetObject(Form("h_NP_%d", i), h_NPr[i]);
    h_NPr[i]->SetDirectory(0);
    inBkg->GetObject(Form("h_SB_%d", i), h_SBr[i]);
    h_SBr[i]->SetDirectory(0);
  }
  inBkg->Close();
  
  // get the fit range from our cosmax(pT)
  ifstream in;
  string dataS;
  in.open("text_output/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", yBins[0]-10, yBins[nBinsY]+10);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);
  
  // bkg fraction in PR SR
  TF1 *f_fNP = new TF1();
  TF1 *f_fSB = new TF1();
  TFile *inFrac = new TFile("files/bkgFrac.root");
  inFrac->GetObject("fit_NP", f_fNP);
  inFrac->GetObject("fit_SB", f_fSB);
  inFrac->Close();
    
  // this part is done for every pT bin
  TFile *fout = new TFile("files/bkgSubRes.root", "recreate");
  TCanvas *c =  new TCanvas("", "", 900, 900);
  double lambda_1[nBinsY], lambda_2[nBinsY], lambda_3[nBinsY], lambda_4[nBinsY];
  double elambda_1[nBinsY], elambda_2[nBinsY], elambda_3[nBinsY], elambda_4[nBinsY];
  double dyBins[nBinsY], vyBins[nBinsY];

  // histos with all the 2d dists to be stored
  TH2D *h_Datas = new TH2D("h_Data", "Data/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_NPs = new TH2D("h_NP", "NP/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_SBs = new TH2D("h_SB", "SB/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_PRs = new TH2D("h_PR", "PR/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_Js = new TH2D("h_J", "Prompt J/#psi/MC", nBinsX, minX, maxX, nBinsY, yBins);

  for(int i = 0; i < nBinsY; i++) {
    double pt_min = yBins[i], pt_max = yBins[i+1];
    dyBins[i] = 0.5* (yBins[i+1] - yBins[i]);
    vyBins[i] = 0.5* (yBins[i+1]+yBins[i]);

    // get number of data events in PR SR
    double N_sig = h_PR2d->Integral(1, nBinsX, i+1, i+1);
    
    // getting the max costh value for the fit, cR
    double cMaxVal = cosMax->Integral(pt_min, pt_max)/(pt_max-pt_min);
    double cR = floor(cMaxVal*10.)/10.;
    if(cMaxVal-cR>0.05) cR += 0.05;
    
    // get the base data and MC 1d projections
    TH1D *h_PR = h_PR2d->ProjectionX(Form("h_PRSR_%d", i), i+1, i+1);
    TH1D *h_MC = h_MC2d->ProjectionX(Form("h_MC_%d", i), i+1, i+1);

    // multiply analytical functions by MC - NP
    TH1D *h_NP = new TH1D();
    h_NP = (TH1D*)h_NPr[i]->Clone(Form("h_NPSR_%d", i));
    h_NP->Sumw2();
    h_NP->Multiply(h_MC);
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
    double f_NP = f_fNP->Eval(vyBins[i]);
    double scFac = f_NP * N_sig;
    // scale the background dist
    h_NP->Scale(scFac);

    // get the proper scaling factor out of the f_bkg
    double f_SB = f_fSB->Eval(vyBins[i]);
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

    // PART 5 - plotting

    h_PR->Sumw2();
    h_PR->Divide(h_MC);
    h_PR->SetTitle(Form("Signal extraction (%.0f < p_{T} < %.0f GeV)", pt_min, pt_max));
    h_PR->SetStats(0);
    h_PR->SetLineColor(kViolet);
    h_PR->SetMarkerColor(kViolet);
    h_PR->SetMinimum(0);
    h_PR->SetMaximum(h_PR->GetBinContent(1)*1.5);
    h_PR->Draw("error");

    h_NP->Divide(h_MC);
    h_NP->SetLineColor(kRed);
    h_NP->SetMarkerColor(kRed);
    h_NP->Draw("same");

    h_SB->Divide(h_MC);
    h_SB->SetLineColor(kGreen);
    h_SB->SetMarkerColor(kGreen);
    h_SB->Draw("same");

    h_purePR->Divide(h_MC);
    h_purePR->SetLineColor(kBlue);
    h_purePR->SetMarkerColor(kBlue);
    h_purePR->Draw("error same");

    TF1 *fcb1 = new TF1("fcb1", "[0]*(1+[1]*x*x)", minX, cR);
    fcb1->SetLineColor(kViolet);
    fcb1->SetLineStyle(kDashed);
    h_PR->Fit("fcb1", "R");
    TF1 *fcb2 = new TF1("fcb2", "[0]*(1+[1]*x*x)", minX, cR);
    h_NP->Fit("fcb2", "R0");
    TF1 *fcb3 = new TF1("fcb3", "[0]*(1+[1]*x*x)", minX, cR);
    h_justPR->Divide(h_MC);
    h_justPR->Fit("fcb3", "R0");
    TF1 *fcb4 = new TF1("fcb4", "[0]*(1+[1]*x*x)", minX, cR);
    fcb4->SetLineColor(kBlue);
    fcb4->SetLineStyle(kDashed);
    h_purePR->Fit("fcb4", "R");

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
  
    lambda_1[i] = fcb1->GetParameter(1);
    elambda_1[i] = fcb1->GetParError(1);
    lambda_2[i] = fcb2->GetParameter(1);
    elambda_2[i] = fcb2->GetParError(1);
    lambda_3[i] = fcb3->GetParameter(1);
    elambda_3[i] = fcb3->GetParError(1);
    lambda_4[i] = fcb4->GetParameter(1);
    elambda_4[i] = fcb4->GetParError(1);

    TLatex lc;
    lc.SetTextSize(0.03);
    lc.DrawLatex(0.1, h_PR->GetMaximum()*0.9, Form("#lambda_{#theta}^{PR SR} = %.3f #pm %.3f", fcb1->GetParameter(1), fcb1->GetParError(1)));
    lc.DrawLatex(0.1, h_PR->GetMaximum()*0.8, Form("#lambda_{#theta}^{prompt J/#psi} = %.3f #pm %.3f", fcb4->GetParameter(1), fcb4->GetParError(1)));
    
    TLine *c_lim = new TLine(cR, 0, cR, h_PR->GetMaximum());
    c_lim->SetLineStyle(kDashed);
    c_lim->SetLineColor(kBlack);
    c_lim->Draw();

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(h_PR, "Total", "pl");
    leg->AddEntry(h_NP, "NP contrib", "pl");
    leg->AddEntry(h_SB, "SB contrib", "pl");
    leg->AddEntry(h_purePR, "prompt J/#psi", "pl");
    leg->Draw();
    
    c->SaveAs(Form("plots/bkgSub/bin_%d.pdf", i));
    c->Clear();

  }
  c->Destructor();
  
  // PART 6 - output
  h_Datas->Write();
  h_NPs->Write();
  h_SBs->Write();
  h_PRs->Write();
  h_Js->Write();
  
  TGraphErrors *g_ld1 = new TGraphErrors(nBinsY, vyBins, lambda_1, dyBins, elambda_1);
  TGraphErrors *g_ld2 = new TGraphErrors(nBinsY, vyBins, lambda_2, dyBins, elambda_2);
  TGraphErrors *g_ld3 = new TGraphErrors(nBinsY, vyBins, lambda_3, dyBins, elambda_3);
  TGraphErrors *g_ld4 = new TGraphErrors(nBinsY, vyBins, lambda_4, dyBins, elambda_4);
  g_ld1->SetName("lambda_total");
  g_ld1->Write();
  g_ld2->SetName("lambda_NP");
  g_ld2->Write();
  g_ld3->SetName("lambda_prompt");
  g_ld3->Write();
  g_ld4->SetName("lambda_J");
  g_ld4->Write();
  fout->Close();

}
