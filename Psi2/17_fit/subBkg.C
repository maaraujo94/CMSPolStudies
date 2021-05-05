// macro to compile everything and do the background subtraction

// mass background model function
double bkg(double m, double A, double B)
{
  double d_m = m-3.686;
    
  return A * ( 1 - 0.05 * d_m / B );
}

// MAIN
void subBkg()
{
  // get the data/MC histogram
  TFile *fin = new TFile("files/ratioHist.root");
  TH2D *dataHist_ab = (TH2D*)fin->Get("ratioHist_ab_S");
  dataHist_ab->SetDirectory(0);
  fin->Close();

  // doing the binning from the data
  int nBinsY = dataHist_ab->GetNbinsY();
  const double* binsY = dataHist_ab->GetYaxis()->GetXbins()->GetArray();
  int nBinsX = dataHist_ab->GetNbinsX();
  double minX = dataHist_ab->GetXaxis()->GetBinLowEdge(1);
  double maxX = dataHist_ab->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/(double)nBinsX;

  // define the pure signal histogram
  TH2D *PSig_ab = new TH2D("PSig_ab", "Data (Pure Signal)", nBinsX, minX, maxX, nBinsY, binsY);

  // get SB costh model fit results (fix l4, free l2)
  TFile *fin_c = new TFile("files/mbSub.root");
  TGraphErrors* g_NL  = (TGraphErrors*)fin_c->Get("graph_NL_2");
  TGraphErrors* g_NR  = (TGraphErrors*)fin_c->Get("graph_NR_2");
  TGraphErrors* g_l2L = (TGraphErrors*)fin_c->Get("graph_l2L_2");
  TGraphErrors* g_l2R = (TGraphErrors*)fin_c->Get("graph_l2R_2");
  TGraphErrors* g_l4L = (TGraphErrors*)fin_c->Get("graph_l4L_2");
  TGraphErrors* g_l4R = (TGraphErrors*)fin_c->Get("graph_l4R_2");
  TGraphErrors* g_AN = (TGraphErrors*)fin_c->Get("fit_A_norm");
  TGraphErrors* g_B = (TGraphErrors*)fin_c->Get("fit_B");
  fin_c->Close();
  
  // get the fit range from our cosmax(pT)
  ifstream in;
  string dataS;
  in.open("text_output/cosMaxFitRes.txt");
  getline(in, dataS);
  getline(in, dataS);
  double maxPar[3], aux;
  in >> maxPar[0] >> aux >> maxPar[1] >> aux >> maxPar[2];
  in.close();
  
  TF1 *cosMax = new TF1("cosMax", "[0]*log([1]+[2]*x)", binsY[0]-10, binsY[nBinsY]+10);
  cosMax->SetParameters(maxPar[0], maxPar[1], maxPar[2]);

  // get data tree, branches, etc
  TFile *infile = new TFile("../Store_data_codes/data17_cos.root");
  TTree *tree = (TTree*)infile->Get("data_cos");

  // define the SB/MC costh model
  TF1 *f_cth = new TF1("f_cth", "[0]*(1+[1]*x*x+[2]*pow(x,4))", minX, maxX);
  f_cth->SetParNames("N", "l_2", "l_4");
  double N, l_2, l_4;

  // define the mass background model
  double m_min[] = {3.4, 3.565, 3.805};
  double m_max[] = {3.52, 3.805, 4.0};
  double mbin = 0.015;
  TF1 *f_mbkg = new TF1("f_mbkg", "bkg(x, [0], [1])*[2]", m_min[0], m_max[2]);
  f_mbkg->SetParNames("A", "B", "bin");
  double A_m, B_m;
  
  // cycle over all pT bins
  TCanvas *c = new TCanvas("", "", 700, 700);
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    double pt_min = binsY[i_pt], pt_max = binsY[i_pt+1];

    // getting the max costh value for the fit, cR
    double cMaxVal = cosMax->Integral(pt_min, pt_max)/(pt_max-pt_min);
    double cR = floor(cMaxVal*10.)/10.;
    if(cMaxVal-cR>0.05) cR += 0.05;

    // define the SB histos
    TH1D *h_LSB = new TH1D(Form("h_LSB_%d", i_pt), "2017 cos#theta", nBinsX, minX, maxX);
    TH1D *h_RSB = new TH1D(Form("h_RSB_%d", i_pt), "2017 cos#theta", nBinsX, minX, maxX);

    // fill histos with the model function
    N = g_NL->GetY()[i_pt];
    l_2 = g_l2L->GetY()[i_pt];
    l_4 = g_l4L->GetY()[i_pt];
    f_cth->SetParameters(N, l_2, l_4);
    for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
      double cos_v = minX+(i_cos+0.5)*dX;
      double f_v = f_cth->Eval(cos_v);
      h_LSB->SetBinContent(i_cos+1, f_v);
      h_LSB->SetBinError(i_cos+1, 0);
    }
    
    N = g_NR->GetY()[i_pt];
    l_2 = g_l2R->GetY()[i_pt];
    l_4 = g_l4R->GetY()[i_pt];
    f_cth->SetParameters(N, l_2, l_4);
    for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
      double cos_v = minX+(i_cos+0.5)*dX;
      double f_v = f_cth->Eval(cos_v);
      h_RSB->SetBinContent(i_cos+1, f_v);
      h_RSB->SetBinError(i_cos+1, 0);
    }
    
    double nLSB = tree->GetEntries(Form("abs(lts)<2.5 && dimPt < %f && dimPt > %f && Mass < %f && Mass > %f", pt_max, pt_min, m_max[0], m_min[0]));
    double nRSB = tree->GetEntries(Form("abs(lts)<2.5 && dimPt < %f && dimPt > %f && Mass < %f && Mass > %f", pt_max, pt_min, m_max[2], m_min[2]));

    // get averaged background
    TH1D *h_costh_bkg = new TH1D(Form("h_bkg_%d", i_pt), "cos#theta", nBinsX, minX, maxX);
    h_costh_bkg->Sumw2();
    h_costh_bkg->Add(h_LSB, h_RSB, nRSB/(nLSB+nRSB), nLSB/(nLSB+nRSB));

    // scale disto to the bkg in SR
    A_m = g_AN->GetY()[i_pt] * (pt_max-pt_min);
    B_m = g_B->GetY()[i_pt] / (1e3);
    f_mbkg->SetParameters(A_m, B_m, mbin*1e3);
    double bkg_in_Sig = f_mbkg->Integral(m_min[1], m_max[1])/mbin; // total background in signal region
    double nBkg = 2.*nLSB*nRSB/(nLSB+nRSB); // total integral of the bkg dist
    h_costh_bkg->Scale(bkg_in_Sig/nBkg); 

    // get the signal region distribution
    TH1D *h_SR = dataHist_ab->ProjectionX(Form("h_SR_%d", i_pt), i_pt+1, i_pt+1);

    // get the pure signal dist
    TH1D *h_PSig = new TH1D(Form("h_PSig_%d", i_pt), "cos#theta", nBinsX, minX, maxX);
    h_PSig->Sumw2();
    h_PSig->Add(h_SR, h_costh_bkg, 1, -1);

    for(int i_c = 0; i_c < nBinsX; i_c++) {
      if(h_SR->GetBinContent(i_c+1, i_pt+1) > 0) {
	PSig_ab->SetBinContent(i_c+1, i_pt+1, h_PSig->GetBinContent(i_c+1));
	PSig_ab->SetBinError(i_c+1, i_pt+1, h_PSig->GetBinError(i_c+1));
      }
      else {
	PSig_ab->SetBinContent(i_c+1, i_pt+1,0);
      }
    }
    
    h_SR->SetTitle(Form("2017 Signal extraction (%.0f < p_{T} < %.0f GeV)", pt_min, pt_max));
    h_SR->SetStats(0);
    h_SR->SetLineColor(kViolet);
    h_SR->SetMarkerColor(kViolet);
    h_SR->SetMinimum(0);
    h_SR->SetMaximum(h_SR->GetMaximum()*1.2);
    h_SR->Draw("error");
    
    h_costh_bkg->SetLineColor(kRed);
    h_costh_bkg->SetMarkerColor(kRed);
    h_costh_bkg->Draw("same");

    h_PSig->SetLineColor(kBlue);
    h_PSig->SetMarkerColor(kBlue);
    h_PSig->Draw("error same");

    TLine *c_lim = new TLine(cR, 0, cR, h_SR->GetMaximum());
    c_lim->SetLineStyle(kDashed);
    c_lim->SetLineColor(kBlack);
    c_lim->Draw();

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(h_SR, "Peak", "pl");
    leg->AddEntry(h_costh_bkg, "Background", "pl");
    leg->AddEntry(h_PSig, "Signal", "pl");
    leg->Draw();
     
    c->SaveAs(Form("plots/bkgSub/sub_bin_%d.pdf", i_pt));
    c->Clear();
  }
  c->Destructor();
  infile->Close();

  TFile *outfile = new TFile("files/ratioHist.root", "UPDATE");
  PSig_ab->Write(0, TObject::kOverwrite);
  outfile->Close();
}
