// macro to fit the mass background

// aux function for fitting around the peak region
double bkg(double m, double A, double B, double mlc, double muc)
{
  if(m > mlc && m < muc)
    {
      TF1::RejectPoint();
    }
  double d_m = m-3.686;
    
  return A * ( 1 - 0.05 * d_m / B );
}

// aux function to assign colors to different plotting sets
double col(double v)
{
  if(v==0) return kRed;
  else if(v==1) return kBlue;
  else if(v==2) return kBlue;
  else if(v==3) return kRed;
  else if(v==4) return kBlack;
  else return 0;
}

// error propagation in division
double sig_div(double a, double b, double s_a, double s_b)
{
  return a/b*sqrt(s_a*s_a/(a*a) + s_b*s_b/(b*b));
}

// MAIN
void binnedSub()
{
  double Mq = 3.686;
  
  // read the data |costh| histograms
  string wname[3] = {"L", "S", "R"};
  
  TH2D **dataHist_ab = new TH2D*[3];

  TFile *fin = new TFile("files/ratioHist.root");
  for(int i_w = 0; i_w < 3; i_w++) {
    dataHist_ab[i_w] = (TH2D*)fin->Get(Form("dataH_ab_%s", wname[i_w].c_str()));
    dataHist_ab[i_w]->SetDirectory(0);
  }
  fin->Close();

  // doing the binning from the data
  int nBinsY = dataHist_ab[0]->GetNbinsY();
  const double* binsY = dataHist_ab[0]->GetYaxis()->GetXbins()->GetArray();
  int nBinsX = dataHist_ab[0]->GetNbinsX();
  double minX = dataHist_ab[0]->GetXaxis()->GetBinLowEdge(1);
  double maxX = dataHist_ab[0]->GetXaxis()->GetBinUpEdge(nBinsX);

  // histo of mass
  double m_min[] = {3.4, 3.565, 3.805};
  double m_max[] = {3.52, 3.805, 4.0};
  
  TH1D **h_PM = new TH1D*[nBinsY];
  double mbin[nBinsY], minH[nBinsY], maxH[nBinsY];
  int nbins[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    maxH[i] = 4.;
    mbin[i] = 0.015;
    minH[i] = 3.355;
    nbins[i] = ceil((maxH[i]-minH[i])/mbin[i]);
  }
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    double pt_min = binsY[i_pt], pt_max = binsY[i_pt+1];
    h_PM[i_pt] = new TH1D(Form("h_PM_%d", i_pt), Form("2017 #psi(2S) Mass (%.0f < p_{T} < %.0f GeV)", pt_min, pt_max), nbins[i_pt], minH[i_pt], maxH[i_pt]);
  }
  TGraph **g_pull = new TGraph*[nBinsY];
  TGraph **g_dev = new TGraph*[nBinsY];
  
  // get data tree, branches, etc
  TFile *infile = new TFile("../Store_data_codes/data17_cos.root");
  TTree *tree = (TTree*)infile->Get("data_cos");
  
  int nEvt = tree->GetEntries();
  
  Double_t data_mass, lts, th, phi, data_pt;
  
  tree->SetBranchAddress("Mass", &data_mass);
  tree->SetBranchAddress("lts", &lts);
  tree->SetBranchAddress("theta", &th);
  tree->SetBranchAddress("phi", &phi);
  tree->SetBranchAddress("dimPt", &data_pt);

  // go through tree
  for(int i = 0; i < nEvt; i++)
    {
      tree->GetEntry(i);
      for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
	double pt_min = binsY[i_pt], pt_max = binsY[i_pt+1];
	if(abs(lts) < 2.5 && data_pt < pt_max && data_pt > pt_min) {	
	  h_PM[i_pt]->Fill(data_mass);
	}
      }
    }
  cout << "Mass histograms filled for all pT bins" << endl;

  TCanvas *c = new TCanvas("", "", 700, 700);

  ofstream fout;
  fout.open("text_output/mbkg_fit.tex");
  fout << "\\begin{tabular}{c|c||c|c|c||c|c|c}\n";
  fout << "$\\pt$ (GeV) & events & $A$ (MeV$^{-1}$) & $B$ (MeV)  & $\\chi^2$/ndf & all evt (Peak) & bkg evt (Peak) & bkg/all (Peak) (\\%)\\\\\n";
  fout << "\\hline\n";
  fout.close();
  
  TFile *tfout = new TFile("files/mbSub.root", "RECREATE");
  tfout->Close();

  // store the bkg/all (sig) fraction for plotting later
  double x_pt[nBinsY], y_fb[nBinsY], ex_pt[nBinsY], ey_fb[nBinsY];
  
  // cycle to go over all pT bins
  double min = m_min[0], max = m_max[2];
  TF1 **massFit = new TF1*[nBinsY];

  for(int i_pt = 0; i_pt < nBinsY; i_pt++)
    {
      double pt_min = binsY[i_pt], pt_max = binsY[i_pt+1];
      
      // fit B mass
      double ymax = h_PM[i_pt]->GetMaximum();
      double A = ymax/200., B = 0.1;
      double mlc = m_max[0], muc = m_min[2];
     	
      cout << "A init: " << A << " ||| B init: " << B << endl << endl; 
      massFit[i_pt] = new TF1(Form("mFit_%d", i_pt), "bkg(x, [0], [1], [2], [3])*[4]", min, max);
      massFit[i_pt]->SetParameters(A, B, mlc, muc, mbin[i_pt]*1e3);
      massFit[i_pt]->SetParNames("A", "B", "mlc", "muc", "bin");
      massFit[i_pt]->FixParameter(2, mlc);
      massFit[i_pt]->FixParameter(3, muc);
      massFit[i_pt]->FixParameter(4, mbin[i_pt]*1e3);

      // plot fit
      c->SetLeftMargin(0.11);
      h_PM[i_pt]->SetMaximum(ymax*1.1);
      h_PM[i_pt]->SetMinimum(0);
      h_PM[i_pt]->SetStats(0);
      h_PM[i_pt]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
      h_PM[i_pt]->GetYaxis()->SetTitle(Form("Events per %.0f MeV", mbin[i_pt]*1000));
      h_PM[i_pt]->GetYaxis()->SetTitleOffset(1.6);
      h_PM[i_pt]->Draw("error");
      TFitResultPtr mF_r = h_PM[i_pt]->Fit(Form("mFit_%d", i_pt), "RS");
      
      TLine *cut1 = new TLine(min, 0, min, ymax*1.1);
      cut1->SetLineColor(kRed);
      cut1->SetLineStyle(kDashed);
      cut1->Draw();
      TLine *cut2 = new TLine(mlc, 0, mlc, ymax*1.1);
      cut2->SetLineColor(kRed);
      cut2->SetLineStyle(kDashed);
      cut2->Draw();
      TLine *cut3 = new TLine(max, 0, max, ymax*1.1);
      cut3->SetLineColor(kRed);
      cut3->SetLineStyle(kDashed);
      cut3->Draw();
      TLine *cut4 = new TLine(muc, 0, muc, ymax*1.1);
      cut4->SetLineColor(kRed);
      cut4->SetLineStyle(kDashed);
      cut4->Draw();

      // calculate number of data events in each region
      double nLSB = tree->GetEntries(Form("abs(lts)<2.5 && dimPt < %f && dimPt > %f && Mass < %f && Mass > %f", pt_max, pt_min, m_max[0], m_min[0]));
      double nSig = tree->GetEntries(Form("abs(lts)<2.5 && dimPt < %f && dimPt > %f && Mass < %f && Mass > %f", pt_max, pt_min, m_max[1], m_min[1]));
      double nRSB = tree->GetEntries(Form("abs(lts)<2.5 && dimPt < %f && dimPt > %f && Mass < %f && Mass > %f", pt_max, pt_min, m_max[2], m_min[2]));

      massFit[i_pt]->SetParameter(2, 4);
      massFit[i_pt]->SetParameter(3, 3);
      double bkg_in_Sig = massFit[i_pt]->Integral(m_min[1], m_max[1])/mbin[i_pt];
 
      fout.open("text_output/mbkg_fit.tex", std::ofstream::app);
      fout << Form("$[%.0f, %.0f]$", pt_min, pt_max) << " & ";
      fout << h_PM[i_pt]->GetEntries() << " & "; 
      fout << Form("%.3f", massFit[i_pt]->GetParameter(0) / (pt_max-pt_min)) << "$\\pm$" << Form("%.3f", massFit[i_pt]->GetParError(0) / (pt_max-pt_min)) << " & ";
      fout << Form("%.1f", massFit[i_pt]->GetParameter(1)*1e3) << "$\\pm$" << Form("%.1f", massFit[i_pt]->GetParError(1)*1e3) << " & ";
      fout << Form("%.0f", massFit[i_pt]->GetChisquare()) << "/" << massFit[i_pt]->GetNDF();
      fout << Form(" & %.0f & %.0f & %.1f ", nSig, bkg_in_Sig, bkg_in_Sig/nSig*100.) << " \\\\\n";
      fout.close();

      x_pt[i_pt] = 0.5*(pt_max+pt_min);
      ex_pt[i_pt] = 0.5*(pt_max-pt_min);
      y_fb[i_pt] = bkg_in_Sig/nSig; 
      double e_BiS = massFit[i_pt]->IntegralError(m_min[1], m_max[1] , mF_r->GetParams(), mF_r->GetCovarianceMatrix().GetMatrixArray() )/mbin[i_pt];
      double e_nS = sqrt(nSig);
      ey_fb[i_pt] = sig_div(bkg_in_Sig, nSig, e_BiS, e_nS);
      
      c->SaveAs(Form("plots/bkgSub/B_fit_%d.pdf", i_pt));
      c->Clear();

      // getting the pulls and relative deviations
      double m_val[nbins[i_pt]], p_val[nbins[i_pt]], d_val[nbins[i_pt]];
      for(int i_b = 0; i_b < nbins[i_pt]; i_b++) {
	double m_l = minH[i_pt]+i_b*mbin[i_pt], m_u = minH[i_pt]+(i_b+1)*mbin[i_pt];
	m_val[i_b] = 0.5*(m_l+m_u);
	// check if it's LSB or RSB
	if((m_l >= min && m_u <= mlc) || (m_l >= muc && m_u <= max)) {
	  double data_val = h_PM[i_pt]->GetBinContent(i_b+1);
	  double data_err = h_PM[i_pt]->GetBinError(i_b+1);
	  double fit_val = massFit[i_pt]->Eval(m_val[i_b]);
	  p_val[i_b] = (data_val-fit_val)/data_err;
	  d_val[i_b] = (data_val-fit_val)/fit_val;
	}
	else {
	  p_val[i_b] = 0;
	  d_val[i_b] = 0;
	}
      } 

      // plotting pulls
      TH1F *fc_p = c->DrawFrame(minH[i_pt], -5, maxH[i_pt], 5);
      fc_p->SetXTitle("M(#mu#mu) (GeV)");
      fc_p->SetYTitle("pulls");
      fc_p->GetYaxis()->SetTitleOffset(1.3);
      fc_p->SetTitle(Form("2017 Bkg fit pulls (%.0f < p_{T} < %.0f GeV)", pt_min, pt_max));

      g_pull[i_pt] = new TGraph(nbins[i_pt], m_val, p_val);
      g_pull[i_pt]->SetName(Form("g_pull_%d", i_pt));
      g_pull[i_pt]->SetMarkerStyle(20);
      g_pull[i_pt]->SetMarkerColor(kBlack);
      g_pull[i_pt]->SetMarkerSize(.75);
      g_pull[i_pt]->SetLineColor(kBlack);
      g_pull[i_pt]->Draw("plsame");

      TLine *cutp1 = new TLine(min, -5, min, 5);
      cutp1->SetLineColor(kRed);
      cutp1->SetLineStyle(kDashed);
      cutp1->Draw();
      TLine *cutp2 = new TLine(mlc, -5, mlc, 5);
      cutp2->SetLineColor(kRed);
      cutp2->SetLineStyle(kDashed);
      cutp2->Draw();
      TLine *cutp3 = new TLine(max, -5, max, 5);
      cutp3->SetLineColor(kRed);
      cutp3->SetLineStyle(kDashed);
      cutp3->Draw();
      TLine *cutp4 = new TLine(muc, -5, muc, 5);
      cutp4->SetLineColor(kRed);
      cutp4->SetLineStyle(kDashed);
      cutp4->Draw();
      
      c->SaveAs(Form("plots/bkgSub/B_pull_%d.pdf", i_pt));
      c->Clear();

      TH1F *fc_d = c->DrawFrame(minH[i_pt], -0.3, maxH[i_pt], 0.3);
      fc_d->SetXTitle("M(#mu#mu) (GeV)");
      fc_d->SetYTitle("Fit deviation");
      fc_d->GetYaxis()->SetTitleOffset(1.3);
      fc_d->SetTitle(Form("2017 Bkg fit deviation (%.0f < p_{T} < %.0f GeV)", pt_min, pt_max));

      g_dev[i_pt] = new TGraph(nbins[i_pt], m_val, d_val);
      g_dev[i_pt]->SetName(Form("g_dev_%d", i_pt));
      g_dev[i_pt]->SetMarkerStyle(20);
      g_dev[i_pt]->SetMarkerColor(kBlack);
      g_dev[i_pt]->SetMarkerSize(.75);
      g_dev[i_pt]->SetLineColor(kBlack);
      g_dev[i_pt]->Draw("plsame");

      TLine *cutd1 = new TLine(min, -0.3, min, 0.3);
      cutd1->SetLineColor(kRed);
      cutd1->SetLineStyle(kDashed);
      cutd1->Draw();
      TLine *cutd2 = new TLine(mlc, -0.3, mlc, 0.3);
      cutd2->SetLineColor(kRed);
      cutd2->SetLineStyle(kDashed);
      cutd2->Draw();
      TLine *cutd3 = new TLine(max, -0.3, max, 0.3);
      cutd3->SetLineColor(kRed);
      cutd3->SetLineStyle(kDashed);
      cutd3->Draw();
      TLine *cutd4 = new TLine(muc, -0.3, muc, 0.3);
      cutd4->SetLineColor(kRed);
      cutd4->SetLineStyle(kDashed);
      cutd4->Draw();
      
      c->SaveAs(Form("plots/bkgSub/B_dev_%d.pdf", i_pt));
      c->Clear();
      
      TFile *tfout = new TFile("files/mbSub.root", "UPDATE");
      h_PM[i_pt]->Write();
      g_pull[i_pt]->Write();
      tfout->Close();
    }

  fout.open("text_output/mbkg_fit.tex", std::ofstream::app);
  fout << "\\end{tabular}\n";
  fout.close();

  TH1F *fc = c->DrawFrame(20, 0, binsY[nBinsY], 1.);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("bkg fraction");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->SetTitle("2017 Fraction of bkg events in peak region");

  TGraphErrors* g_frac = new TGraphErrors(nBinsY, x_pt, y_fb, ex_pt, ey_fb);
  g_frac->SetLineColor(kBlack);
  g_frac->SetMarkerColor(kBlack);
  g_frac->SetMarkerStyle(20);
  g_frac->SetMarkerSize(.75);
  g_frac->Draw("psame");

  c->SaveAs("plots/bkgSub/bkg_frac.pdf");
  c->Clear();

  double Av[nBinsY], eAv[nBinsY], Bv[nBinsY], eBv[nBinsY], aNv[nBinsY], eaNv[nBinsY];
  double chiP[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    // normalize A by width of mass binning
    Av[i] = massFit[i]->GetParameter(0);
    Bv[i] = massFit[i]->GetParameter(1)*1e3;
    eAv[i] = massFit[i]->GetParError(0);
    eBv[i] = massFit[i]->GetParError(1)*1e3;
    aNv[i] = Av[i] / (binsY[i+1]-binsY[i]);
    eaNv[i] = eAv[i] / (binsY[i+1]-binsY[i]);

    chiP[i] = TMath::Prob(massFit[i]->GetChisquare(), massFit[i]->GetNDF());
  }

  c->SetLogy();
  TH1F *fcAN = c->DrawFrame(20, 1e-1, binsY[nBinsY], 2e2);
  fcAN->SetXTitle("p_{T} (GeV)");
  fcAN->SetYTitle("A (MeV^{-1}) / pT bin width");
  fcAN->GetYaxis()->SetTitleOffset(1.3);
  fcAN->SetTitle("2017 A/pT bin width (p_{T})");

  TGraphErrors* g_AN = new TGraphErrors(nBinsY, x_pt, aNv, ex_pt, eaNv);
  g_AN->SetLineColor(kBlack);
  g_AN->SetMarkerColor(kBlack);
  g_AN->SetMarkerStyle(20);
  g_AN->SetMarkerSize(.75);
  g_AN->Draw("psame");

  
  c->SaveAs("plots/bkgSub/fit_AN.pdf");
  c->Clear();

  c->SetLogy(0);
  TH1F *fcB = c->DrawFrame(20, 50, binsY[nBinsY], 150);
  fcB->SetXTitle("p_{T} (GeV)");
  fcB->SetYTitle("B (MeV)");
  fcB->GetYaxis()->SetTitleOffset(1.3);
  fcB->SetTitle("2017 B (p_{T})");

  TGraphErrors* g_B = new TGraphErrors(nBinsY, x_pt, Bv, ex_pt, eBv);
  g_B->SetLineColor(kBlack);
  g_B->SetMarkerColor(kBlack);
  g_B->SetMarkerStyle(20);
  g_B->SetMarkerSize(.75);
  g_B->Draw("psame");

  c->SaveAs("plots/bkgSub/fit_B.pdf");
  c->Clear();

  c->SetLogy(0);
  TH1F *fcchi = c->DrawFrame(20, 0, binsY[nBinsY], 1);
  fcchi->SetXTitle("p_{T} (GeV)");
  fcchi->SetYTitle("P(#chi^{2}, ndf)");
  fcchi->GetYaxis()->SetTitleOffset(1.3);
  fcchi->SetTitle("2017 P(#chi^{2}, ndf)");

  TGraph* g_chi = new TGraph(nBinsY, x_pt, chiP);
  g_chi->SetLineColor(kBlack);
  g_chi->SetMarkerColor(kBlack);
  g_chi->SetMarkerStyle(20);
  g_chi->SetMarkerSize(1.5);
  g_chi->Draw("psame");

  c->SaveAs("plots/bkgSub/fit_m_chiP.pdf");
  c->Clear();

  TFile *tfout2 = new TFile("files/mbSub.root", "UPDATE");
  g_frac->SetName("bkg_frac");
  g_frac->Write();
  g_AN->SetName("fit_A_norm");
  g_AN->Write();
  g_B->SetName("fit_B");
  g_B->Write();
  tfout2->Close();
}
