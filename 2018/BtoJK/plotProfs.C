double lambda_theta(double cos, double lambda_N)
{
  return -lambda_N*(1-3*cos)/(2+lambda_N*(1-cos));
}
double lambda_phi(double s1, double cos, double lambda_N)
{
  return lambda_N*s1/(2+lambda_N*(1-cos));
}
double lambda_thetaphi(double s2, double cos, double lambda_N)
{
  return lambda_N*s2/(2+lambda_N*(1-cos));
}

void plotProfs()
{
  TFile *infile = new TFile("../../Store_data_codes/2018/btopk_cos.root");
  TTree *tree = (TTree*)infile->Get("data_cos");
  
  int nEvt = tree->GetEntries();
  cout << nEvt << " data events after cuts" << endl;
  cout << tree->GetEntries("lts > 2.5") << " events with lts cut" << endl;
  cout << tree->GetEntries("lts > 2.5 && BMass > 5.24 && BMass < 5.32") << " events in signal region" << endl;

  // prepare binning and histograms for plots
  const int nPtBins = 45;
  double ptBins[nPtBins+1];
  for(int i=0; i<15; i++) ptBins[i] = (i+25.);
  for(int i=0; i<31; i++) ptBins[i+15] = (40.+2.*i);
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // making histograms
  // split into pT bins
  TH1D **h_costh_pt = new  TH1D*[nPtBins];
  TH1D **h_cosTH_pt = new  TH1D*[nPtBins];
  TH1D **h_phi_pt   = new  TH1D*[nPtBins];
  TH1D **h_PHI_pt   = new  TH1D*[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    h_costh_pt[i] = new TH1D(Form("h_costh_pt_%d", i), Form("cos#theta p_{T} bin %d: [%.0f, %.0f] GeV", i+1, ptBins[i], ptBins[i+1]), 100, -1.,   1.  );
    h_cosTH_pt[i] = new TH1D(Form("h_cosTH_pt_%d", i), Form("cos#theta p_{T} bin %d: [%.0f, %.0f] GeV", i+1, ptBins[i], ptBins[i+1]), 100, -1.,   1.  );
    h_phi_pt[i]   = new TH1D(Form("h_phi_pt_%d",   i), Form("#phi p_{T} bin %d: [%.0f, %.0f] GeV", i+1, ptBins[i], ptBins[i+1]),      100, -180., 180.);
    h_PHI_pt[i]   = new TH1D(Form("h_PHI_pt_%d",   i), Form("#phi p_{T} bin %d: [%.0f, %.0f] GeV", i+1, ptBins[i], ptBins[i+1]),      100, -180., 180.);
  }
  // phi diff
  TH1D *h_pdiff = new TH1D("h_pdif", "#phi - #Phi", 200, -180, 180);
  TProfile *p_pdiff = new TProfile("p_pdif", "#phi - #Phi vs p_{T}", nPtBins, ptBins, -180, 180);
  // profiles for lambdas
  TProfile *p_angC = new TProfile("p_angC", "cos^{2}#Theta vs p_{T}", nPtBins, ptBins, -1, 1);
  TProfile *p_angS1 = new TProfile("p_angS1", "sin^{2}#Thetacos(2#Phi) vs p_{T}", nPtBins, ptBins, -1, 1);
  TProfile *p_angS2 = new TProfile("p_angS2", "sin(2#Theta)cos#Phi vs p_{T}", nPtBins, ptBins, -1, 1);
 
  Double_t th, TH, phi, PHI, data_pt, data_mass, lts;

  tree->SetBranchAddress("th", &th);
  tree->SetBranchAddress("phi", &phi);
  tree->SetBranchAddress("TH", &TH);
  tree->SetBranchAddress("PHI", &PHI);
  tree->SetBranchAddress("JpsiPt", &data_pt);
  tree->SetBranchAddress("BMass", &data_mass);
  tree->SetBranchAddress("lts", &lts);

  double m_min[] = {5, 5.24, 5.36};
  double m_max[] = {5.2, 5.32, 5.6};

  double pdiff;
  // cycle over data, fill the costh histogram acc to B mass
  for(int i = 0; i < nEvt; i++)
    {
      tree->GetEntry(i);
      if(lts > 2.5) {
	if(data_mass > m_min[1] && data_mass < m_max[1]) {
	  pdiff = phi - PHI;
	  
	  if(pdiff > 180) pdiff -= 360;
	  if(pdiff < -180) pdiff += 360;

	  if(pdiff < -180 || pdiff > 180) cout << pdiff << endl;
	  
	  h_pdiff->Fill(pdiff);
	  p_pdiff->Fill(data_pt, pdiff,1);

	  p_angC->Fill(data_pt, cos(TH)*cos(TH), 1);
	  p_angS1->Fill(data_pt, sin(TH)*sin(TH)*cos(2*PHI), 1);
	  p_angS2->Fill(data_pt, sin(2*TH)*cos(PHI), 1);
	  for(int j = 0; j < nPtBins; j++) {
	    if(data_pt < ptBins[j+1] && data_pt > ptBins[j]) {	      
	      h_costh_pt[j]->Fill(cos(th));
	      h_cosTH_pt[j]->Fill(cos(TH));
	      h_phi_pt[j]->Fill(phi);
	      h_PHI_pt[j]->Fill(PHI);	      
	    }
	  }
	}
      }
    }
  
  TCanvas *c = new TCanvas("", "", 700, 700);
  double norm;
    
  for(int i = 0; i < nPtBins; i++) {
    h_cosTH_pt[i]->SetStats(0);
    h_cosTH_pt[i]->GetXaxis()->SetTitle("cos#theta");
    h_cosTH_pt[i]->SetLineColor(kBlack);
    h_cosTH_pt[i]->SetMarkerColor(kBlack);
    h_cosTH_pt[i]->SetMinimum(0);
    h_cosTH_pt[i]->Draw("hist");

    h_costh_pt[i]->SetStats(0);
    h_costh_pt[i]->GetXaxis()->SetTitle("cos#theta");
    h_costh_pt[i]->SetLineColor(kBlue);
    h_costh_pt[i]->SetMarkerColor(kBlue);
    h_costh_pt[i]->Draw("hist same");

    TLegend *leg_pt = new TLegend(0.2, 0.6, 0.4, 0.8);
    leg_pt->SetTextSize(0.03);
    leg_pt->AddEntry(h_cosTH_pt[i], "cos#Theta", "l");
    leg_pt->AddEntry(h_costh_pt[i], "cos#theta", "l");
    leg_pt->Draw();

    //c->SaveAs(Form("plots/costh_pt%d.pdf", i+1));
    c->Clear();
  }


  for(int i = 0; i < nPtBins; i++) {    
    h_phi_pt[i]->SetStats(0);
    h_phi_pt[i]->GetXaxis()->SetTitle("#phi");
    h_phi_pt[i]->SetLineColor(kBlue);
    h_phi_pt[i]->SetMarkerColor(kBlue);
    h_phi_pt[i]->SetMinimum(0);
    h_phi_pt[i]->Draw("hist");

    h_PHI_pt[i]->SetStats(0);
    h_PHI_pt[i]->GetXaxis()->SetTitle("#phi");
    h_PHI_pt[i]->SetLineColor(kBlack);
    h_PHI_pt[i]->SetMarkerColor(kBlack);
    h_PHI_pt[i]->Draw("hist same");

    TLegend *leg_pt = new TLegend(0.2, 0.2, 0.4, 0.4);
    leg_pt->SetTextSize(0.03);
    leg_pt->AddEntry(h_PHI_pt[i], "#Phi", "l");
    leg_pt->AddEntry(h_phi_pt[i], "#phi", "l");
    leg_pt->Draw();

    //c->SaveAs(Form("plots/phi_pt%d.pdf", i+1));
    c->Clear();
  }

  h_pdiff->SetStats(0);
  h_pdiff->GetXaxis()->SetTitle("#phi-#Phi");
  h_pdiff->SetMinimum(0);
  h_pdiff->Draw("hist");

  c->SaveAs("plots/phi_diff.pdf");
  c->Clear();

  p_pdiff->SetStats(0);
  p_pdiff->GetYaxis()->SetTitle("#phi-#Phi");
  p_pdiff->GetXaxis()->SetTitle("p_{T}");
  p_pdiff->Draw("error");

  c->SaveAs("plots/phi_diff_prof.pdf");
  c->Clear();

  p_angC->SetStats(0);
  p_angC->GetYaxis()->SetTitle("<cos^{2}#Theta>");
  p_angC->GetYaxis()->SetTitleOffset(1.2);
  p_angC->GetXaxis()->SetTitle("p_{T} (GeV)");
  p_angC->SetMinimum(-1);
  p_angC->SetMaximum(1);
  p_angC->Draw("error");

  c->SaveAs("plots/prof_cos.pdf");
  c->Clear();
  
  p_angS1->SetStats(0);
  p_angS1->GetYaxis()->SetTitle("<sin^{2}#Thetacos(2#Phi)>");
  p_angS1->GetYaxis()->SetTitleOffset(1.2);
  p_angS1->GetXaxis()->SetTitle("p_{T} (GeV)");
  p_angS1->SetMinimum(-1);
  p_angS1->SetMaximum(1);
  p_angS1->Draw("error");

  c->SaveAs("plots/prof_sin1.pdf");
  c->Clear();

  p_angS2->SetStats(0);
  p_angS2->GetYaxis()->SetTitle("<sin(2#Theta)cos#Phi>");
  p_angS2->GetYaxis()->SetTitleOffset(1.2);
  p_angS2->GetXaxis()->SetTitle("p_{T} (GeV)");
  p_angS2->SetMinimum(-1);
  p_angS2->SetMaximum(1);
  p_angS2->Draw("error");

  c->SaveAs("plots/prof_sin2.pdf");
  c->Clear();

  // calculate lambdas
  double ln[] = {-0.5, -0.75, -1};

  double xBins[nPtBins], exBins[nPtBins], zero[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    xBins[i] = 0.5*(ptBins[i]+ptBins[i+1]);
    exBins[i] = 0.5*(ptBins[i+1]-ptBins[i]);
    zero[i] = 0;
  }

  for(int i = 0; i < 3; i++) {
    double lambda_N = ln[i];
    double lambda_th[nPtBins], lambda_th_up[nPtBins], lambda_th_down[nPtBins];
    double lambda_ph[nPtBins], lambda_ph_up[nPtBins], lambda_ph_down[nPtBins];
    double lambda_tp[nPtBins], lambda_tp_up[nPtBins], lambda_tp_down[nPtBins];
    for(int j = 0; j < nPtBins; j++) {
      double cos = p_angC->GetBinContent(j+1);
      double e_cos = p_angC->GetBinError(j+1);
      
      double s1 = p_angS1->GetBinContent(j+1);
      double e_s1 = p_angS1->GetBinError(j+1);

      double s2 = p_angS2->GetBinContent(j+1);
      double e_s2 = p_angS2->GetBinError(j+1);

      lambda_th[j] = lambda_theta(cos, lambda_N);
      double e1 = lambda_theta(cos+e_cos, lambda_N);
      double e2 = lambda_theta(cos-e_cos, lambda_N);

      double max_val = max(e1, e2);
      double min_val = min(e1, e2);
      
      lambda_th_up[j] = max_val-lambda_th[j];
      lambda_th_down[j] = lambda_th[j]-min_val;

      lambda_ph[j] = lambda_phi(s1, cos, lambda_N);

      e1 = lambda_phi(s1+e_s1, cos+e_cos, lambda_N);
      e2 = lambda_phi(s1+e_s1, cos-e_cos, lambda_N);
      double e3 = lambda_phi(s1-e_s1, cos+e_cos, lambda_N);
      double e4 = lambda_phi(s1-e_s1, cos-e_cos, lambda_N);

      max_val = max(e1, e2);
      max_val = max(max_val, e3);
      max_val = max(max_val, e4);
      min_val = min(e1, e2);
      min_val = min(min_val, e3);
      min_val = min(min_val, e4);

      lambda_ph_up[j] = max_val - lambda_ph[j];
      lambda_ph_down[j] = lambda_ph[j] - min_val;
      
      lambda_tp[j] = lambda_thetaphi(s2, cos, lambda_N);

      e1 = lambda_thetaphi(s2+e_s2, cos+e_cos, lambda_N);
      e2 = lambda_thetaphi(s2+e_s2, cos-e_cos, lambda_N);
      e3 = lambda_thetaphi(s2-e_s2, cos+e_cos, lambda_N);
      e4 = lambda_thetaphi(s2-e_s2, cos-e_cos, lambda_N);

      max_val = max(e1, e2);
      max_val = max(max_val, e3);
      max_val = max(max_val, e4);
      min_val = min(e1, e2);
      min_val = min(min_val, e3);
      min_val = min(min_val, e4);

      lambda_tp_up[j] = max_val - lambda_tp[j];
      lambda_tp_down[j] = lambda_tp[j] - min_val;
}
    
    TGraphAsymmErrors *l_th = new TGraphAsymmErrors(nPtBins, xBins, lambda_th, exBins, exBins, lambda_th_down, lambda_th_up);
    TGraphAsymmErrors *l_ph = new TGraphAsymmErrors(nPtBins, xBins, lambda_ph, exBins, exBins, lambda_ph_down, lambda_ph_up);
    TGraphAsymmErrors *l_tp = new TGraphAsymmErrors(nPtBins, xBins, lambda_tp, exBins, exBins, lambda_tp_down, lambda_tp_up);

    TH1F *fc = c->DrawFrame(0, -1.0, 215, 1.0);
    fc->SetXTitle("p_{T}/M");
    fc->SetYTitle("#lambda");
    fc->GetYaxis()->SetTitleOffset(1);
    fc->SetTitle(Form("#lambda_{N} = %.2f", ln[i]));
    
    l_th->SetLineColor(kBlack);
    l_th->SetMarkerColor(kBlack);
    l_th->Draw("P");
    l_ph->SetLineColor(kBlue);
    l_ph->SetMarkerColor(kBlue);
    l_ph->Draw("P");
    l_tp->SetLineColor(kViolet);
    l_tp->SetMarkerColor(kViolet);
    l_tp->Draw("P");

    TLegend *leg_l = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_l->SetTextSize(0.03);
    leg_l->AddEntry(l_th, "#lambda_{#theta}", "pl");
    leg_l->AddEntry(l_ph, "#lambda_{#phi}", "pl");
    leg_l->AddEntry(l_tp, "#lambda_{#theta#phi}", "pl");
    leg_l->Draw();

    c->SaveAs(Form("plots/lambda_N_%d.pdf", i));
    c->Clear();
  }
  
  infile->Close();

}
