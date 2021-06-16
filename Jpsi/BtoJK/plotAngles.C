void plotAngles()
{
  TFile *infile = new TFile("../Store_data_codes/btopk_cos.root");
  TTree *tree = (TTree*)infile->Get("data_cos");
  
  int nEvt = tree->GetEntries();
  cout << nEvt << " data events after cuts" << endl;
  cout << tree->GetEntries("lts > 2.5") << " events with lts cut" << endl;
  cout << tree->GetEntries("lts > 2.5 && BMass > 5.24 && BMass < 5.32") << " events in signal region" << endl;

  // prepare binning and histograms for plots
  const int nPtBins = 21;
  double ptBins[nPtBins+1];
  for(int i=0; i<13; i++) ptBins[i] = 24.+2.*i;
  for(int i=0; i<4; i++) ptBins[i+13] = (50.+2.5*i);
  for(int i=0; i<2; i++) ptBins[i+17] = (60. + 5.*i);
  for(int i=0; i<1; i++) ptBins[i+19] = (70. + 10.*i);
  for(int i=0; i<2; i++) ptBins[i+20] = (80. + 20.*i);
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // making histograms
  // angular, dividing LSB, RSB, signal
  TH1D **h_costh = new TH1D*[3];
  TH1D **h_cosTH = new TH1D*[3];
  TH1D **h_phi   = new TH1D*[3];
  TH1D **h_PHI   = new TH1D*[3];
  for(int i = 0; i < 3; i++) {
    h_costh[i] = new TH1D(Form("h_costh_%d", i), "cos#theta (norm)", 200, -1.,   1.  );
    h_cosTH[i] = new TH1D(Form("h_cosTH_%d", i), "cos#Theta (norm)", 200, -1.,   1.  );
    h_phi[i]   = new TH1D(Form("h_phi_%d",   i), "#phi (norm)",      200, -180., 180.);
    h_PHI[i]   = new TH1D(Form("h_PHI_%d",   i), "#Phi (norm)",      200, -180., 180.);
  }
  // B mass and Jpsi pt
  TH1D *h_BM = new TH1D("h_BM", "B Mass (NP)", 200, 5, 5.6);
  TH1D *h_BM_nocut = new TH1D("h_BM_nocut", "B Mass (all)", 200, 5, 5.6);  
  TH1D *h_JPt = new TH1D("h_JPt", "J/#psi p_{T}", 200, 24, 100);
  // 2d of angular vs angular
  TH2D *h_cos = new TH2D("h_cos", "cos#theta:cos#Theta", 200, -1, 1, 200, -1, 1);
  TH2D *h_ph = new TH2D("h_ph", "#phi:#Phi", 200, -180, 180, 200, -180, 180);
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

  Double_t th, TH, phi, PHI, data_pt, data_mass, lts;

  tree->SetBranchAddress("th", &th);
  tree->SetBranchAddress("phi", &phi);
  tree->SetBranchAddress("TH", &TH);
  tree->SetBranchAddress("PHI", &PHI);
  tree->SetBranchAddress("JpsiPt", &data_pt);
  tree->SetBranchAddress("BMass", &data_mass);
  tree->SetBranchAddress("lts", &lts);

  double m_min[] = {5.17, 5.24, 5.36};
  double m_max[] = {5.2, 5.32, 5.5};

  double pdiff;
  //nEvt = 100;
  // cycle over data, fill the costh histogram acc to B mass
  for(int i = 0; i < nEvt; i++)
    {
      tree->GetEntry(i);
      h_BM_nocut->Fill(data_mass);
      if(lts > 2.5) {
	h_BM->Fill(data_mass);
	for(int j = 0; j < 3; j++) {
	  if(data_mass > m_min[j]  && data_mass < m_max[j]) {
	    h_costh[j]->Fill(cos(th));
	    h_cosTH[j]->Fill(cos(TH));
	    h_phi[j]->Fill(phi);
	    h_PHI[j]->Fill(PHI);
	  }
	}
	if(data_mass > m_min[1] && data_mass < m_max[1]) {
	  h_cos->Fill(cos(th), cos(TH));
	  h_ph->Fill(phi, PHI);
	  h_JPt->Fill(data_pt);
	  pdiff = phi-PHI;
	  if(pdiff > 180) pdiff -= 360;
	  if(pdiff < -180) pdiff += 360;

	  if(pdiff < -180 || pdiff > 180) cout << pdiff << endl;
	  
	  h_pdiff->Fill(pdiff);
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
  
  for(int i = 0; i < 3; i++) {
    norm = h_costh[i]->Integral();
    h_costh[i]->Scale(1./norm);
    h_costh[i]->SetStats(0);
    h_costh[i]->GetXaxis()->SetTitle("cos#theta");
    h_costh[i]->SetLineColor(i+1);
    h_costh[i]->SetMarkerColor(i+1);
    h_costh[i]->SetMinimum(0);
    h_costh[i]->SetMaximum(0.02);
    if(i == 0) h_costh[i]->Draw("hist");
    else h_costh[i]->Draw("hist same");
  }

  TLegend *leg = new TLegend(0.2, 0.6, 0.4, 0.8);
  leg->SetTextSize(0.03);
  for(int i = 0; i < 3; i++) {
    string name = i == 0 ? "LSB" : i == 1 ? "Signal" : "RSB";
    leg->AddEntry(h_costh[i], name.c_str(), "l");
  }
  leg->Draw();
  c->SaveAs("plots/costh.pdf");
  c->Clear();

  for(int i = 0; i < 3; i++) {
    norm = h_cosTH[i]->Integral();
    h_cosTH[i]->Scale(1./norm);
    h_cosTH[i]->SetStats(0);
    h_cosTH[i]->GetXaxis()->SetTitle("cos#Theta");
    h_cosTH[i]->SetLineColor(i+1);
    h_cosTH[i]->SetMarkerColor(i+1);
    h_cosTH[i]->SetMaximum(0.05);
    h_cosTH[i]->SetMinimum(0);
    if(i == 0) h_cosTH[i]->Draw("hist");
    else h_cosTH[i]->Draw("hist same");
  }

  leg->Draw();
  c->SaveAs("plots/cosTH.pdf");
  c->Clear();

  for(int i = 0; i < 3; i++) {
    norm = h_phi[i]->Integral();
    h_phi[i]->Scale(1./norm);
    h_phi[i]->SetStats(0);
    h_phi[i]->GetXaxis()->SetTitle("#phi");
    h_phi[i]->SetLineColor(i+1);
    h_phi[i]->SetMarkerColor(i+1);
    h_phi[i]->SetMaximum(0.008);
    h_phi[i]->SetMinimum(0);
    if(i == 0) h_phi[i]->Draw("hist");
    else h_phi[i]->Draw("hist same");
  }

  TLegend *leg_p = new TLegend(0.2, 0.2, 0.4, 0.4);
  leg_p->SetTextSize(0.03);
  for(int i = 0; i < 3; i++) {
    string name = i == 0 ? "LSB" : i == 1 ? "Signal" : "RSB";
    leg_p->AddEntry(h_phi[i], name.c_str(), "l");
  }
  leg_p->Draw();
  c->SaveAs("plots/phi.pdf");
  c->Clear();

  for(int i = 0; i < 3; i++) {
    norm = h_PHI[i]->Integral();
    h_PHI[i]->Scale(1./norm);
    h_PHI[i]->SetStats(0);
    h_PHI[i]->GetXaxis()->SetTitle("#Phi");
    h_PHI[i]->SetLineColor(i+1);
    h_PHI[i]->SetMarkerColor(i+1);
    h_PHI[i]->SetMaximum(0.008);
    h_PHI[i]->SetMinimum(0);
    if(i == 0) h_PHI[i]->Draw("hist");
    else h_PHI[i]->Draw("hist same");
  }

  leg_p->Draw();
  c->SaveAs("plots/PHI.pdf");
  c->Clear();

  h_cos->SetStats(0);
  h_cos->GetXaxis()->SetTitle("cos#theta");
  h_cos->GetYaxis()->SetTitle("cos#Theta");
  h_cos->GetYaxis()->SetTitleOffset(1.2);
  h_cos->Draw("COLZ");
  c->SaveAs("plots/cos_comp.pdf");

  h_ph->SetStats(0);
  h_ph->GetXaxis()->SetTitle("#phi");
  h_ph->GetYaxis()->SetTitle("#Phi");
  h_ph->GetYaxis()->SetTitleOffset(1.2);
  h_ph->Draw("COLZ");
  c->SaveAs("plots/phi_comp.pdf");

  h_BM->SetStats(0);
  h_BM->GetXaxis()->SetTitle("M^{B} (GeV)");
  h_BM->SetMinimum(0);
  h_BM->SetMaximum(5e4);
  h_BM->Draw("hist");
  
  TLine *LSB_1 = new TLine(m_max[0], h_BM->GetMinimum(), m_max[0], h_BM->GetMaximum());
  LSB_1->SetLineStyle(kDashed);
  LSB_1->SetLineColor(kRed);
  LSB_1->Draw();
  TLine *LSB_2 = new TLine(m_min[0], h_BM->GetMinimum(), m_min[0], h_BM->GetMaximum());
  LSB_2->SetLineStyle(kDashed);
  LSB_2->SetLineColor(kRed);
  LSB_2->Draw();
  TLine *RSB_1 = new TLine(m_min[2], h_BM->GetMinimum(), m_min[2], h_BM->GetMaximum());
  RSB_1->SetLineStyle(kDashed);
  RSB_1->SetLineColor(kRed);
  RSB_1->Draw();
  TLine *RSB_2 = new TLine(m_max[2], h_BM->GetMinimum(), m_max[2], h_BM->GetMaximum());
  RSB_2->SetLineStyle(kDashed);
  RSB_2->SetLineColor(kRed);
  RSB_2->Draw();
  TLine *Sig_1 = new TLine(m_min[1], h_BM->GetMinimum(), m_min[1], h_BM->GetMaximum());
  Sig_1->SetLineStyle(kDashed);
  Sig_1->Draw();
  TLine *Sig_2 = new TLine(m_max[1], h_BM->GetMinimum(), m_max[1], h_BM->GetMaximum());
  Sig_2->SetLineStyle(kDashed);
  Sig_2->Draw();

  c->SaveAs("plots/mass_cut.pdf");
  c->Clear();

  h_BM_nocut->SetStats(0);
  h_BM_nocut->GetXaxis()->SetTitle("M^{B} (GeV)");
  h_BM_nocut->SetMinimum(0);
  h_BM_nocut->SetMaximum(5e4);
  h_BM_nocut->Draw("hist");
  
  TLine *LSB_nocut_1 = new TLine(m_max[0], h_BM_nocut->GetMinimum(), m_max[0], h_BM_nocut->GetMaximum());
  LSB_nocut_1->SetLineStyle(kDashed);
  LSB_nocut_1->SetLineColor(kRed);  
  LSB_nocut_1->Draw();
  TLine *LSB_nocut_2 = new TLine(m_min[0], h_BM_nocut->GetMinimum(), m_min[0], h_BM_nocut->GetMaximum());
  LSB_nocut_2->SetLineStyle(kDashed);
  LSB_nocut_2->SetLineColor(kRed);  
  LSB_nocut_2->Draw();
  TLine *RSB_nocut_1 = new TLine(m_min[2], h_BM_nocut->GetMinimum(), m_min[2], h_BM_nocut->GetMaximum());
  RSB_nocut_1->SetLineStyle(kDashed);
  RSB_nocut_1->SetLineColor(kRed);
  RSB_nocut_1->Draw();
  TLine *RSB_nocut_2 = new TLine(m_max[2], h_BM_nocut->GetMinimum(), m_max[2], h_BM_nocut->GetMaximum());
  RSB_nocut_2->SetLineStyle(kDashed);
  RSB_nocut_2->SetLineColor(kRed);
  RSB_nocut_2->Draw();
  TLine *Sig_1_nocut = new TLine(m_min[1], h_BM_nocut->GetMinimum(), m_min[1], h_BM_nocut->GetMaximum());
  Sig_1_nocut->SetLineStyle(kDashed);
  Sig_1_nocut->Draw();
  TLine *Sig_2_nocut = new TLine(m_max[1], h_BM_nocut->GetMinimum(), m_max[1], h_BM_nocut->GetMaximum());
  Sig_2_nocut->SetLineStyle(kDashed);
  Sig_2_nocut->Draw();

  c->SaveAs("plots/mass_cut_all.pdf");
  c->Clear();

  c->SetLogy();
  h_JPt->SetStats(0);
  h_JPt->GetXaxis()->SetTitle("J/#psi p_{T} (GeV)");
  h_JPt->SetMinimum(3e1);
  h_JPt->SetMaximum(3e4);
  h_JPt->Draw("hist");
  c->SaveAs("plots/Jpsi_pt.pdf");
  c->Clear();
  c->SetLogy(0);
  
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

    c->SaveAs(Form("plots/costh_pt%d.pdf", i+1));
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

    c->SaveAs(Form("plots/phi_pt%d.pdf", i+1));
    c->Clear();
  }

  h_pdiff->SetStats(0);
  h_pdiff->GetXaxis()->SetTitle("#phi-#Phi");
  //h_pdiff->SetMaximum(0.008);
  h_pdiff->SetMinimum(0);
  h_pdiff->Draw("hist");

  c->SaveAs("plots/phi_diff.pdf");
  c->Clear();
  
  infile->Close();

}
