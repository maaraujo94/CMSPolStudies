// code to plot the cut variables for the MC
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y, ct/cterr
*/

void cutvarMC()
{
  // tree for MC
  TChain *tree = new TChain("jpsi_tuple");

  tree->Add("/eos/user/m/maaraujo/Jpsi2012/MC/flat_tuple_onia2MuMu_tree_validation_skim.root");

  // creating desired vars and setting branch address
  Double_t muPPt, muNPt, muPEta, muNEta, JpsiMass, JpsiPt, JpsiRap, Jpsict, JpsictErr;
  Long64_t trigger;

  tree->SetBranchAddress("muPPt", &muPPt);
  tree->SetBranchAddress("muNPt", &muNPt);
  tree->SetBranchAddress("muPEta", &muPEta);
  tree->SetBranchAddress("muNEta", &muNEta);
  tree->SetBranchAddress("JpsiMass", &JpsiMass);
  tree->SetBranchAddress("JpsiPt", &JpsiPt);
  tree->SetBranchAddress("JpsiRap", &JpsiRap);
  tree->SetBranchAddress("Jpsict", &Jpsict);
  tree->SetBranchAddress("JpsictErr", &JpsictErr);
  tree->SetBranchAddress("trigger", &trigger);

  // aux vars for reading tree
  int nEvt = tree->GetEntries();
  int perc = nEvt / 100;

  // preparing histograms to be filled
  TH1D *h_muPpT  = new TH1D("h_muPpT",  "MC muon pT",            100, 0, 70);
  TH1D *h_muNpT  = new TH1D("h_muNpT",  "muonN pT",              100, 0, 70);
  TH1D *h_muPEta = new TH1D("h_muPEta", "MC muon eta",           100, -3, 3);
  TH1D *h_muNEta = new TH1D("h_muNEta", "muonN eta",             100, -3, 3);
  TH1D *h_JMass  = new TH1D("h_JMass",  "MC Jpsi mass",          100, 2.5, 3.5);
  TH1D *h_JPt    = new TH1D("h_JPt",    "MC Jpsi pT",            100, 0, 90);
  TH1D *h_Jy     = new TH1D("h_Jy",     "MC Jpsi y",             100, -2.5, 2.5);
  TH1D *h_Jlts   = new TH1D("h_Jlts",   "MC Jpsi lifetime sig",  100, 0, 25);

  // reading tree, filling histograms
  for(int i = 0; i < nEvt; i++) {
    tree->GetEntry(i);

    if(trigger == 1)
      {
	h_muPpT->Fill(muPPt);
	h_muNpT->Fill(muNPt);
	h_muPEta->Fill(muPEta);
	h_muNEta->Fill(muNEta);
	h_JMass->Fill(JpsiMass);
	h_JPt->Fill(JpsiPt);
	h_Jy->Fill(JpsiRap);
	h_Jlts->Fill(abs(Jpsict/JpsictErr));
      }

    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with MC" << endl; 
  }

  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLogy();
  double norm;

  norm = h_muPpT->Integral();
  h_muPpT->Scale(1./norm);
  h_muNpT->Scale(1./norm);
  h_muPpT->GetYaxis()->SetRangeUser(1e-6, 3e-1);
  h_muPpT->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
  h_muPpT->SetLineColor(kRed);
  h_muPpT->Draw("hist");
  h_muNpT->SetLineColor(kBlue);
  h_muNpT->Draw("hist same");
  TLine *mupT = new TLine(6, 1e-6, 6, 3e-1);
  mupT->SetLineStyle(kDashed);
  mupT->SetLineColor(kBlack);
  mupT->Draw("lsame");
  c->SaveAs("plots/MC_muon_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_muPEta->Integral();
  h_muPEta->Scale(1./norm);
  h_muNEta->Scale(1./norm);
  h_muPEta->GetYaxis()->SetRangeUser(0, 0.022);
  h_muPEta->GetXaxis()->SetTitle("#eta(#mu)");
  h_muPEta->SetLineColor(kRed);
  h_muPEta->Draw("hist");
  h_muNEta->SetLineColor(kBlue);
  h_muNEta->Draw("hist same");
  TLine *muEta_1 = new TLine(-2, 0, -2, 0.022);
  muEta_1->SetLineStyle(kDashed);
  muEta_1->SetLineColor(kBlack);
  muEta_1->Draw("lsame");
  TLine *muEta_2 = new TLine(2, 0, 2, 0.022);
  muEta_2->SetLineStyle(kDashed);
  muEta_2->SetLineColor(kBlack);
  muEta_2->Draw("lsame");
  c->SaveAs("plots/MC_muon_eta.pdf");
  c->Clear();

  c->SetLogy();
  norm = h_JMass->Integral();
  h_JMass->Scale(1./norm);
  h_JMass->GetYaxis()->SetRangeUser(1e-6, 3e-1);
  h_JMass->GetXaxis()->SetTitle("M(J/#psi) (GeV)");
  h_JMass->Draw("hist");
  TLine *JMass_1 = new TLine(3, 1e-6, 3, 3e-1);
  JMass_1->SetLineStyle(kDashed);
  JMass_1->SetLineColor(kBlack);
  JMass_1->Draw("lsame");
  TLine *JMass_2 = new TLine(3.2, 1e-6, 3.2, 3e-1);
  JMass_2->SetLineStyle(kDashed);
  JMass_2->SetLineColor(kBlack);
  JMass_2->Draw("lsame");
  c->SaveAs("plots/MC_jpsi_mass.pdf");
  c->Clear();

  norm = h_JPt->Integral();
  h_JPt->Scale(1./norm);
  h_JPt->GetYaxis()->SetRangeUser(1e-6, 3e-1);
  h_JPt->GetXaxis()->SetTitle("p_{T}(J/#psi) (GeV)");
  h_JPt->Draw("hist");
  TLine *JPt_1 = new TLine(12, 1e-6, 12, 3e-1);
  JPt_1->SetLineStyle(kDashed);
  JPt_1->SetLineColor(kBlack);
  JPt_1->Draw("lsame");
  TLine *JPt_2 = new TLine(70, 1e-6, 70, 3e-1);
  JPt_2->SetLineStyle(kDashed);
  JPt_2->SetLineColor(kBlack);
  JPt_2->Draw("lsame");
  c->SaveAs("plots/MC_jpsi_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_Jy->Integral();
  h_Jy->Scale(1./norm);
   h_Jy->GetYaxis()->SetRangeUser(0, 0.019);
  h_Jy->GetXaxis()->SetTitle("y(J/#psi)");
  h_Jy->Draw("hist");
  TLine *JRap_1 = new TLine(-1.5, 0, -1.5, 0.019);
  JRap_1->SetLineStyle(kDashed);
  JRap_1->SetLineColor(kBlack);
  JRap_1->Draw("lsame");
  TLine *JRap_2 = new TLine(1.5, 0, 1.5, 0.019);
  JRap_2->SetLineStyle(kDashed);
  JRap_2->SetLineColor(kBlack);
  JRap_2->Draw("lsame");
  c->SaveAs("plots/MC_jpsi_rap.pdf");
  c->Clear();

  c->SetLogy();
  norm = h_Jlts->Integral();
  h_Jlts->Scale(1./norm);
  h_Jlts->GetYaxis()->SetRangeUser(1e-4, 3e-1);
  h_Jlts->GetXaxis()->SetTitle("|c#tau|/#sigma_{c#tau}");
  h_Jlts->Draw("hist");
  TLine *Jlts = new TLine(2, 1e-4, 2, 3e-1);
  Jlts->SetLineStyle(kDashed);
  Jlts->SetLineColor(kBlack);
  Jlts->Draw("lsame");
  c->SaveAs("plots/MC_jpsi_lts.pdf");

}
