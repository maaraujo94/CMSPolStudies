// code to plot the cut variables for the Data and MC Psi(2S) sample
/* variables to plot
   - single muon pT, eta
   - dimuon mass, pT, y, ct/cterr
   - vProb
*/

void plotCutVar()
{
  // each of the sources of cut vars
  string s_source[3] = {"data", "MC", "MC_hpt"};
  string s_name[3] = {"Data", "MC (low p_{T})", "MC (high p_{T})"};

  // cut var ranges and cuts
  double mup_range[2] = {1e1, 2e6};
  double mup_cut = 5.6;
  double mue_range[2] = {0, 1.5e5};
  double mue_cut = 1.4;
  double mass_range[2] = {0, 3.5e5};
  double mass_cut[2] = {3.565, 3.805};
  double dmp_range[2] = {2e1, 5e6};
  double dmp_cut[3][2] = {{25, 100}, {25, 46}, {46, 100}};
  double dmy_range[2] = {0, 1.5e5};
  double dmy_cut = 1.2;
  double dml_range[2] = {1e4, 1e6};
  double dml_cut[2] = {2.5, 4};
  double vp_range[2] = {1e4, 2e5};
  double vp_cut = 0.01;
  
  // cycle over sources
  for(int i_i = 0; i_i < 3; i_i++) {
    // 2017 plots
    TFile *fin7 = new TFile(Form("../Store_data_codes/Psi2_17_%s_cuts.root", s_source[i_i].c_str()));
    TH1D *h7_muPpT = (TH1D*)fin7->Get("h7_muPpT");
    TH1D *h7_muNpT = (TH1D*)fin7->Get("h7_muNpT");
    TH1D *h7_muPEta = (TH1D*)fin7->Get("h7_muPEta");
    TH1D *h7_muNEta = (TH1D*)fin7->Get("h7_muNEta");
    TH1D *h7_JMass = (TH1D*)fin7->Get("h7_JMass");
    TH1D *h7_JPt = (TH1D*)fin7->Get("h7_JPt");
    TH1D *h7_Jy = (TH1D*)fin7->Get("h7_Jy");
    TH1D *h7_Jlts = (TH1D*)fin7->Get("h7_Jlts");
    TH1D *h7_vP;
    if(i_i == 0)
      fin7->GetObject("h7_vP", h7_vP);
  
    TCanvas *c = new TCanvas("", "", 700, 700);
    c->SetLogy();

    h7_muPpT->GetYaxis()->SetRangeUser(mup_range[0], mup_range[1]);
    h7_muPpT->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
    h7_muPpT->SetLineColor(kRed);
    h7_muPpT->SetMarkerColor(kRed);
    h7_muPpT->SetTitle(Form("2017 %s muon p_{T}", s_name[i_i].c_str()));
    h7_muPpT->Draw("error");
    h7_muNpT->SetLineColor(kBlue);
    h7_muNpT->SetMarkerColor(kBlue);
    h7_muNpT->Draw("error same");
    TLine *mupT = new TLine(mup_cut, mup_range[0], mup_cut, mup_range[1]);
    mupT->SetLineStyle(kDashed);
    mupT->SetLineColor(kBlack);
    mupT->Draw("lsame");
    c->SaveAs(Form("plots/%s_17_muon_pt.pdf", s_source[i_i].c_str()));
    c->Clear();

    c->SetLogy(0);
    h7_muPEta->GetYaxis()->SetRangeUser(mue_range[0], mue_range[1]);
    h7_muPEta->GetXaxis()->SetTitle("#eta(#mu)");
    h7_muPEta->SetLineColor(kRed);
    h7_muPEta->SetMarkerColor(kRed);
    h7_muPEta->SetTitle(Form("2017 %s muon #eta", s_name[i_i].c_str()));
    h7_muPEta->Draw("error");
    h7_muNEta->SetLineColor(kBlue);
    h7_muNEta->SetMarkerColor(kBlue);
    h7_muNEta->Draw("error same");
    TLine *muEta_1 = new TLine(-mue_cut, mue_range[0], -mue_cut, mue_range[1]);
    muEta_1->SetLineStyle(kDashed);
    muEta_1->SetLineColor(kBlack);
    muEta_1->Draw("lsame");
    TLine *muEta_2 = new TLine(mue_cut, mue_range[0], mue_cut, mue_range[1]);
    muEta_2->SetLineStyle(kDashed);
    muEta_2->SetLineColor(kBlack);
    muEta_2->Draw("lsame");
    c->SaveAs(Form("plots/%s_17_muon_eta.pdf", s_source[i_i].c_str()));
    c->Clear();

    c->SetLogy(0);
    h7_JMass->GetYaxis()->SetRangeUser(mass_range[0], mass_range[1]);
    h7_JMass->GetXaxis()->SetTitle("M(#psi(2S)) (GeV)");
    h7_JMass->SetTitle(Form("2017 %s #psi(2S) mass", s_name[i_i].c_str()));
    h7_JMass->Draw("error");
    TLine *JMass_1 = new TLine(mass_cut[0], mass_range[0], mass_cut[0], mass_range[1]);
    JMass_1->SetLineStyle(kDashed);
    JMass_1->SetLineColor(kBlack);
    JMass_1->Draw("lsame");
    TLine *JMass_2 = new TLine(mass_cut[1], mass_range[0], mass_cut[1], mass_range[1]);
    JMass_2->SetLineStyle(kDashed);
    JMass_2->SetLineColor(kBlack);
    JMass_2->Draw("lsame");
    c->SaveAs(Form("plots/%s_17_psi2_mass.pdf", s_source[i_i].c_str()));
    c->Clear();

    c->SetLogy();
    h7_JPt->GetYaxis()->SetRangeUser(dmp_range[0], dmp_range[1]);
    h7_JPt->GetXaxis()->SetTitle("p_{T}(#psi(2S)) (GeV)");
    h7_JPt->SetTitle(Form("2017 %s #psi(2S) p_{T}", s_name[i_i].c_str()));
    h7_JPt->Draw("error");
    TLine *JPt_1 = new TLine(dmp_cut[i_i][0], dmp_range[0], dmp_cut[i_i][0], dmp_range[1]);
    JPt_1->SetLineStyle(kDashed);
    JPt_1->SetLineColor(kBlack);
    JPt_1->Draw("lsame");
    TLine *JPt_2 = new TLine(dmp_cut[i_i][1], dmp_range[0], dmp_cut[i_i][1], dmp_range[1]);
    JPt_2->SetLineStyle(kDashed);
    JPt_2->SetLineColor(kBlack);
    JPt_2->Draw("lsame");
    c->SaveAs(Form("plots/%s_17_psi2_pt.pdf", s_source[i_i].c_str()));
    c->Clear();

    c->SetLogy(0);
    h7_Jy->GetYaxis()->SetRangeUser(dmy_range[0], dmy_range[1]);
    h7_Jy->GetXaxis()->SetTitle("y(#psi(2S))");
    h7_Jy->SetTitle(Form("2017 %s #psi(2S) y", s_name[i_i].c_str()));
    h7_Jy->Draw("error");
    TLine *JRap_1 = new TLine(-dmy_cut, dmy_range[0], -dmy_cut, dmy_range[1]);
    JRap_1->SetLineStyle(kDashed);
    JRap_1->SetLineColor(kBlack);
    JRap_1->Draw("lsame");
    TLine *JRap_2 = new TLine(dmy_cut, dmy_range[0], dmy_cut, dmy_range[1]);
    JRap_2->SetLineStyle(kDashed);
    JRap_2->SetLineColor(kBlack);
    JRap_2->Draw("lsame");
    c->SaveAs(Form("plots/%s_17_psi2_rap.pdf", s_source[i_i].c_str()));
    c->Clear();

    c->SetLogy();
    h7_Jlts->GetYaxis()->SetRangeUser(dml_range[0], dml_range[1]);
    h7_Jlts->GetXaxis()->SetTitle("|c#tau|/#sigma_{c#tau}");
    h7_Jlts->SetTitle(Form("2017 %s #psi(2S) lts", s_name[i_i].c_str()));
    h7_Jlts->Draw("error");
    TLine *Jlts_0 = new TLine(-dml_cut[0], dml_range[0], -dml_cut[0], dml_range[1]);
    Jlts_0->SetLineStyle(kDashed);
    Jlts_0->SetLineColor(kBlack);
    Jlts_0->Draw("lsame");
    TLine *Jlts_1 = new TLine(dml_cut[0], dml_range[0], dml_cut[0], dml_range[1]);
    Jlts_1->SetLineStyle(kDashed);
    Jlts_1->SetLineColor(kBlack);
    Jlts_1->Draw("lsame");
    TLine *Jlts_2 = new TLine(dml_cut[1], dml_range[0], dml_cut[1], dml_range[1]);
    Jlts_2->SetLineStyle(kDashed);
    Jlts_2->SetLineColor(kBlack);
    Jlts_2->Draw("lsame");
    c->SaveAs(Form("plots/%s_17_psi2_lts.pdf", s_source[i_i].c_str()));
    c->Clear();

    TLine *vP_0 = new TLine(vp_cut, vp_range[0], vp_cut, vp_range[1]);
    if(i_i == 0) {
      c->SetLogy(0);
      h7_vP->GetYaxis()->SetRangeUser(vp_range[0], vp_range[1]);
      h7_vP->GetXaxis()->SetTitle("vProb");
      h7_vP->SetTitle(Form("2017 %s #psi(2S) vProb", s_name[i_i].c_str()));
      h7_vP->Draw("error");
      vP_0->SetLineStyle(kDashed);
      vP_0->SetLineColor(kBlack);
      vP_0->Draw("lsame");
      c->SaveAs(Form("plots/%s_17_psi2_vP.pdf", s_source[i_i].c_str()));
    }
    
    fin7->Close();

    // 2018 plots
    TFile *fin8 = new TFile(Form("../Store_data_codes/Psi2_18_%s_cuts.root", s_source[i_i].c_str()));
    TH1D *h8_muPpT = (TH1D*)fin8->Get("h8_muPpT");
    TH1D *h8_muNpT = (TH1D*)fin8->Get("h8_muNpT");
    TH1D *h8_muPEta = (TH1D*)fin8->Get("h8_muPEta");
    TH1D *h8_muNEta = (TH1D*)fin8->Get("h8_muNEta");
    TH1D *h8_JMass = (TH1D*)fin8->Get("h8_JMass");
    TH1D *h8_JPt = (TH1D*)fin8->Get("h8_JPt");
    TH1D *h8_Jy = (TH1D*)fin8->Get("h8_Jy");
    TH1D *h8_Jlts = (TH1D*)fin8->Get("h8_Jlts");
    TH1D *h8_vP;
    if(i_i == 0)
      fin8->GetObject("h8_vP", h8_vP);
  
    c->SetLogy();
    
    h8_muPpT->GetYaxis()->SetRangeUser(mup_range[0], mup_range[1]);
    h8_muPpT->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
    h8_muPpT->SetLineColor(kRed);
    h8_muPpT->SetMarkerColor(kRed);
    h8_muPpT->SetTitle(Form("2018 %s muon p_{T}", s_name[i_i].c_str()));
    h8_muPpT->Draw("error");
    h8_muNpT->SetLineColor(kBlue);
    h8_muNpT->SetMarkerColor(kBlue);
    h8_muNpT->Draw("error same");
    mupT->Draw("lsame");
    c->SaveAs(Form("plots/%s_18_muon_pt.pdf", s_source[i_i].c_str()));
    c->Clear();

    c->SetLogy(0);
    h8_muPEta->GetYaxis()->SetRangeUser(mue_range[0], mue_range[1]);
    h8_muPEta->GetXaxis()->SetTitle("#eta(#mu)");
    h8_muPEta->SetLineColor(kRed);
    h8_muPEta->SetMarkerColor(kRed);
    h8_muPEta->SetTitle(Form("2018 %s muon #eta", s_name[i_i].c_str()));
    h8_muPEta->Draw("error");
    h8_muNEta->SetLineColor(kBlue);
    h8_muNEta->SetMarkerColor(kBlue);
    h8_muNEta->Draw("error same");
    muEta_1->Draw("lsame");
    muEta_2->Draw("lsame");
    c->SaveAs(Form("plots/%s_18_muon_eta.pdf", s_source[i_i].c_str()));
    c->Clear();

    c->SetLogy(0);
    h8_JMass->GetYaxis()->SetRangeUser(mass_range[0], mass_range[1]);
    h8_JMass->GetXaxis()->SetTitle("M(#psi(2S)) (GeV)");
    h8_JMass->SetTitle(Form("2018 %s #psi(2S) mass", s_name[i_i].c_str()));
    h8_JMass->Draw("error");
    JMass_1->Draw("lsame");
    JMass_2->Draw("lsame");
    c->SaveAs(Form("plots/%s_18_psi2_mass.pdf", s_source[i_i].c_str()));
    c->Clear();

    c->SetLogy();
    h8_JPt->GetYaxis()->SetRangeUser(dmp_range[0], dmp_range[1]);
    h8_JPt->GetXaxis()->SetTitle("p_{T}(#psi(2S)) (GeV)");
    h8_JPt->SetTitle(Form("2018 %s #psi(2S) p_{T}", s_name[i_i].c_str()));
    h8_JPt->Draw("error");
    JPt_1->Draw("lsame");
    JPt_2->Draw("lsame");
    c->SaveAs(Form("plots/%s_18_psi2_pt.pdf", s_source[i_i].c_str()));
    c->Clear();

    c->SetLogy(0);
    h8_Jy->GetYaxis()->SetRangeUser(dmy_range[0], dmy_range[1]);
    h8_Jy->GetXaxis()->SetTitle("y(#psi(2S))");
    h8_Jy->SetTitle(Form("2018 %s #psi(2S) y", s_name[i_i].c_str()));
    h8_Jy->Draw("error");
    JRap_1->Draw("lsame");
    JRap_2->Draw("lsame");
    c->SaveAs(Form("plots/%s_18_psi2_rap.pdf", s_source[i_i].c_str()));
    c->Clear();

    c->SetLogy();
    h8_Jlts->GetYaxis()->SetRangeUser(dml_range[0], dml_range[1]);
    h8_Jlts->GetXaxis()->SetTitle("|c#tau|/#sigma_{c#tau}");
    h8_Jlts->SetTitle(Form("2018 %s #psi(2S) lts", s_name[i_i].c_str()));
    h8_Jlts->Draw("error");
    Jlts_0->Draw("lsame");
    Jlts_1->Draw("lsame");
    Jlts_2->Draw("lsame");
    c->SaveAs(Form("plots/%s_18_psi2_lts.pdf", s_source[i_i].c_str()));
    c->Clear();

    if(i_i == 0) {
      c->SetLogy(0);
      h8_vP->GetYaxis()->SetRangeUser(vp_range[0], vp_range[1]);
      h8_vP->GetXaxis()->SetTitle("vProb");
      h8_vP->SetTitle(Form("2018 %s #psi(2S) vProb", s_name[i_i].c_str()));
      h8_vP->Draw("error");
      vP_0->Draw("lsame");
      c->SaveAs(Form("plots/%s_18_psi2_vP.pdf", s_source[i_i].c_str()));
    }
    
    fin8->Close();

    c->Destructor();
   }
}
