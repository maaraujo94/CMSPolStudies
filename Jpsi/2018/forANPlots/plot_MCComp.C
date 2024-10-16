// macro to take the stored histos and plot nicely

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


void plot_MCComp()
{
  // PART 1 : reading the histograms

  // 6 dists of everything
  TH1D **h_pT = new TH1D*[6]; 
  TH1D **h_y = new TH1D*[6];
  TH1D **h_m = new TH1D*[6]; 
  TH1D **hP_pT = new TH1D*[6]; 
  TH1D **hM_pT = new TH1D*[6]; 
  TH1D **hP_eta = new TH1D*[6]; 
  TH1D **hM_eta = new TH1D*[6];

  double pT_i[] = {41, 47, 67};
  double pT_f[] = {45, 51, 120};
  
  TFile *fin = new TFile("files/store_MCdists.root");
 
  // read the dists
  string lbl_s[] = {"lowPtMC0", "midPtMC0",
		    "midPtMC1", "highpTMC1",
		    "highpTMC2", "vhighpTMC2"};
  for(int i = 0; i < 6; i++) {
    h_pT[i] = (TH1D*)fin->Get(Form("h_pT_%s", lbl_s[i].c_str()));
    h_y[i] = (TH1D*)fin->Get(Form("h_y_%s", lbl_s[i].c_str()));
    h_m[i] = (TH1D*)fin->Get(Form("h_m_%s", lbl_s[i].c_str()));
    hP_pT[i] = (TH1D*)fin->Get(Form("hP_pT_%s", lbl_s[i].c_str()));
    hM_pT[i] = (TH1D*)fin->Get(Form("hM_pT_%s", lbl_s[i].c_str()));
    hP_eta[i] = (TH1D*)fin->Get(Form("hP_eta_%s", lbl_s[i].c_str()));
    hM_eta[i] = (TH1D*)fin->Get(Form("hM_eta_%s", lbl_s[i].c_str()));
  }

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  
  // plot pT
  c->SetLogy();
  c->SetLeftMargin(0.11);
  int colpt[] = {kBlack, kRed+1, kViolet+2, kGreen+3, kBlue};
		
  for(int i = 0; i < 6; i++) {

    h_pT[i]->SetStats(0);
    h_pT[i]->SetLineColor(colpt[i/2+1]);
    h_pT[i]->SetMinimum(1.5e3);
    h_pT[i]->SetMaximum(7e5);
    h_pT[i]->GetXaxis()->SetRangeUser(35,125);
    h_pT[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV)");
    h_pT[i]->GetXaxis()->SetTitleOffset(1.2);
    h_pT[i]->GetXaxis()->CenterTitle(true);
    h_pT[i]->GetXaxis()->SetLabelOffset(0.01);
    h_pT[i]->GetYaxis()->SetLabelOffset(0.01);
    h_pT[i]->GetYaxis()->SetTitle("dN/d#it{p}_{T}");
    h_pT[i]->SetTitle("");
    if(i == 0) {
      h_pT[i]->Draw("histo");
    }
    else if(i%2 == 1) {
      h_pT[i]->SetLineStyle(kDashed);
      h_pT[i]->Scale(h_pT[i-1]->Integral()/h_pT[i]->Integral());
      h_pT[i]->Draw("histo same");
    }
    else {
      h_pT[i]->Draw("histo same");
    }
  }

  h_m[0]->SetLineColor(kBlack);
  h_m[1]->SetLineColor(kBlack);
  h_m[1]->SetLineStyle(kDashed);

  TLegend *legcpt = new TLegend(0.575, 0.85, 0.875, 0.975);
  legcpt->SetTextSize(0.04);
  legcpt->SetBorderSize(0);
  legcpt->SetFillColorAlpha(kWhite,0);
  legcpt->AddEntry(h_m[0], "lower-#it{p}_{T} MC", "l");
  legcpt->AddEntry(h_m[1], "higher-#it{p}_{T} MC", "l");
  legcpt->Draw();

  TLatex lcpt;
  double yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.91, 1);
  lcpt.SetTextSize(0.04);
  lcpt.DrawLatex(45, yl, "#bf{2018 J/#psi}");
  lcpt.SetTextColor(colpt[1]);
  yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.78, 1);
  lcpt.DrawLatex(90, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[0], pT_f[0]));
  lcpt.SetTextColor(colpt[2]);
  yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.7, 1);
  lcpt.DrawLatex(90, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[1], pT_f[1]));
  lcpt.SetTextColor(colpt[3]);
  yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.62, 1);
  lcpt.DrawLatex(90, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[2], pT_f[2]));

  c->SaveAs(Form("plots/MCdists/pt_scale.pdf"));
  c->Clear();

  // plotting y
  for(int i = 0; i < 6; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(colpt[i/2+1]);
    h_y[i]->SetMinimum(2e3);
    h_y[i]->SetMaximum(4e4);
    h_y[i]->GetXaxis()->SetTitle("#it{y}");
    h_y[i]->GetXaxis()->SetTitleOffset(1.2);
    h_y[i]->GetXaxis()->CenterTitle(true);
    h_y[i]->GetXaxis()->SetLabelOffset(0.01);
    h_y[i]->GetYaxis()->SetLabelOffset(0.01);
    h_y[i]->GetYaxis()->SetTitle("dN/d#it{y}");
    h_y[i]->SetTitle("");
    if(i == 0) {
      h_y[i]->Draw("histo");
    }
    else if(i%2 == 1) {
      h_y[i]->SetLineStyle(kDashed);
      h_y[i]->Scale(h_y[i-1]->Integral()/h_y[i]->Integral());
      h_y[i]->Draw("histo same");
    }
    else {
      h_y[i]->Draw("histo same");
    }
  }

  TLegend *legcy = new TLegend(0.55, 0.35, 0.85, 0.475);
  legcy->SetTextSize(0.04);
  legcy->SetBorderSize(0);
  legcy->SetFillColorAlpha(kWhite,0);
  legcy->AddEntry(h_m[0], "lower-#it{p}_{T} MC", "l");
  legcy->AddEntry(h_m[1], "higher-#it{p}_{T} MC", "l");
  legcy->Draw();
  
  TLatex lcy;
  lcy.SetTextSize(0.04);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.91, 1);
  lcy.DrawLatex(-1.1, yl, "#bf{2018 J/#psi}");
  lcy.SetTextColor(colpt[1]);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.4, 1);
  lcy.DrawLatex(-1, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[0], pT_f[0]));
  lcy.SetTextColor(colpt[2]);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.32, 1);
  lcy.DrawLatex(-1, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[1], pT_f[1]));
  lcy.SetTextColor(colpt[3]);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.24, 1);
  lcy.DrawLatex(-1, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[2], pT_f[2]));

  c->SaveAs(Form("plots/MCdists/y_scale.pdf"));
  c->Clear();

  // y ratio data/MC
  c->SetLogy(0);
  double sc_fc[] = {-0.2, 0, 0.2};
  TH1D **h_yr = new TH1D*[3]; 
  for(int i = 0; i < 3; i++) {
    h_yr[i] = (TH1D*)h_y[2*i]->Clone(Form("rH_%d", i));
    h_yr[i]->Sumw2();
    h_yr[i]->Divide(h_y[2*i+1]);
    // scaling up and down
    for(int j = 0; j < h_yr[i]->GetNbinsX(); j++) {
      h_yr[i]->SetBinContent(j+1, h_yr[i]->GetBinContent(j+1)+sc_fc[i]);
    }
    
    h_yr[i]->SetStats(0);
    h_yr[i]->SetLineColor(colpt[i+1]);
    h_yr[i]->SetMinimum(0.51);
    h_yr[i]->SetMaximum(1.49);
    h_yr[i]->GetXaxis()->SetTitle("#it{y}");
    h_yr[i]->GetXaxis()->SetTitleOffset(1.2);
    h_yr[i]->GetXaxis()->CenterTitle(true);
    h_yr[i]->GetXaxis()->SetLabelOffset(0.01);
    h_yr[i]->GetYaxis()->SetLabelOffset(0.01);
    h_yr[i]->GetYaxis()->SetTitle("MC ratio");
    h_yr[i]->SetTitle("");
    if(i == 0) {
      h_yr[i]->Draw("histo");
    }
    else {
      h_yr[i]->Draw("histo same");
    }
  }
  
  TLatex lcyr;
  lcyr.SetTextSize(0.04);
  yl = getPos(h_yr[0]->GetMinimum(), h_yr[0]->GetMaximum(), 0.91, 0);
  lcyr.DrawLatex(-1.1, yl, "#bf{2018 J/#psi}");
  lcyr.SetTextColor(colpt[1]);
  yl = getPos(h_yr[0]->GetMinimum(), h_yr[0]->GetMaximum(), 0.92, 0);
  lcyr.DrawLatex(0.45, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[0], pT_f[0]));
  lcyr.SetTextColor(colpt[2]);
  yl = getPos(h_yr[0]->GetMinimum(), h_yr[0]->GetMaximum(), 0.84, 0);
  lcyr.DrawLatex(0.45, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[1], pT_f[1]));
  lcyr.SetTextColor(colpt[3]);
  yl = getPos(h_yr[0]->GetMinimum(), h_yr[0]->GetMaximum(), 0.76, 0);
  lcyr.DrawLatex(0.45, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[2], pT_f[2]));

  TLine **ly = new TLine*[3];
  for(int i = 0; i < 3; i++) {
    ly[i] = new TLine(-1.5, 1+sc_fc[i], 1.5, 1+sc_fc[i]);
    ly[i]->SetLineStyle(kDashed);
    ly[i]->SetLineColor(colpt[i+1]);
    ly[i]->Draw();
  }
  
  c->SaveAs(Form("plots/MCdists/y_ratio.pdf"));
  c->Clear();

  // plotting m
  c->SetLogy(0);
  c->SetLeftMargin(0.12);
  for(int i = 0; i < 6; i++) {

    h_m[i]->SetStats(0);
    h_m[i]->SetLineColor(colpt[i/2+1]);
    h_m[i]->Scale(1./h_m[i]->Integral());
    h_m[i]->SetMinimum(0);
    h_m[i]->SetMaximum(0.085);
    h_m[i]->GetXaxis()->SetTitle("#it{m} (GeV)");
    h_m[i]->GetXaxis()->SetTitleOffset(1.2);
    h_m[i]->GetXaxis()->CenterTitle(true);
    h_m[i]->GetXaxis()->SetLabelOffset(0.01);
    h_m[i]->GetYaxis()->SetLabelOffset(0.01);
    h_m[i]->GetYaxis()->SetTitle("dN/d#it{m} (normalized)");
    h_m[i]->SetTitle("");
    if(i == 0) {
      h_m[i]->Draw("histo");
    }
    else if(i%2 == 1) {
      h_m[i]->SetLineStyle(kDashed);
      h_m[i]->Scale(h_m[i-1]->Integral()/h_m[i]->Integral());
      h_m[i]->Draw("histo same");
    }
    else {
      h_m[i]->Draw("histo same");
    }
  }

  h_pT[0]->SetLineColor(kBlack);
  h_pT[1]->SetLineColor(kBlack);
  h_pT[1]->SetLineStyle(kDashed);

  TLegend *legcm = new TLegend(0.625, 0.825, 0.925, 0.95);
  legcm->SetTextSize(0.04);
  legcm->SetBorderSize(0);
  legcm->SetFillColorAlpha(kWhite,0);
  legcm->AddEntry(h_pT[0], "lower-#it{p}_{T} MC", "l");
  legcm->AddEntry(h_pT[1], "higher-#it{p}_{T} MC", "l");
  legcm->Draw();

  TLatex lcm;
  lcm.SetTextSize(0.04);
  yl = getPos(h_m[0]->GetMinimum(), h_m[0]->GetMaximum(), 0.91, 0);
  lcm.DrawLatex(2.95, yl, "#bf{2018 J/#psi}");
  lcm.SetTextColor(colpt[1]);
  yl = getPos(h_m[0]->GetMinimum(), h_m[0]->GetMaximum(), 0.6, 0);
  lcm.DrawLatex(3.175, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[0], pT_f[0]));
  lcm.SetTextColor(colpt[2]);
  yl = getPos(h_m[0]->GetMinimum(), h_m[0]->GetMaximum(), 0.52, 0);
  lcm.DrawLatex(3.175, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[1], pT_f[1]));
  lcm.SetTextColor(colpt[3]);
  yl = getPos(h_m[0]->GetMinimum(), h_m[0]->GetMaximum(), 0.44, 0);
  lcm.DrawLatex(3.175, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[2], pT_f[2]));

  c->SaveAs(Form("plots/MCdists/m_scale.pdf"));
  c->Clear();

  // plotting mu+ pT
  c->SetLogy();
  c->SetLeftMargin(0.11);
  for(int i = 0; i < 6; i++) {

    hP_pT[i]->SetStats(0);
    hP_pT[i]->SetLineColor(colpt[i/2+1]);
    hP_pT[i]->SetMinimum(1e2);
    hP_pT[i]->SetMaximum(5e4);
    hP_pT[i]->GetXaxis()->SetRangeUser(0,125);
    hP_pT[i]->GetXaxis()->SetTitle("#it{p}_{T}(#mu^{+}) (GeV)");
    hP_pT[i]->GetXaxis()->SetTitleOffset(1.2);
    hP_pT[i]->GetXaxis()->CenterTitle(true);
    hP_pT[i]->GetXaxis()->SetLabelOffset(0.01);
    hP_pT[i]->GetYaxis()->SetLabelOffset(0.01);
    hP_pT[i]->GetYaxis()->SetTitle("dN/d#it{p}_{T}");
    hP_pT[i]->SetTitle("");
    if(i == 0) {
      hP_pT[i]->Draw("histo");
    }
    else if(i%2 == 1) {
      hP_pT[i]->SetLineStyle(kDashed);
      hP_pT[i]->Scale(hP_pT[i-1]->Integral()/hP_pT[i]->Integral());
      hP_pT[i]->Draw("histo same");
    }
    else {
      hP_pT[i]->Draw("histo same");
    }
  }

  TLegend *legcPPt = new TLegend(0.65, 0.85, 0.95, 0.975);
  legcPPt->SetTextSize(0.04);
  legcPPt->SetBorderSize(0);
  legcPPt->SetFillColorAlpha(kWhite,0);
  legcPPt->AddEntry(h_pT[0], "lower-#it{p}_{T} MC", "l");
  legcPPt->AddEntry(h_pT[1], "higher-#it{p}_{T} MC", "l");
  legcPPt->Draw();
  
  TLatex lcPPt;
  lcPPt.SetTextSize(0.04);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.925, 1);
  lcPPt.DrawLatex(45, yl, "#bf{2018 J/#psi}");
  lcPPt.SetTextColor(colpt[1]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.76, 1);
  lcPPt.DrawLatex(90, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[0], pT_f[0]));
  lcPPt.SetTextColor(colpt[2]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.68, 1);
  lcPPt.DrawLatex(90, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[1], pT_f[1]));
  lcPPt.SetTextColor(colpt[3]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.6, 1);
  lcPPt.DrawLatex(90, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[2], pT_f[2]));

  c->SaveAs(Form("plots/MCdists/PPt_scale.pdf"));
  c->Clear();

  // plotting mu- pT
  for(int i = 0; i < 6; i++) {

    hM_pT[i]->SetStats(0);
    hM_pT[i]->SetLineColor(colpt[i/2+1]);
    hM_pT[i]->SetMinimum(1e2);
    hM_pT[i]->SetMaximum(5e4);
    hM_pT[i]->GetXaxis()->SetRangeUser(0,125);
    hM_pT[i]->GetXaxis()->SetTitle("#it{p}_{T}(#mu^{-}) (GeV)");
    hM_pT[i]->GetXaxis()->SetTitleOffset(1.2);
    hM_pT[i]->GetXaxis()->CenterTitle(true);
    hM_pT[i]->GetXaxis()->SetLabelOffset(0.01);
    hM_pT[i]->GetYaxis()->SetLabelOffset(0.01);
    hM_pT[i]->GetYaxis()->SetTitle("dN/d#it{p}_{T}");
    hM_pT[i]->SetTitle("");
    if(i == 0) {
      hM_pT[i]->Draw("histo");
    }
    else if(i%2 == 1) {
      hM_pT[i]->SetLineStyle(kDashed);
      hM_pT[i]->Scale(hM_pT[i-1]->Integral()/hM_pT[i]->Integral());
      hM_pT[i]->Draw("histo same");
    }
    else {
      hM_pT[i]->Draw("histo same");
    }
  }

  legcPPt->Draw();
  
  TLatex lcMPt;
  lcMPt.SetTextSize(0.04);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.925, 1);
  lcMPt.DrawLatex(45, yl, "#bf{2018 J/#psi}");
  lcMPt.SetTextColor(colpt[1]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.76, 1);
  lcMPt.DrawLatex(90, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[0], pT_f[0]));
  lcMPt.SetTextColor(colpt[2]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.68, 1);
  lcMPt.DrawLatex(90, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[1], pT_f[1]));
  lcMPt.SetTextColor(colpt[3]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.6, 1);
  lcMPt.DrawLatex(90, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[2], pT_f[2]));

  c->SaveAs(Form("plots/MCdists/MPt_scale.pdf"));
  c->Clear();

  // plotting mu+ eta
  for(int i = 0; i < 6; i++) {

    hP_eta[i]->SetStats(0);
    hP_eta[i]->SetLineColor(colpt[i/2+1]);
    hP_eta[i]->SetMinimum(1e3);
    hP_eta[i]->SetMaximum(5e4);
    hP_eta[i]->GetXaxis()->SetTitle("#eta(#mu^{+})");
    hP_eta[i]->GetXaxis()->SetTitleOffset(1.2);
    hP_eta[i]->GetXaxis()->CenterTitle(true);
    hP_eta[i]->GetXaxis()->SetLabelOffset(0.01);
    hP_eta[i]->GetYaxis()->SetLabelOffset(0.01);
    hP_eta[i]->GetYaxis()->SetTitle("dN/d#eta");
    hP_eta[i]->SetTitle("");
    if(i == 0) {
      hP_eta[i]->Draw("histo");
    }
    else if(i%2 == 1) {
      hP_eta[i]->SetLineStyle(kDashed);
      hP_eta[i]->Scale(hP_eta[i-1]->Integral()/hP_eta[i]->Integral());
      hP_eta[i]->Draw("histo same");
    }
    else {
      hP_eta[i]->Draw("histo same");
    }
  }

  TLegend *legcPEta = new TLegend(0.55, 0.35, 0.85, 0.475);
  legcPEta->SetTextSize(0.04);
  legcPEta->SetBorderSize(0);
  legcPEta->SetFillColorAlpha(kWhite,0);
  legcPEta->AddEntry(h_pT[0], "lower-#it{p}_{T} MC", "l");
  legcPEta->AddEntry(h_pT[1], "higher-#it{p}_{T} MC", "l");
  legcPEta->Draw();

  TLatex lcPEta;
  lcPEta.SetTextSize(0.04);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.91, 1);
  lcPEta.DrawLatex(-1.2, yl, "#bf{2018 J/#psi}");
  lcPEta.SetTextColor(colpt[1]);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.42, 1);
  lcPEta.DrawLatex(-1, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[0], pT_f[0]));
  lcPEta.SetTextColor(colpt[2]);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.34, 1);
  lcPEta.DrawLatex(-1, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[1], pT_f[1]));
  lcPEta.SetTextColor(colpt[3]);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.26, 1);
  lcPEta.DrawLatex(-1, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[2], pT_f[2]));

  c->SaveAs(Form("plots/MCdists/PEta_scale.pdf"));
  c->Clear();

    // eta(mu+) ratio data/MC
  c->SetLogy(0);
  TH1D **hP_etar = new TH1D*[3]; 
  for(int i = 0; i < 3; i++) {
    hP_etar[i] = (TH1D*)hP_eta[2*i]->Clone(Form("rHp_%d", i));
    hP_etar[i]->Sumw2();
    hP_etar[i]->Divide(hP_eta[2*i+1]);
    // scaling up and down
    for(int j = 0; j < hP_etar[i]->GetNbinsX(); j++) {
      hP_etar[i]->SetBinContent(j+1, hP_etar[i]->GetBinContent(j+1)+sc_fc[i]);
    }
    
    hP_etar[i]->SetStats(0);
    hP_etar[i]->SetLineColor(colpt[i+1]);
    hP_etar[i]->SetMinimum(0.51);
    hP_etar[i]->SetMaximum(1.49);
    hP_etar[i]->GetXaxis()->SetTitle("#eta(#mu^{+})");
    hP_etar[i]->GetXaxis()->SetTitleOffset(1.2);
    hP_etar[i]->GetXaxis()->CenterTitle(true);
    hP_etar[i]->GetXaxis()->SetLabelOffset(0.01);
    hP_etar[i]->GetYaxis()->SetLabelOffset(0.01);
    hP_etar[i]->GetYaxis()->SetTitle("MC ratio");
    hP_etar[i]->SetTitle("");
    if(i == 0) {
      hP_etar[i]->Draw("histo");
    }
    else {
      hP_etar[i]->Draw("histo same");
    }
  }
  
  TLatex lcPEtar;
  lcPEtar.SetTextSize(0.04);
  yl = getPos(hP_etar[0]->GetMinimum(), hP_etar[0]->GetMaximum(), 0.91, 0);
  lcPEtar.DrawLatex(-1.1, yl, "#bf{2018 J/#psi}");
  lcPEtar.SetTextColor(colpt[1]);
  yl = getPos(hP_etar[0]->GetMinimum(), hP_etar[0]->GetMaximum(), 0.92, 0);
  lcPEtar.DrawLatex(0.45, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[0], pT_f[0]));
  lcPEtar.SetTextColor(colpt[2]);
  yl = getPos(hP_etar[0]->GetMinimum(), hP_etar[0]->GetMaximum(), 0.84, 0);
  lcPEtar.DrawLatex(0.45, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[1], pT_f[1]));
  lcPEtar.SetTextColor(colpt[3]);
  yl = getPos(hP_etar[0]->GetMinimum(), hP_etar[0]->GetMaximum(), 0.76, 0);
  lcPEtar.DrawLatex(0.45, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[2], pT_f[2]));

  for(int i = 0; i < 3; i++) {
    ly[i]->Draw();
  }
  
  c->SaveAs(Form("plots/MCdists/PEta_ratio.pdf"));
  c->Clear();

  // plotting mu- eta
  for(int i = 0; i < 6; i++) {

    hM_eta[i]->SetStats(0);
    hM_eta[i]->SetLineColor(colpt[i/2+1]);
    hM_eta[i]->SetMinimum(1e3);
    hM_eta[i]->SetMaximum(5e4);
    hM_eta[i]->GetXaxis()->SetTitle("#eta(#mu^{-})");
    hM_eta[i]->GetXaxis()->SetTitleOffset(1.2);
    hM_eta[i]->GetXaxis()->CenterTitle(true);
    hM_eta[i]->GetXaxis()->SetLabelOffset(0.01);
    hM_eta[i]->GetYaxis()->SetLabelOffset(0.01);
    hM_eta[i]->GetYaxis()->SetTitle("dN/d#eta");
    hM_eta[i]->SetTitle("");
    if(i == 0) {
      hM_eta[i]->Draw("histo");
    }
    else if(i%2 == 1) {
      hM_eta[i]->SetLineStyle(kDashed);
      hM_eta[i]->Scale(hM_eta[i-1]->Integral()/hM_eta[i]->Integral());
      hM_eta[i]->Draw("histo same");
    }
    else {
      hM_eta[i]->Draw("histo same");
    }
  }

  legcPEta->Draw();

  TLatex lcMEta;
  lcMEta.SetTextSize(0.04);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.91, 1);
  lcMEta.DrawLatex(-1.2, yl, "#bf{2018 J/#psi}");
  lcMEta.SetTextColor(colpt[1]);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.42, 1);
  lcMEta.DrawLatex(-1, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[0], pT_f[0]));
  lcMEta.SetTextColor(colpt[2]);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.34, 1);
  lcMEta.DrawLatex(-1, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[1], pT_f[1]));
  lcMEta.SetTextColor(colpt[3]);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.26, 1);
  lcMEta.DrawLatex(-1, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[2], pT_f[2]));

  c->SaveAs(Form("plots/MCdists/MEta_scale.pdf"));
  c->Clear();

    // eta(mu-) ratio data/MC
  c->SetLogy(0);
  TH1D **hM_etar = new TH1D*[3]; 
  for(int i = 0; i < 3; i++) {
    hM_etar[i] = (TH1D*)hM_eta[2*i]->Clone(Form("rHm_%d", i));
    hM_etar[i]->Sumw2();
    hM_etar[i]->Divide(hM_eta[2*i+1]);
    // scaling up and down
    for(int j = 0; j < hM_etar[i]->GetNbinsX(); j++) {
      hM_etar[i]->SetBinContent(j+1, hM_etar[i]->GetBinContent(j+1)+sc_fc[i]);
    }
    
    hM_etar[i]->SetStats(0);
    hM_etar[i]->SetLineColor(colpt[i+1]);
    hM_etar[i]->SetMinimum(0.51);
    hM_etar[i]->SetMaximum(1.49);
    hM_etar[i]->GetXaxis()->SetTitle("#eta(#mu^{-})");
    hM_etar[i]->GetXaxis()->SetTitleOffset(1.2);
    hM_etar[i]->GetXaxis()->CenterTitle(true);
    hM_etar[i]->GetXaxis()->SetLabelOffset(0.01);
    hM_etar[i]->GetYaxis()->SetLabelOffset(0.01);
    hM_etar[i]->GetYaxis()->SetTitle("MC ratio");
    hM_etar[i]->SetTitle("");
    if(i == 0) {
      hM_etar[i]->Draw("histo");
    }
    else {
      hM_etar[i]->Draw("histo same");
    }
  }
  
  TLatex lcMEtar;
  lcMEtar.SetTextSize(0.04);
  yl = getPos(hM_etar[0]->GetMinimum(), hM_etar[0]->GetMaximum(), 0.91, 0);
  lcMEtar.DrawLatex(-1.1, yl, "#bf{2018 J/#psi}");
  lcMEtar.SetTextColor(colpt[1]);
  yl = getPos(hM_etar[0]->GetMinimum(), hM_etar[0]->GetMaximum(), 0.92, 0);
  lcMEtar.DrawLatex(0.45, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[0], pT_f[0]));
  lcMEtar.SetTextColor(colpt[2]);
  yl = getPos(hM_etar[0]->GetMinimum(), hM_etar[0]->GetMaximum(), 0.84, 0);
  lcMEtar.DrawLatex(0.45, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[1], pT_f[1]));
  lcMEtar.SetTextColor(colpt[3]);
  yl = getPos(hM_etar[0]->GetMinimum(), hM_etar[0]->GetMaximum(), 0.76, 0);
  lcMEtar.DrawLatex(0.45, yl, Form("#bf{[%.0f,%.0f] GeV}", pT_i[2], pT_f[2]));

  for(int i = 0; i < 3; i++) {
    ly[i]->Draw();
  }
  
  c->SaveAs(Form("plots/MCdists/MEta_ratio.pdf"));
  c->Clear();

  c->Destructor();
  fin->Close();
}
