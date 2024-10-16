// macro to take the stored histos and plot nicely (single muon)

int col_data(int inp) {
  if(inp < 3) return kBlack;
  else return inp-1;
}

int col_three(int inp) {
  if(inp < 4) return kBlack;
  else if(inp < 8) return kViolet;
  else return inp-7;
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


void plot_muDists()
{
  // PART 1 : reading the histograms

  // 8 pT dists: PRSR data + the MC (4 pT regions)
  TH1D **hP_pT = new TH1D*[8]; 
  TH1D **hM_pT = new TH1D*[8]; 
  // 8 eta dists: PRSR data + the MC (4 pT regions)
  TH1D **hP_eta = new TH1D*[8]; 
  TH1D **hM_eta = new TH1D*[8]; 
  
  TFile *fin = new TFile("files/store_muDists.root");
 
  // store the dists
  string lbl[] = {"lowPtData", "midPtData", "highPtData", "highestPtData", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC"};
  for(int i = 0; i < 8; i++) {
    hP_pT[i] = (TH1D*)fin->Get(Form("hP_pT_%s", lbl[i].c_str()));
    hM_pT[i] = (TH1D*)fin->Get(Form("hM_pT_%s", lbl[i].c_str()));
    hP_eta[i] = (TH1D*)fin->Get(Form("hP_eta_%s", lbl[i].c_str()));
    hM_eta[i] = (TH1D*)fin->Get(Form("hM_eta_%s", lbl[i].c_str()));
  }

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  
  // plot pT (mu+) - don't scale
  c->SetLogy();
  c->SetLeftMargin(0.11);
  int col_v[] = {kRed+1, kViolet+2, kGreen+3, kBlue};

  for(int i = 0; i < 8; i++) {

    hP_pT[i]->SetStats(0);
    hP_pT[i]->SetLineColor(col_v[i%4]);
    hP_pT[i]->SetMinimum(1e1);
    hP_pT[i]->SetMaximum(4e6);
    hP_pT[i]->GetXaxis()->SetTitle("#it{p}_{T} (#mu^{+}) (GeV)");
    hP_pT[i]->GetXaxis()->SetTitleOffset(1.1);
    hP_pT[i]->GetXaxis()->CenterTitle(true);
    hP_pT[i]->GetXaxis()->SetLabelOffset(0.01);
    hP_pT[i]->GetYaxis()->SetLabelOffset(0.01);
    hP_pT[i]->GetYaxis()->SetTitle("dN/d#it{p}_{T}");
    hP_pT[i]->SetTitle("");

    if(i > 3) {
      //    hP_pT[i]->Scale(hP_pT[i-4]->Integral() / hP_pT[i]->Integral());
      hP_pT[i]->SetLineStyle(kDashed);
    }
    if(i == 0) hP_pT[i]->Draw("histo");
    else hP_pT[i]->Draw("histo same");  

  }

  // just formatting another histo for the y legend
  hP_eta[1]->SetLineColor(kBlack);
  hP_eta[1]->SetLineStyle(kSolid);
  hP_eta[5]->SetLineColor(kBlack);
  hP_eta[5]->SetLineStyle(kDashed);
  
  TLegend *legPPts = new TLegend(0.65, 0.85, 0.95, 0.975);
  legPPts->SetTextSize(0.04);
  legPPts->SetBorderSize(0);
  legPPts->SetFillColorAlpha(kWhite,0);
  legPPts->AddEntry(hP_eta[1], "PRS Data", "l");
  legPPts->AddEntry(hP_eta[5], "MC", "l");
  legPPts->Draw();

  TLatex lcPPts;
  lcPPts.SetTextSize(0.04);
  double yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.9, 1);
  lcPPts.DrawLatex(35, yl, "#bf{2018 J/#psi}");
  lcPPts.SetTextColor(col_v[0]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.76, 1);
  lcPPts.DrawLatex(95, yl, "#bf{[25,45] GeV}");
  lcPPts.SetTextColor(col_v[1]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.68, 1);
  lcPPts.DrawLatex(95, yl, "#bf{[45,50] GeV}");
  lcPPts.SetTextColor(col_v[2]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.6, 1);
  lcPPts.DrawLatex(95, yl, "#bf{[50,70] GeV}");
  lcPPts.SetTextColor(col_v[3]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.52, 1);
  lcPPts.DrawLatex(95, yl, "#bf{[70,120] GeV}");
  
  c->SaveAs(Form("plots/muDists/muPPt_all.pdf"));
  c->Clear();

  // plot pT (mu-) - don't scale
  for(int i = 0; i < 8; i++) {

    hM_pT[i]->SetStats(0);
    hM_pT[i]->SetLineColor(col_v[i%4]);
    hM_pT[i]->SetMinimum(1e1);
    hM_pT[i]->SetMaximum(4e6);
    hM_pT[i]->GetXaxis()->SetTitle("#it{p}_{T} (#mu^{-}) (GeV)");
    hM_pT[i]->GetXaxis()->SetTitleOffset(1.1);
    hM_pT[i]->GetXaxis()->CenterTitle(true);
    hM_pT[i]->GetXaxis()->SetLabelOffset(0.01);
    hM_pT[i]->GetYaxis()->SetLabelOffset(0.01);
    hM_pT[i]->GetYaxis()->SetTitle("dN/d#it{p}_{T}");
    hM_pT[i]->SetTitle("");

    if(i > 3) {
      //    hM_pT[i]->Scale(hM_pT[i-4]->Integral() / hM_pT[i]->Integral());
      hM_pT[i]->SetLineStyle(kDashed);
    }
    if(i == 0) hM_pT[i]->Draw("histo");
    else hM_pT[i]->Draw("histo same");  

  }

  legPPts->Draw();

  TLatex lcMPts;
  lcMPts.SetTextSize(0.04);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.9, 1);
  lcMPts.DrawLatex(35, yl, "#bf{2018 J/#psi}");
  lcMPts.SetTextColor(col_v[0]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.76, 1);
  lcMPts.DrawLatex(95, yl, "#bf{[25,45] GeV}");
  lcMPts.SetTextColor(col_v[1]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.68, 1);
  lcMPts.DrawLatex(95, yl, "#bf{[45,50] GeV}");
  lcMPts.SetTextColor(col_v[2]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.6, 1);
  lcMPts.DrawLatex(95, yl, "#bf{[50,70] GeV}");
  lcMPts.SetTextColor(col_v[3]);
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.52, 1);
  lcMPts.DrawLatex(95, yl, "#bf{[70,120] GeV}");
  
  c->SaveAs(Form("plots/muDists/muMPt_all.pdf"));
  c->Clear();

  // now the etas
  // plot eta(mu+) - scaling MC to data
  for(int i = 0; i < 8; i++) {

    hP_eta[i]->SetStats(0);
    hP_eta[i]->SetLineColor(col_v[i%4]);
    hP_eta[i]->SetMinimum(1e1);
    hP_eta[i]->SetMaximum(9e5);
    hP_eta[i]->GetXaxis()->SetTitle("#eta (#mu^{+})");
    hP_eta[i]->GetXaxis()->SetTitleOffset(1.1);
    hP_eta[i]->GetXaxis()->CenterTitle(true);
    hP_eta[i]->GetXaxis()->SetLabelOffset(0.01);
    hP_eta[i]->GetYaxis()->SetLabelOffset(0.01);
    hP_eta[i]->GetYaxis()->SetTitle("dN/d#eta");
    hP_eta[i]->SetTitle("");

    if(i > 3) {
      hP_eta[i]->Scale(hP_eta[i-4]->Integral() / hP_eta[i]->Integral());
      hP_eta[i]->SetLineStyle(kDashed);
    }
    if(i == 0) hP_eta[i]->Draw("histo");
    else hP_eta[i]->Draw("histo same");  

  }

  // just formatting another histo for the y legend
  hP_pT[1]->SetLineColor(kBlack);
  hP_pT[1]->SetLineStyle(kSolid);
  hP_pT[5]->SetLineColor(kBlack);
  hP_pT[5]->SetLineStyle(kDashed);
  
  TLegend *legPEtas = new TLegend(0.6, 0.25, 0.9, 0.375);
  legPEtas->SetTextSize(0.04);
  legPEtas->SetBorderSize(0);
  legPEtas->SetFillColorAlpha(kWhite,0);
  legPEtas->AddEntry(hP_pT[1], "PRS Data", "l");
  legPEtas->AddEntry(hP_pT[5], "Scaled MC", "l");
  legPEtas->Draw();

  TLatex lcPEtas;
  lcPEtas.SetTextSize(0.04);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.9, 1);
  lcPEtas.DrawLatex(-1.25, yl, "#bf{2018 J/#psi}");
  lcPEtas.SetTextColor(col_v[0]);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.32, 1);
  lcPEtas.DrawLatex(-1.1, yl, "#bf{[25,45] GeV}");
  lcPEtas.SetTextColor(col_v[1]);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.24, 1);
  lcPEtas.DrawLatex(-1.1, yl, "#bf{[45,50] GeV}");
  lcPEtas.SetTextColor(col_v[2]);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.16, 1);
  lcPEtas.DrawLatex(-1.1, yl, "#bf{[50,70] GeV}");
  lcPEtas.SetTextColor(col_v[3]);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.08, 1);
  lcPEtas.DrawLatex(-1.1, yl, "#bf{[70,120] GeV}");
  
  c->SaveAs(Form("plots/muDists/muPEta_scale.pdf"));
  c->Clear();

    // eta ratio data/MC
  c->SetLogy(0);
  TH1D **hP_etar = new TH1D*[4]; 
  for(int i = 0; i < 4; i++) {
    hP_etar[i] = (TH1D*)hP_eta[i]->Clone(Form("rHp_%d",i));
    hP_etar[i]->Sumw2();
    hP_etar[i]->Divide(hP_eta[i+4]);

    hP_etar[i]->SetStats(0);
    hP_etar[i]->SetLineColor(col_v[i]);
    hP_etar[i]->SetMinimum(0.51);
    hP_etar[i]->SetMaximum(1.49);
    hP_etar[i]->GetXaxis()->SetTitle("#eta (#mu^{+})");
    hP_etar[i]->GetXaxis()->SetTitleOffset(1.1);
    hP_etar[i]->GetXaxis()->CenterTitle(true);
    hP_etar[i]->GetXaxis()->SetLabelOffset(0.01);
    hP_etar[i]->GetYaxis()->SetLabelOffset(0.01);
    hP_etar[i]->GetYaxis()->SetTitle("Data/MC");
    hP_etar[i]->SetTitle("");

    if(i == 0) hP_etar[i]->Draw("histo");
    else hP_etar[i]->Draw("histo same");  

  }
  
  TLatex lcPEtar;
  lcPEtar.SetTextSize(0.04);
  yl = getPos(hP_etar[0]->GetMinimum(), hP_etar[0]->GetMaximum(), 0.9, 0);
  lcPEtar.DrawLatex(-1.22, yl, "#bf{2018 J/#psi}");
  lcPEtar.SetTextColor(col_v[0]);
  yl = getPos(hP_etar[0]->GetMinimum(), hP_etar[0]->GetMaximum(), 0.9, 0);
  lcPEtar.DrawLatex(0.3, yl, "#bf{[25,45] GeV}");
  lcPEtar.SetTextColor(col_v[1]);
  yl = getPos(hP_etar[0]->GetMinimum(), hP_etar[0]->GetMaximum(), 0.82, 0);
  lcPEtar.DrawLatex(0.3, yl, "#bf{[45,50] GeV}");
  lcPEtar.SetTextColor(col_v[2]);
  yl = getPos(hP_etar[0]->GetMinimum(), hP_etar[0]->GetMaximum(), 0.74, 0);
  lcPEtar.DrawLatex(0.3, yl, "#bf{[50,70] GeV}");
  lcPEtar.SetTextColor(col_v[3]);
  yl = getPos(hP_etar[0]->GetMinimum(), hP_etar[0]->GetMaximum(), 0.66, 0);
  lcPEtar.DrawLatex(0.3, yl, "#bf{[70,120] GeV}");
  
  c->SaveAs(Form("plots/muDists/muPEta_ratio.pdf"));
  c->Clear();

  // plot eta(mu-) - scaling MC to data
  for(int i = 0; i < 8; i++) {

    hM_eta[i]->SetStats(0);
    hM_eta[i]->SetLineColor(col_v[i%4]);
    hM_eta[i]->SetMinimum(1e1);
    hM_eta[i]->SetMaximum(9e5);
    hM_eta[i]->GetXaxis()->SetTitle("#eta (#mu^{-})");
    hM_eta[i]->GetXaxis()->SetTitleOffset(1.1);
    hM_eta[i]->GetXaxis()->CenterTitle(true);
    hM_eta[i]->GetXaxis()->SetLabelOffset(0.01);
    hM_eta[i]->GetYaxis()->SetLabelOffset(0.01);
    hM_eta[i]->GetYaxis()->SetTitle("dN/d#eta");
    hM_eta[i]->SetTitle("");

    if(i > 3) {
      hM_eta[i]->Scale(hM_eta[i-4]->Integral() / hM_eta[i]->Integral());
      hM_eta[i]->SetLineStyle(kDashed);
    }
    if(i == 0) hM_eta[i]->Draw("histo");
    else hM_eta[i]->Draw("histo same");  

  }
  
  TLegend *legMEtas = new TLegend(0.6, 0.25, 0.9, 0.375);
  legMEtas->SetTextSize(0.04);
  legMEtas->SetBorderSize(0);
  legMEtas->SetFillColorAlpha(kWhite,0);
  legMEtas->AddEntry(hP_pT[1], "PRS Data", "l");
  legMEtas->AddEntry(hP_pT[5], "Scaled MC", "l");
  legMEtas->Draw();

  TLatex lcMEtas;
  lcMEtas.SetTextSize(0.04);
  yl = getPos(hM_eta[0]->GetMinimum(), hM_eta[0]->GetMaximum(), 0.9, 1);
  lcMEtas.DrawLatex(-1.25, yl, "#bf{2018 J/#psi}");
  lcMEtas.SetTextColor(col_v[0]);
  yl = getPos(hM_eta[0]->GetMinimum(), hM_eta[0]->GetMaximum(), 0.32, 1);
  lcMEtas.DrawLatex(-1.1, yl, "#bf{[25,45] GeV}");
  lcMEtas.SetTextColor(col_v[1]);
  yl = getPos(hM_eta[0]->GetMinimum(), hM_eta[0]->GetMaximum(), 0.24, 1);
  lcMEtas.DrawLatex(-1.1, yl, "#bf{[45,50] GeV}");
  lcMEtas.SetTextColor(col_v[2]);
  yl = getPos(hM_eta[0]->GetMinimum(), hM_eta[0]->GetMaximum(), 0.16, 1);
  lcMEtas.DrawLatex(-1.1, yl, "#bf{[50,70] GeV}");
  lcMEtas.SetTextColor(col_v[3]);
  yl = getPos(hM_eta[0]->GetMinimum(), hM_eta[0]->GetMaximum(), 0.08, 1);
  lcMEtas.DrawLatex(-1.1, yl, "#bf{[70,120] GeV}");
  
  c->SaveAs(Form("plots/muDists/muMEta_scale.pdf"));
  c->Clear();

    // eta ratio data/MC
  c->SetLogy(0);
  TH1D **hM_etar = new TH1D*[4]; 
  for(int i = 0; i < 4; i++) {
    hM_etar[i] = (TH1D*)hM_eta[i]->Clone(Form("rHp_%d",i));
    hM_etar[i]->Sumw2();
    hM_etar[i]->Divide(hM_eta[i+4]);

    hM_etar[i]->SetStats(0);
    hM_etar[i]->SetLineColor(col_v[i]);
    hM_etar[i]->SetMinimum(0.51);
    hM_etar[i]->SetMaximum(1.49);
    hM_etar[i]->GetXaxis()->SetTitle("#eta (#mu^{-})");
    hM_etar[i]->GetXaxis()->SetTitleOffset(1.1);
    hM_etar[i]->GetXaxis()->CenterTitle(true);
    hM_etar[i]->GetXaxis()->SetLabelOffset(0.01);
    hM_etar[i]->GetYaxis()->SetLabelOffset(0.01);
    hM_etar[i]->GetYaxis()->SetTitle("Data/MC");
    hM_etar[i]->SetTitle("");

    if(i == 0) hM_etar[i]->Draw("histo");
    else hM_etar[i]->Draw("histo same");  

  }
  
  TLatex lcMEtar;
  lcMEtar.SetTextSize(0.04);
  yl = getPos(hM_etar[0]->GetMinimum(), hM_etar[0]->GetMaximum(), 0.9, 0);
  lcMEtar.DrawLatex(-1.22, yl, "#bf{2018 J/#psi}");
  lcMEtar.SetTextColor(col_v[0]);
  yl = getPos(hM_etar[0]->GetMinimum(), hM_etar[0]->GetMaximum(), 0.9, 0);
  lcMEtar.DrawLatex(0.3, yl, "#bf{[25,45] GeV}");
  lcMEtar.SetTextColor(col_v[1]);
  yl = getPos(hM_etar[0]->GetMinimum(), hM_etar[0]->GetMaximum(), 0.82, 0);
  lcMEtar.DrawLatex(0.3, yl, "#bf{[45,50] GeV}");
  lcMEtar.SetTextColor(col_v[2]);
  yl = getPos(hM_etar[0]->GetMinimum(), hM_etar[0]->GetMaximum(), 0.74, 0);
  lcMEtar.DrawLatex(0.3, yl, "#bf{[50,70] GeV}");
  lcMEtar.SetTextColor(col_v[3]);
  yl = getPos(hM_etar[0]->GetMinimum(), hM_etar[0]->GetMaximum(), 0.66, 0);
  lcMEtar.DrawLatex(0.3, yl, "#bf{[70,120] GeV}");
  
  c->SaveAs(Form("plots/muDists/muMEta_ratio.pdf"));
  c->Clear();

  c->Destructor();
  fin->Close();
}
