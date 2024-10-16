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

  // 2 pT dists: PRSR data + the MC (4 pT regions)
  TH1D **hP_pT = new TH1D*[2]; 
  TH1D **hM_pT = new TH1D*[2]; 
  // 2 eta dists: PRSR data + the MC (4 pT regions)
  TH1D **hP_eta = new TH1D*[2]; 
  TH1D **hM_eta = new TH1D*[2]; 
  
  TFile *fin = new TFile("files/store_muDists.root");
 
  // store the dists
  string lbl[] = {"Data", "MC"};
  for(int i = 0; i < 2; i++) {
    hP_pT[i] = (TH1D*)fin->Get(Form("hP_pT_%s", lbl[i].c_str()));
    hM_pT[i] = (TH1D*)fin->Get(Form("hM_pT_%s", lbl[i].c_str()));
    hP_eta[i] = (TH1D*)fin->Get(Form("hP_eta_%s", lbl[i].c_str()));
    hM_eta[i] = (TH1D*)fin->Get(Form("hM_eta_%s", lbl[i].c_str()));
  }

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  
  // plot pT (mu+) - scaling MC to data
  c->SetLogy();
  c->SetLeftMargin(0.11);
  int col_v[] = {kBlack};

  for(int i = 0; i < 2; i++) {

    hP_pT[i]->SetStats(0);
    hP_pT[i]->SetLineColor(col_v[0]);
    hP_pT[i]->SetMinimum(1e0);
    hP_pT[i]->SetMaximum(2e6);
    hP_pT[i]->GetXaxis()->SetTitle("#it{p}_{T} (#mu^{+}) (GeV)");
    hP_pT[i]->GetXaxis()->SetTitleOffset(1.1);
    hP_pT[i]->GetXaxis()->CenterTitle(true);
    hP_pT[i]->GetXaxis()->SetLabelOffset(0.01);
    hP_pT[i]->GetYaxis()->SetLabelOffset(0.01);
    hP_pT[i]->GetYaxis()->SetTitle("dN/d#it{p}_{T}");
    hP_pT[i]->SetTitle("");

    if(i > 0) {
      //      hP_pT[i]->Scale(hP_pT[i-1]->Integral() / hP_pT[i]->Integral());
      hP_pT[i]->SetLineStyle(kDashed);
    }
    if(i == 0) hP_pT[i]->Draw("histo");
    else hP_pT[i]->Draw("histo same");  

  }

  TLegend *legPPts = new TLegend(0.65, 0.825, 0.95, 0.95);
  legPPts->SetTextSize(0.04);
  legPPts->SetBorderSize(0);
  legPPts->SetFillColorAlpha(kWhite,0);
  legPPts->AddEntry(hP_pT[0], "PRS Data", "l");
  legPPts->AddEntry(hP_pT[1], "MC", "l");
  legPPts->Draw();

  TLatex lcPPts;
  lcPPts.SetTextSize(0.04);
  double yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.9, 1);
  lcPPts.DrawLatex(25, yl, "#bf{2017 #psi(2S)}");
  yl = getPos(hP_pT[0]->GetMinimum(), hP_pT[0]->GetMaximum(), 0.7, 1);
  lcPPts.DrawLatex(70, yl, "#bf{[20,100] GeV}");

  c->SaveAs(Form("plots/muDists/muPPt_all.pdf"));
  c->Clear();

  // plot pT (mu-) - scaling MC to data
  for(int i = 0; i < 2; i++) {

    hM_pT[i]->SetStats(0);
    hM_pT[i]->SetLineColor(col_v[0]);
    hM_pT[i]->SetMinimum(1e0);
    hM_pT[i]->SetMaximum(2e6);
    hM_pT[i]->GetXaxis()->SetTitle("#it{p}_{T} (#mu^{-}) (GeV)");
    hM_pT[i]->GetXaxis()->SetTitleOffset(1.1);
    hM_pT[i]->GetXaxis()->CenterTitle(true);
    hM_pT[i]->GetXaxis()->SetLabelOffset(0.01);
    hM_pT[i]->GetYaxis()->SetLabelOffset(0.01);
    hM_pT[i]->GetYaxis()->SetTitle("dN/d#it{p}_{T}");
    hM_pT[i]->SetTitle("");

    if(i > 0) {
      //      hM_pT[i]->Scale(hM_pT[i-1]->Integral() / hM_pT[i]->Integral());
      hM_pT[i]->SetLineStyle(kDashed);
    }
    if(i == 0) hM_pT[i]->Draw("histo");
    else hM_pT[i]->Draw("histo same");  

  }

  TLegend *legMPts = new TLegend(0.65, 0.825, 0.95, 0.95);
  legMPts->SetTextSize(0.04);
  legMPts->SetBorderSize(0);
  legMPts->SetFillColorAlpha(kWhite,0);
  legMPts->AddEntry(hM_pT[0], "PRS Data", "l");
  legMPts->AddEntry(hM_pT[1], "MC", "l");
  legMPts->Draw();

  TLatex lcMPts;
  lcMPts.SetTextSize(0.04);
  yl = getPos(hM_pT[0]->GetMinimum(), hM_pT[0]->GetMaximum(), 0.9, 1);
  lcMPts.DrawLatex(25, yl, "#bf{2017 #psi(2S)}");
  yl = getPos(hM_pT[0]->GetMinimum(), hM_pT[0]->GetMaximum(), 0.7, 1);
  lcMPts.DrawLatex(70, yl, "#bf{[20,100] GeV}");
  
  c->SaveAs(Form("plots/muDists/muMPt_all.pdf"));
  c->Clear();

  // now the etas
  // plot eta(mu+) - scaling MC to data
  for(int i = 0; i < 2; i++) {

    hP_eta[i]->SetStats(0);
    hP_eta[i]->SetLineColor(col_v[0]);
    hP_eta[i]->SetMinimum(8e2);
    hP_eta[i]->SetMaximum(9e4);
    hP_eta[i]->GetXaxis()->SetTitle("#eta (#mu^{+})");
    hP_eta[i]->GetXaxis()->SetTitleOffset(1.1);
    hP_eta[i]->GetXaxis()->CenterTitle(true);
    hP_eta[i]->GetXaxis()->SetLabelOffset(0.01);
    hP_eta[i]->GetYaxis()->SetLabelOffset(0.01);
    hP_eta[i]->GetYaxis()->SetTitle("dN/d#eta");
    hP_eta[i]->SetTitle("");

    if(i > 0) {
      hP_eta[i]->Scale(hP_eta[i-1]->Integral() / hP_eta[i]->Integral());
      hP_eta[i]->SetLineStyle(kDashed);
    }
    if(i == 0) hP_eta[i]->Draw("histo");
    else hP_eta[i]->Draw("histo same");  

  }

  TLegend *legPEtas = new TLegend(0.6, 0.275, 0.9, 0.4);
  legPEtas->SetTextSize(0.04);
  legPEtas->SetBorderSize(0);
  legPEtas->SetFillColorAlpha(kWhite,0);
  legPEtas->AddEntry(hP_eta[0], "PRS Data", "l");
  legPEtas->AddEntry(hP_eta[1], "Scaled MC", "l");
  legPEtas->Draw();

  TLatex lcPEtas;
  lcPEtas.SetTextSize(0.04);
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.9, 1);
  lcPEtas.DrawLatex(-1.25, yl, "#bf{2017 #psi(2S)}");
  yl = getPos(hP_eta[0]->GetMinimum(), hP_eta[0]->GetMaximum(), 0.28, 1);
  lcPEtas.DrawLatex(-1.1, yl, "#bf{[20,100] GeV}");

  c->SaveAs(Form("plots/muDists/muPEta_scale.pdf"));
  c->Clear();

  // eta ratio data/MC
  c->SetLogy(0);
  TH1D **hP_etar = new TH1D*[1]; 
  for(int i = 0; i < 1; i++) {
    hP_etar[i] = (TH1D*)hP_eta[i]->Clone(Form("rHp_%d",i));
    hP_etar[i]->Sumw2();
    hP_etar[i]->Divide(hP_eta[i+1]);

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
  lcPEtar.DrawLatex(-1.22, yl, "#bf{2017 #psi(2S)}");
  yl = getPos(hP_etar[0]->GetMinimum(), hP_etar[0]->GetMaximum(), 0.9, 0);
  lcPEtar.DrawLatex(0.3, yl, "#bf{[20,100] GeV}");
  
  c->SaveAs(Form("plots/muDists/muPEta_ratio.pdf"));
  c->Clear();

  // plot eta(mu-) - scaling MC to data
  for(int i = 0; i < 2; i++) {

    hM_eta[i]->SetStats(0);
    hM_eta[i]->SetLineColor(col_v[0]);
    hM_eta[i]->SetMinimum(8e2);
    hM_eta[i]->SetMaximum(9e4);
    hM_eta[i]->GetXaxis()->SetTitle("#eta (#mu^{-})");
    hM_eta[i]->GetXaxis()->SetTitleOffset(1.1);
    hM_eta[i]->GetXaxis()->CenterTitle(true);
    hM_eta[i]->GetXaxis()->SetLabelOffset(0.01);
    hM_eta[i]->GetYaxis()->SetLabelOffset(0.01);
    hM_eta[i]->GetYaxis()->SetTitle("dN/d#eta");
    hM_eta[i]->SetTitle("");

    if(i > 0) {
      hM_eta[i]->Scale(hM_eta[i-1]->Integral() / hM_eta[i]->Integral());
      hM_eta[i]->SetLineStyle(kDashed);
    }
    if(i == 0) hM_eta[i]->Draw("histo");
    else hM_eta[i]->Draw("histo same");  

  }
  
  TLegend *legMEtas = new TLegend(0.6, 0.275, 0.9, 0.4);
  legMEtas->SetTextSize(0.04);
  legMEtas->SetBorderSize(0);
  legMEtas->SetFillColorAlpha(kWhite,0);
  legMEtas->AddEntry(hP_pT[0], "PRS Data", "l");
  legMEtas->AddEntry(hP_pT[1], "Scaled MC", "l");
  legMEtas->Draw();

  TLatex lcMEtas;
  lcMEtas.SetTextSize(0.04);
  yl = getPos(hM_eta[0]->GetMinimum(), hM_eta[0]->GetMaximum(), 0.9, 1);
  lcMEtas.DrawLatex(-1.25, yl, "#bf{2017 #psi(2S)}");
  yl = getPos(hM_eta[0]->GetMinimum(), hM_eta[0]->GetMaximum(), 0.28, 1);
  lcMEtas.DrawLatex(-1.1, yl, "#bf{[20,100] GeV}");
  
  c->SaveAs(Form("plots/muDists/muMEta_scale.pdf"));
  c->Clear();

  // eta ratio data/MC
  c->SetLogy(0);
  TH1D **hM_etar = new TH1D*[1]; 
  for(int i = 0; i < 1; i++) {
    hM_etar[i] = (TH1D*)hM_eta[i]->Clone(Form("rHm_%d",i));
    hM_etar[i]->Sumw2();
    hM_etar[i]->Divide(hM_eta[i+1]);

    hM_etar[i]->SetStats(0);
    hM_etar[i]->SetLineColor(col_v[i]);
    hM_etar[i]->SetMinimum(0.51);
    hM_etar[i]->SetMaximum(1.49);
    hM_etar[i]->GetXaxis()->SetTitle("#eta (#mu^{+})");
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
  lcMEtar.DrawLatex(-1.22, yl, "#bf{2017 #psi(2S)}");
  yl = getPos(hM_etar[0]->GetMinimum(), hM_etar[0]->GetMaximum(), 0.9, 0);
  lcMEtar.DrawLatex(0.3, yl, "#bf{[20,100] GeV}");
  
  c->SaveAs(Form("plots/muDists/muMEta_ratio.pdf"));
  c->Clear();

  c->Destructor();
  fin->Close();
}
