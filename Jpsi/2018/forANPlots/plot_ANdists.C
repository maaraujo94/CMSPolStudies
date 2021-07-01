// macro to take the stored histos and plot nicely

int col_data(int inp) {
  if(inp < 3) return kBlack;
  else return inp-1;
}

int col_three(int inp) {
  if(inp < 4) return kBlack;
  else if(inp < 8) return kViolet;
  else return inp-7;
}

void plot_ANdists()
{
  // PART 1 : reading the histograms

  // four pT dists: PRSR data + the 3 MC
  TH1D **h_pT = new TH1D*[4]; 
  // 6 y dists: PRSR data + MC over 3 pT regions
  TH1D **h_y = new TH1D*[6]; 
  // 6 M dists: PR data + MC over 3 pT regions
  TH1D **h_m = new TH1D*[6]; 
  // 4 lifetime dists: data in 4 pT regions
  TH1D **h_lt = new TH1D*[4];
  // 12 costh dists: PR data + NP data + MC over 4 pT regions
  TH1D **h_cos = new TH1D*[12]; 
  
  TFile *fin = new TFile("store_ANdists.root");
 
  // read the pT dists
  string lbl_pt[] = {"Data", "lowPtMC", "midPtMC", "highPtMC"};
  for(int i = 0; i < 4; i++) {
    h_pT[i] = (TH1D*)fin->Get(Form("h_pT_%s", lbl_pt[i].c_str()));
  }

  // store the y dists
  string lbl_y[] = {"lowPtData", "midPtData", "highPtData", "lowPtMC", "midPtMC", "highPtMC"};
  for(int i = 0; i < 6; i++) {
    h_y[i] = (TH1D*)fin->Get(Form("h_y_%s", lbl_y[i].c_str()));
  }

  // store the M dists
  string lbl_m[] = {"lowPtData", "midPtData", "highPtData", "lowPtMC", "midPtMC", "highPtMC"};
  for(int i = 0; i < 6; i++) {
    h_m[i] = (TH1D*)fin->Get(Form("h_m_%s", lbl_m[i].c_str()));
  }

  // store the lifetime dists
  string lbl_lt[] = {"full", "lowPt", "midPt", "highPt"};
  for(int i = 0; i < 4; i++) {
    h_lt[i] = (TH1D*)fin->Get(Form("h_lt_%s", lbl_lt[i].c_str()));
  }

  // store the cos dists
  string lbl_cos[] = {"fullPR", "lowPtPR", "midPtPR", "highPtPR", "fullNP", "lowPtNP", "midPtNP", "highPtNP", "fullMC", "lowPtMC", "midPtMC", "highPtMC"};
  for(int i = 0; i < 12; i++) {
    h_cos[i] = (TH1D*)fin->Get(Form("h_cos_%s", lbl_cos[i].c_str()));
  }

  TCanvas *c = new TCanvas("", "", 900, 900);

  // plot pT
  c->SetLogy();
  c->SetLeftMargin(0.11);

  for(int i = 0; i < 4; i++) {

    h_pT[i]->SetStats(0);
    h_pT[i]->SetLineColor(i+1);
    h_pT[i]->SetMinimum(1e2);
    h_pT[i]->SetMaximum(1e7);
    h_pT[i]->GetXaxis()->SetTitle("p_{T}(#mu#mu) (GeV)");
    h_pT[i]->GetXaxis()->SetTitleOffset(1.1);
    h_pT[i]->GetYaxis()->SetTitle("dN/dp_{T}");
    h_pT[i]->SetTitle(Form("p_{T} (%s)", lbl_pt[i].c_str()));
    h_pT[i]->Draw("histo");
    
    c->SaveAs(Form("plots/ANdists/pt_%s.pdf", lbl_pt[i].c_str()));
    c->Clear();
  
  }

  // plot y
  c->SetLogy();

  for(int i = 0; i < 6; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(col_data(i));
    h_y[i]->SetMinimum(1e3);
    h_y[i]->SetMaximum(6e5);
    h_y[i]->GetXaxis()->SetTitle("y(#mu#mu)");
    h_y[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y[i]->GetYaxis()->SetTitle("dN/dy");
    h_y[i]->SetTitle(Form("y (%s)", lbl_y[i].c_str()));
    h_y[i]->Draw("histo");
    
    c->SaveAs(Form("plots/ANdists/y_%s.pdf", lbl_y[i].c_str()));
    c->Clear();
  
  }

  // plot lifetime
  c->SetLogy();

  for(int i = 0; i < 4; i++) {

    h_lt[i]->SetStats(0);
    h_lt[i]->SetLineColor(i+1);
    h_lt[i]->SetMinimum(1e3);
    h_lt[i]->SetMaximum(3e6);
    h_lt[i]->GetXaxis()->SetTitle("c#tau (#mum)");
    h_lt[i]->GetXaxis()->SetTitleOffset(1.1);
    h_lt[i]->GetYaxis()->SetTitle("dN/dc#tau");
    h_lt[i]->SetTitle(Form("lifetime (%s)", lbl_lt[i].c_str()));
    h_lt[i]->Draw("histo");

    TLine *l_PRm = new TLine(-100, h_lt[i]->GetMinimum(), -100, h_lt[i]->GetMaximum());
    l_PRm->SetLineStyle(kDashed);
    l_PRm->SetLineColor(kBlack);
    l_PRm->Draw();
    TLine *l_PRp = new TLine(100, h_lt[i]->GetMinimum(), 100, h_lt[i]->GetMaximum());
    l_PRp->SetLineStyle(kDashed);
    l_PRp->SetLineColor(kBlack);
    l_PRp->Draw();

    TLine *l_NPm = new TLine(140, h_lt[i]->GetMinimum(), 140, h_lt[i]->GetMaximum());
    l_NPm->SetLineStyle(kDashed);
    l_NPm->SetLineColor(kBlack);
    l_NPm->Draw();
    TLine *l_NPp = new TLine(500, h_lt[i]->GetMinimum(), 500, h_lt[i]->GetMaximum());
    l_NPp->SetLineStyle(kDashed);
    l_NPp->SetLineColor(kBlack);
    l_NPp->Draw();

    c->SaveAs(Form("plots/ANdists/lt_%s.pdf", lbl_lt[i].c_str()));
    c->Clear();
  
  }

  // plot M
  c->SetLogy(0);
  c->SetLeftMargin(0.12);

  for(int i = 0; i < 6; i++) {

    h_m[i]->SetStats(0);
    h_m[i]->Scale(1./h_m[i]->Integral());
    h_m[i]->SetLineColor(col_data(i));
    h_m[i]->SetMinimum(0);
    h_m[i]->SetMaximum(0.085);
    h_m[i]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
    h_m[i]->GetXaxis()->SetTitleOffset(1.1);
    h_m[i]->GetYaxis()->SetTitle("dN/dM (normalized)");
    h_m[i]->SetTitle(Form("M (%s)", lbl_y[i].c_str()));
    h_m[i]->Draw("histo");
    
    c->SaveAs(Form("plots/ANdists/m_%s.pdf", lbl_y[i].c_str()));
    c->Clear();
  
  }

  // plot costh
  c->SetLogy(0);

  for(int i = 0; i < 12; i++) {

    h_cos[i]->SetStats(0);
    h_cos[i]->Scale(1./h_cos[i]->Integral());
    h_cos[i]->SetLineColor(col_three(i));
    h_cos[i]->SetMinimum(0);
    h_cos[i]->SetMaximum(0.09);
    h_cos[i]->GetXaxis()->SetTitle("|cos#theta|");
    h_cos[i]->GetXaxis()->SetTitleOffset(1.);
    h_cos[i]->GetYaxis()->SetTitle("dN/d|cos#theta| (normalized)");
    h_cos[i]->SetTitle(Form("cos#theta (%s)", lbl_cos[i].c_str()));
    h_cos[i]->Draw("histo");
    
    c->SaveAs(Form("plots/ANdists/cos_%s.pdf", lbl_cos[i].c_str()));
    c->Clear();
  
  }

  fin->Close();
}
