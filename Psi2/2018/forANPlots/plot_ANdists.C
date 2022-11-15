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
  TH1D **h_pT = new TH1D*[2]; 
  // 6 y dists: PRSR data + MC over 3 pT regions
  TH1D **h_y = new TH1D*[2]; 
  // 7 M dists: PR data + MC over 3 pT regions + full pT data
  TH1D **h_m = new TH1D*[2]; 
  // 4 lifetime dists: data in 4 pT regions
  TH1D **h_lt = new TH1D*[1];
  // 12 costh dists: PR data + NP data + MC over 4 pT regions
  TH1D **h_cos = new TH1D*[3]; 
  
  TFile *fin = new TFile("files/store_ANdists.root");
 
  // read the pT dists
  string lbl_pt[] = {"Data", "MC"};
  for(int i = 0; i < 2; i++) {
    h_pT[i] = (TH1D*)fin->Get(Form("h_pT_%s", lbl_pt[i].c_str()));
  }

  // store the y dists
  string lbl_y[] = {"Data", "MC"};
  for(int i = 0; i < 2; i++) {
    h_y[i] = (TH1D*)fin->Get(Form("h_y_%s", lbl_y[i].c_str()));
  }

  // store the M dists
  string lbl_m[] = {"Data", "MC"};
  for(int i = 0; i < 2; i++) {
    h_m[i] = (TH1D*)fin->Get(Form("h_m_%s", lbl_m[i].c_str()));
  }

  // store the lifetime dists
  string lbl_lt[] = {"full"};
  for(int i = 0; i < 1; i++) {
    h_lt[i] = (TH1D*)fin->Get(Form("h_lt_%s", lbl_lt[i].c_str()));
  }

  // store the cos dists
  string lbl_cos[] = {"PR", "NP", "MC"};
  for(int i = 0; i < 3; i++) {
    h_cos[i] = (TH1D*)fin->Get(Form("h_cos_%s", lbl_cos[i].c_str()));
  }

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  
  // plot pT
  c->SetLogy();
  c->SetLeftMargin(0.11);
  int colpt[] = {kBlack, kRed+1};

  for(int i = 0; i < 2; i++) {

    h_pT[i]->SetStats(0);
    h_pT[i]->SetLineColor(colpt[i]);
    h_pT[i]->SetMinimum(5e1);
    h_pT[i]->SetMaximum(5e6);
    h_pT[i]->GetXaxis()->SetTitle("p_{T}(#mu#mu) (GeV)");
    h_pT[i]->GetXaxis()->SetTitleOffset(1.1);
    h_pT[i]->GetYaxis()->SetTitle("dN/dp_{T}");
    h_pT[i]->SetTitle("");
    if(i == 0) h_pT[i]->Draw("histo");
    else {
      h_pT[i]->Draw("histo same");
    }
  }
  
  TLatex lcpt;
  lcpt.SetTextSize(0.04);
  lcpt.SetTextColor(colpt[0]);
  lcpt.DrawLatex(120, 2e6, "Peak data");
  lcpt.SetTextColor(colpt[1]);
  lcpt.DrawLatex(120, 7e5, "MC");

  TLine *ptL = new TLine(100, 0, 100,  exp(0.5*(log(h_pT[0]->GetMaximum())+log(h_pT[0]->GetMinimum()))));
  ptL->SetLineStyle(kDashed);
  ptL->Draw();

  c->SaveAs(Form("plots/ANdists/pt_all.pdf"));
  c->Clear();

  // plot y
  c->SetLogy();
  int coly[] = {kBlack};
  for(int i = 0; i < 2; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(coly[i%1]);
    h_y[i]->SetMinimum(5e3);
    h_y[i]->SetMaximum(2e5);
    h_y[i]->GetXaxis()->SetTitle("y(#mu#mu)");
    h_y[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y[i]->GetYaxis()->SetTitle("dN/dy");
    h_y[i]->SetTitle("");
    if(i > 0) h_y[i]->SetLineStyle(kDashed);
    if(i == 0) h_y[i]->Draw("histo");
    else h_y[i]->Draw("histo same");  
  }

  TLatex lcy;
  lcy.SetTextSize(0.03);
  lcy.DrawLatex(0, 1e5, "#minus Peak data");
  lcy.DrawLatex(0, 7e4, "-- MC");
  
  c->SaveAs(Form("plots/ANdists/y_all.pdf"));
  c->Clear();

  // y scaling MC to data
  for(int i = 0; i < 2; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(coly[i%1]);
    h_y[i]->SetMinimum(5e3);
    h_y[i]->SetMaximum(2e5);
    h_y[i]->GetXaxis()->SetTitle("y(#mu#mu)");
    h_y[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y[i]->GetYaxis()->SetTitle("dN/dy");
    h_y[i]->SetTitle("");
    if(i > 0) {
      h_y[i]->Scale(h_y[i-1]->Integral() / h_y[i]->Integral());
      h_y[i]->SetLineStyle(kDashed);
    }
    if(i == 0) h_y[i]->Draw("histo");
    else h_y[i]->Draw("histo same");  
  }

  TLatex lcys;
  lcys.SetTextSize(0.03);
  lcys.DrawLatex(0, 6e4, "#minus Peak data");
  lcys.DrawLatex(0, 3.5e4, "-- scaled MC");

  c->SaveAs(Form("plots/ANdists/y_scale.pdf"));
  c->Clear();

  // plot lifetime
  c->SetLogy();

  // plot just full lt
  h_lt[0]->SetStats(0);
  h_lt[0]->SetLineColor(colpt[0]);
  h_lt[0]->SetMinimum(1e3);
  h_lt[0]->SetMaximum(2e5);
  h_lt[0]->GetXaxis()->SetTitle("c#tau (#mum)");
  h_lt[0]->GetXaxis()->SetTitleOffset(1.1);
  h_lt[0]->GetYaxis()->SetTitle("dN/dc#tau");
  h_lt[0]->SetTitle("");
  h_lt[0]->Draw("histo");
  
  TLine *l_PRm = new TLine(-50, h_lt[0]->GetMinimum(), -50, h_lt[0]->GetMaximum());
  l_PRm->SetLineStyle(kDashed);
  l_PRm->SetLineColor(kBlack);
  l_PRm->Draw();
  TLine *l_PRp = new TLine(50, h_lt[0]->GetMinimum(), 50, h_lt[0]->GetMaximum());
  l_PRp->SetLineStyle(kDashed);
  l_PRp->SetLineColor(kBlack);
  l_PRp->Draw();

  TLine *l_NPm = new TLine(100, h_lt[0]->GetMinimum(), 100, h_lt[0]->GetMaximum());
  l_NPm->SetLineStyle(kDotted);
  l_NPm->SetLineColor(kBlack);
  l_NPm->Draw();
  TLine *l_NPp = new TLine(500, h_lt[0]->GetMinimum(), 500, h_lt[0]->GetMaximum());
  l_NPp->SetLineStyle(kDotted);
  l_NPp->SetLineColor(kBlack);
  l_NPp->Draw();

  TLatex lcltf;
  lcltf.SetTextSize(0.04);
  lcltf.SetTextColor(colpt[0]);
  lcltf.DrawLatex(300, 1e5, "2018 Data");

  c->SaveAs(Form("plots/ANdists/lt_full.pdf"));
  c->Clear();

  // plot M
  c->SetLogy(0);
  c->SetLeftMargin(0.12);

  h_m[0]->SetStats(0);
  h_m[0]->SetLineColor(kBlack);
  h_m[0]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
  h_m[0]->GetXaxis()->SetTitleOffset(1.1);
  h_m[0]->GetYaxis()->SetTitle("dN/dM");
  h_m[0]->SetMaximum(130e3);
  h_m[0]->SetTitle("");
  h_m[0]->Draw("histo");  

  TLatex lcmf;
  lcmf.SetTextSize(0.04);
  lcmf.SetTextColor(colpt[0]);
  lcmf.DrawLatex(3.8, 1.2e5, "2018 Data");
  
  c->SaveAs(Form("plots/ANdists/m_full.pdf"));
  c->Clear();

  for(int i = 0; i < 2; i++) {
    h_m[i]->SetStats(0);
    h_m[i]->Scale(1./h_m[i]->Integral());
    h_m[i]->Scale(h_m[i%1]->GetMaximum()/h_m[i]->GetMaximum());
    h_m[i]->SetLineColor(coly[i%1]);
    h_m[i]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
    h_m[i]->GetXaxis()->SetTitleOffset(1.1);
    h_m[i]->GetYaxis()->SetTitle("dN/dM (normalized)");
    h_m[i]->SetTitle("");
    if(i > 0) h_m[i]->SetLineStyle(kDashed);
  }
  
  for(int i = 0; i < 2; i++) {
    h_m[i]->SetMinimum(0);
    h_m[i]->SetMaximum(0.07);
    if(i==0) h_m[i]->Draw("histo");
    else h_m[i]->Draw("histo same");
  }

  TLatex lcm;
  lcm.SetTextSize(0.04);
  lcm.SetTextColor(kBlack);
  lcm.DrawLatex(3.8, 0.066, "#minus Data");
  lcm.DrawLatex(3.8, 0.059, "--MC");

  c->SaveAs(Form("plots/ANdists/m_scale.pdf"));
  c->Clear();

  // plot costh
  c->SetLogy(0);

  int colc[] = {kViolet-1, kRed, kBlack};

  for(int i = 0; i < 3; i++) {

    h_cos[i]->SetStats(0);
    h_cos[i]->Scale(0.075/h_cos[i]->GetBinContent(1));
    h_cos[i]->SetLineColor(colc[i]);
    h_cos[i]->SetMinimum(0);
    h_cos[i]->SetMaximum(0.09);
    h_cos[i]->GetXaxis()->SetTitle("|cos#theta|");
    h_cos[i]->GetXaxis()->SetTitleOffset(1.);
    h_cos[i]->GetYaxis()->SetTitle("dN/d|cos#theta| (a.u.)");
    h_cos[i]->SetTitle("");
    if(i == 0) h_cos[i]->Draw("histo");
    else h_cos[i]->Draw("histo same"); 
  }

  TLatex lccf;
  lccf.SetTextSize(0.04);
  lccf.SetTextColor(colc[0]);
  lccf.DrawLatex(0.1, 0.04, "Data Peak");
  lccf.SetTextColor(colc[2]);
  lccf.DrawLatex(0.1, 0.033, "MC");
  lccf.SetTextColor(colc[1]);
  lccf.DrawLatex(0.1, 0.026, "Data NP");

  lccf.SetTextColor(colc[2]);
  lccf.SetTextSize(0.03);
  lccf.DrawLatex(0.08, 0.01, "20<p_{T}<100 GeV");

  c->SaveAs(Form("plots/ANdists/cos_full.pdf"));
  c->Clear();
  
  c->Destructor();
  fin->Close();
}
