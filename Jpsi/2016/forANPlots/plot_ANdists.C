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

  // 3 pT dists: PRSR data + the 2 MC
  TH1D **h_pT = new TH1D*[3]; 
  // 3 y dists: PRSR data + the 2 MC
  TH1D **h_y = new TH1D*[3]; 
  // 3 M dists: PR data + the 2 MC
  TH1D **h_m = new TH1D*[3]; 
  // 1 lifetime dist: data
  TH1D **h_lt = new TH1D*[1];
  // 4 costh dists: PR data + NP data + the 2 MC
  TH1D **h_cos = new TH1D*[4]; 
  
  TFile *fin = new TFile("files/store_ANdists.root");
 
  // read the pT dists
  string lbl_pt[] = {"Data", "MC1", "MC2"};
  for(int i = 0; i < 3; i++) {
    h_pT[i] = (TH1D*)fin->Get(Form("h_pT_%s", lbl_pt[i].c_str()));
  }

  // store the y dists
  string lbl_y[] = {"Data", "MC1", "MC2"};
  for(int i = 0; i < 3; i++) {
    h_y[i] = (TH1D*)fin->Get(Form("h_y_%s", lbl_y[i].c_str()));
  }

  // store the M dists
  string lbl_m[] = {"Data", "MC1", "MC2"};
  for(int i = 0; i < 3; i++) {
    h_m[i] = (TH1D*)fin->Get(Form("h_m_%s", lbl_m[i].c_str()));
  }

  // store the lifetime dists
  string lbl_lt[] = {"full"};
  for(int i = 0; i < 1; i++) {
    h_lt[i] = (TH1D*)fin->Get(Form("h_lt_%s", lbl_lt[i].c_str()));
  }

  // store the cos dists
  string lbl_cos[] = {"PR", "NP", "MC1", "MC2"};
  for(int i = 0; i < 4; i++) {
    h_cos[i] = (TH1D*)fin->Get(Form("h_cos_%s", lbl_cos[i].c_str()));
  }

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  
  // plot pT
  c->SetLogy();
  c->SetLeftMargin(0.11);
  int colpt[] = {kBlack, kRed+1, kViolet+2, kGreen+3, kBlue};

  for(int i = 0; i < 3; i++) {

    h_pT[i]->SetStats(0);
    h_pT[i]->SetLineColor(colpt[i]);
    h_pT[i]->SetMinimum(1e1);
    h_pT[i]->SetMaximum(1e7);
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
  lcpt.DrawLatex(120, 3.5e6, "Peak data");
  lcpt.SetTextColor(colpt[1]);
  lcpt.DrawLatex(120, 1.5e6, "MC 1");
  lcpt.SetTextColor(colpt[2]);
  lcpt.DrawLatex(120, 6e5, "MC 2");

  TLine *ptL = new TLine(80, 0, 80,  exp(0.5*(log(h_pT[0]->GetMaximum())+log(h_pT[0]->GetMinimum()))));
  ptL->SetLineStyle(kDashed);
  ptL->Draw();

  c->SaveAs(Form("plots/ANdists/pt_all.pdf"));
  c->Clear();
  
  // plot y
  c->SetLogy();
  int coly[] = {kRed+1, kViolet+2, kGreen+3, kBlue};
  for(int i = 0; i < 3; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(coly[i]);
    h_y[i]->SetMinimum(1e3);
    h_y[i]->SetMaximum(1e6);
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
  for(int i = 0; i < 3; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(coly[i]);
    h_y[i]->SetMinimum(1e3);
    h_y[i]->SetMaximum(1e6);
    h_y[i]->GetXaxis()->SetTitle("y(#mu#mu)");
    h_y[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y[i]->GetYaxis()->SetTitle("dN/dy");
    h_y[i]->SetTitle("");
    if(i > 0) {
      h_y[i]->Scale(h_y[0]->Integral() / h_y[i]->Integral());
      h_y[i]->SetLineStyle(kDashed);
    }
    if(i == 0) h_y[i]->Draw("histo");
    else h_y[i]->Draw("histo same");  
  }

  TLatex lcys;
  lcys.SetTextSize(0.03);
  lcys.DrawLatex(0, 6e4, "#minus Peak data");
  lcys.DrawLatex(0, 3.5e4, "-- scaled MC");

  lcys.SetTextColor(coly[1]);
  lcys.DrawLatex(-1, 6e4, "MC 1");
  lcys.SetTextColor(coly[2]);
  lcys.DrawLatex(-1, 3.5e4, "MC 2");

  c->SaveAs(Form("plots/ANdists/y_scale.pdf"));
  c->Clear();
  
  // plot lifetime
  c->SetLogy();

  // plot just full lt
  h_lt[0]->SetStats(0);
  h_lt[0]->SetLineColor(colpt[0]);
  h_lt[0]->SetMinimum(1e3);
  h_lt[0]->SetMaximum(3e6);
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
  lcltf.DrawLatex(300, 1.5e6, "2016 Data");

  c->SaveAs(Form("plots/ANdists/lt_full.pdf"));
  c->Clear();
  
  // plot M
  c->SetLogy(0);
  c->SetLeftMargin(0.12);

  double n_m = h_m[0]->Integral();  
  for(int i = 0; i < 3; i++) {
    h_m[i]->SetStats(0);
    h_m[i]->Scale(1./h_m[i]->Integral());
    h_m[i]->Scale(h_m[0]->GetMaximum()/h_m[i]->GetMaximum());
    h_m[i]->SetLineColor(coly[i]);
    h_m[i]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
    h_m[i]->GetXaxis()->SetTitleOffset(1.1);
    h_m[i]->GetYaxis()->SetTitle("dN/dM (normalized)");
    h_m[i]->SetTitle("");
    if(i > 0) h_m[i]->SetLineStyle(kDashed);
    if(i==0) h_m[i]->Draw("histo");
    else h_m[i]->Draw("histo same");
  }
  
  TLatex lcm;
  lcm.SetTextSize(0.04);
  lcm.SetTextColor(kBlack);
  lcm.DrawLatex(2.925, 0.073, "#minus Data");
  lcm.DrawLatex(2.925, 0.066, "--MC");

  c->SaveAs(Form("plots/ANdists/m_scale.pdf"));
  c->Clear();

  h_m[0]->Scale(n_m);
  h_m[0]->SetStats(0);
  h_m[0]->SetLineColor(kBlack);
  h_m[0]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
  h_m[0]->GetXaxis()->SetTitleOffset(1.1);
  h_m[0]->GetYaxis()->SetTitle("dN/dM");
  //h_m[0]->SetMaximum(720e3);
  h_m[0]->SetTitle("");
  h_m[0]->Draw("histo");  

  TLatex lcmf;
  lcmf.SetTextSize(0.04);
  lcmf.SetTextColor(colpt[0]);
  lcmf.DrawLatex(3.15, 6.5e5, "2016 Data");

  
  c->SaveAs(Form("plots/ANdists/m_full.pdf"));
  c->Clear();
  
  // plot costh
  c->SetLogy(0);

  for(int i = 2; i < 4; i++) {

    h_cos[i]->SetStats(0);
    h_cos[i]->Scale(0.06/h_cos[i]->GetBinContent(1));
    h_cos[i]->SetLineColor(coly[i-1]);
    h_cos[i]->SetMinimum(0);
    h_cos[i]->SetMaximum(0.09);
    h_cos[i]->GetXaxis()->SetTitle("|cos#theta|");
    h_cos[i]->GetXaxis()->SetTitleOffset(1.);
    h_cos[i]->GetYaxis()->SetTitle("dN/d|cos#theta| (a.u.)");
    h_cos[i]->SetTitle("");
    if(i == 2) h_cos[i]->Draw("histo");
    else h_cos[i]->Draw("histo same"); 
  }

  TLatex lcc;
  lcc.SetTextSize(0.04);
  lcc.SetTextColor(coly[1]);
  lcc.DrawLatex(0.05, 0.03, "MC 1");
  lcc.SetTextColor(coly[2]);
  lcc.DrawLatex(0.05, 0.023, "MC 2");
  
  c->SaveAs(Form("plots/ANdists/cos_MC.pdf"));
  c->Clear();
  /*
  int colc[] = {kViolet-1, kRed, kBlack};

  for(int i = 0; i < 3; i++) {

    int j = (i+1)*5-1;
    
    h_cos[j]->SetStats(0);
    h_cos[j]->Scale(0.075/h_cos[j]->GetBinContent(1));
    h_cos[j]->SetLineColor(colc[i]);
    h_cos[j]->SetMinimum(0);
    h_cos[j]->SetMaximum(0.09);
    h_cos[j]->GetXaxis()->SetTitle("|cos#theta|");
    h_cos[j]->GetXaxis()->SetTitleOffset(1.);
    h_cos[j]->GetYaxis()->SetTitle("dN/d|cos#theta| (a.u.)");
    h_cos[j]->SetTitle("");
    if(i == 0) h_cos[j]->Draw("histo");
    else h_cos[j]->Draw("histo same"); 
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
  lccf.DrawLatex(0.08, 0.01, "25<p_{T}<120 GeV");

  c->SaveAs(Form("plots/ANdists/cos_full.pdf"));
  c->Clear();
  */
  c->Destructor();
  fin->Close();
}
