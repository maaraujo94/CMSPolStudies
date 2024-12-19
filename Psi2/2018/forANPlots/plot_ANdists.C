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

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
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
  string lbl_lt[] = {"Data"};
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
  int colpt[] = {kBlack, kRed+1, kViolet+2, kGreen+3, kBlue};

  for(int i = 0; i < 2; i++) {

    h_pT[i]->SetStats(0);
    h_pT[i]->SetLineColor(colpt[i]);
    h_pT[i]->SetMarkerColor(colpt[i]);
    h_pT[i]->SetMarkerStyle(20);
    h_pT[i]->SetMarkerSize(.5);
    h_pT[i]->SetMinimum(5e1);
    h_pT[i]->SetMaximum(8e6);
    h_pT[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV)");
    h_pT[i]->GetXaxis()->SetTitleOffset(1.2);
    h_pT[i]->GetXaxis()->CenterTitle(true);
    h_pT[i]->GetXaxis()->SetLabelOffset(0.01);
    h_pT[i]->GetYaxis()->SetLabelOffset(0.01);
    h_pT[i]->GetXaxis()->SetRangeUser(0,170);
    h_pT[i]->GetYaxis()->SetTitle("dN/d#it{p}_{T}");
    h_pT[i]->SetTitle("");
    if(i == 0) {
      h_pT[i]->Draw("error");
      //      h_pT[i]->SetLineWidth(2);
      // h_pT[i]->SetFillColorAlpha(kYellow,0.2);
    }
    else {
      h_pT[i]->Draw("error same");
    }
  }
  
  TLatex lcpt;
  double yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.91, 1);
  lcpt.SetTextSize(0.04);
  lcpt.DrawLatex(35, yl, "#bf{2018 #psi(2S)}");
  lcpt.SetTextColor(colpt[0]);
  lcpt.DrawLatex(90, yl, "#bf{PRS data}");
  yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.83, 1);
  lcpt.SetTextColor(colpt[1]);
  lcpt.DrawLatex(90, yl, "#bf{Prompt MC sample}");

  TLine *ptL = new TLine(100, 0, 100,  2e3);
  ptL->SetLineStyle(kDashed);
  ptL->Draw();

  c->SaveAs(Form("plots/ANdists/pt_all.pdf"));
  c->Clear();

  // plot y
  c->SetLogy();
  int coly[] = {kBlack};
  for(int i = 0; i < 2; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(colpt[i]);
    h_y[i]->SetMarkerColor(colpt[i]);
    h_y[i]->SetMarkerStyle(20);
    h_y[i]->SetMarkerSize(.5);
    h_y[i]->SetMinimum(1e3);
    h_y[i]->SetMaximum(2e5);
    h_y[i]->GetXaxis()->SetTitle("#it{y}");
    h_y[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y[i]->GetXaxis()->CenterTitle(true);
    h_y[i]->GetXaxis()->SetLabelOffset(0.01);
    h_y[i]->GetYaxis()->SetLabelOffset(0.01);
    h_y[i]->GetYaxis()->SetTitle("dN/d#it{y}");
    h_y[i]->SetTitle("");
    if(i > 0) {
      //      h_y[i]->SetLineStyle(kDashed);
      h_y[i]->SetMarkerStyle(24);
    }
    if(i == 0) h_y[i]->Draw("error");
    else h_y[i]->Draw("error same");  
  }

  TLegend *legcy = new TLegend(0.6, 0.725, 0.9, 0.85);
  legcy->SetTextSize(0.04);
  legcy->SetBorderSize(0);
  legcy->SetFillColorAlpha(kWhite,0);
  legcy->AddEntry(h_y[0], "PRS Data", "pl");
  legcy->AddEntry(h_y[1], "MC", "pl");
  legcy->Draw();

  TLatex lcy;
  lcy.SetTextSize(0.04);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.76, 1);
  lcy.DrawLatex(-1, yl, "#bf{2018 #psi(2S)}");

  
  //  c->SaveAs(Form("plots/ANdists/y_all.pdf"));
  c->Clear();

  // y scaling MC to data
  for(int i = 0; i < 2; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(colpt[i]);
    h_y[i]->SetMinimum(1e3);
    h_y[i]->SetMaximum(2e5);
    h_y[i]->SetTitle("");
    if(i > 0) {
      h_y[i]->Scale(h_y[i-1]->Integral() / h_y[i]->Integral());
    }
    if(i == 0) h_y[i]->Draw("error");
    else h_y[i]->Draw("error same");  
  }

  TLegend *legcys = new TLegend(0.575, 0.475, 0.875, 0.6);
  legcys->SetTextSize(0.04);
  legcys->SetBorderSize(0);
  legcys->SetFillColorAlpha(kWhite,0);
  legcys->AddEntry(h_y[0], "PRS Data", "pl");
  legcys->AddEntry(h_y[1], "Scaled MC", "pl");
  legcys->Draw();

  TLatex lcys;
  lcys.SetTextSize(0.04);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.9, 1);
  lcys.DrawLatex(-1, yl, "#bf{2018 #psi(2S)}");
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.475, 1);
  lcys.DrawLatex(-1, yl, "#bf{[20,100] GeV}");

  c->SaveAs(Form("plots/ANdists/y_scale.pdf"));
  c->Clear();

    // plot ratio data/MC
  c->SetLogy(0);
  TH1D **h_yr = new TH1D*[1]; 
  for(int i = 0; i < 1; i++) {
    h_yr[i] = (TH1D*)h_y[i]->Clone(Form("rH_%d",i));
    h_yr[i]->Sumw2();
    h_yr[i]->Divide(h_y[i+1]);

    h_yr[i]->SetMinimum(0.51);
    h_yr[i]->SetMaximum(1.49);
    h_yr[i]->SetStats(0);
    h_yr[i]->SetLineColor(colpt[i]);
    h_yr[i]->GetXaxis()->SetTitle("#it{y}");
    h_yr[i]->GetXaxis()->SetTitleOffset(1.1);
    h_yr[i]->GetXaxis()->CenterTitle(true);
    h_yr[i]->GetXaxis()->SetLabelOffset(0.01);
    h_yr[i]->GetYaxis()->SetLabelOffset(0.01);
    h_yr[i]->GetYaxis()->SetTitle("Data/MC");
    h_yr[i]->SetTitle("");

    if(i == 0) h_yr[i]->Draw("error");
    else h_yr[i]->Draw("error same");  
  }

  TLatex lcyr;
  lcyr.SetTextSize(0.04);
  yl = getPos(h_yr[0]->GetMinimum(), h_yr[0]->GetMaximum(), 0.9, 0);
  lcyr.DrawLatex(-1.1, yl, "#bf{2018 #psi(2S)}");
  yl = getPos(h_yr[0]->GetMinimum(), h_yr[0]->GetMaximum(), 0.9, 0);
  lcyr.DrawLatex(0.35, yl, "#bf{[20,100] GeV}");

  //  c->SaveAs(Form("plots/ANdists/y_ratio.pdf"));
  c->Clear();

  // plot lifetime
  c->SetLogy();

  for(int i = 0; i < 1; i++) {
    h_lt[i]->SetStats(0);
    h_lt[i]->SetLineColor(coly[i]);
    h_lt[i]->SetMinimum(1e3);
    h_lt[i]->SetMaximum(4e5);
    h_lt[i]->GetXaxis()->SetTitle("c#tau (#mum)");
    h_lt[i]->GetXaxis()->SetTitleOffset(1.1);
    h_lt[i]->GetXaxis()->CenterTitle(true);
    h_lt[i]->GetXaxis()->SetLabelOffset(0.01);
    h_lt[i]->GetYaxis()->SetLabelOffset(0.01);
    h_lt[i]->GetYaxis()->SetTitle("dN/d c#tau");
    h_lt[i]->SetTitle("");
    if(i == 0) h_lt[i]->Draw("error");
    else h_lt[i]->Draw("error same");  
  }

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
  TLine *l_NPp = new TLine(800, h_lt[0]->GetMinimum(), 800, h_lt[0]->GetMaximum());
  l_NPp->SetLineStyle(kDotted);
  l_NPp->SetLineColor(kBlack);
  l_NPp->Draw();

  TLatex lclt;
  lclt.SetTextSize(0.04);
  lclt.SetTextColor(colpt[0]);
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.9, 1);
  lclt.DrawLatex(125, yl, "#bf{2018 #psi(2S)}");
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.35, 1);
  lclt.DrawLatex(125, yl, "#bf{Signal Region Data}");
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.28, 1);
  lclt.SetTextAlign(11);
  lclt.DrawLatex(225, yl, "#bf{[20,100] GeV}");

  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.1, 1);
  double xl = getPos(-50, 50, 0.5, 0);
  lclt.SetTextAlign(21);
  lclt.DrawLatex(xl, yl, "#bf{PR}");
  xl = getPos(100, 800, 0.5, 0);
  lclt.DrawLatex(xl, yl, "#bf{NP}");

  c->SaveAs(Form("plots/ANdists/lt_all.pdf"));
  c->Clear();

  // plot just full lt
  h_lt[0]->SetStats(0);
  h_lt[0]->SetLineColor(colpt[0]);
  h_lt[0]->SetMinimum(1e3);
  h_lt[0]->SetMaximum(4e6);
  h_lt[0]->GetXaxis()->SetTitle("c#tau (#mum)");
  h_lt[0]->GetXaxis()->SetTitleOffset(1.1);
  h_lt[0]->GetYaxis()->SetTitle("dN/dc#tau");
  h_lt[0]->SetTitle("");
  h_lt[0]->Draw("error");
  
  l_PRm->Draw();
  l_PRp->Draw();
  l_NPm->Draw();
  l_NPp->Draw();

  TLatex lcltf;
  lcltf.SetTextSize(0.04);
  lcltf.SetTextColor(colpt[0]);
  lcltf.DrawLatex(300, 1.5e6, "2018 Data");

  //  c->SaveAs(Form("plots/ANdists/lt_full.pdf"));
  c->Clear();

  // lifetime scaling to low pT
  for(int i = 0; i < 1; i++) {
    h_lt[i]->SetStats(0);
    h_lt[i]->SetLineColor(coly[i]);
    h_lt[i]->Scale(h_lt[0]->Integral()/h_lt[i]->Integral());
    h_lt[i]->SetMinimum(1e3);
    h_lt[i]->SetMaximum(4e6);
    h_lt[i]->GetXaxis()->SetTitle("c#tau (#mum)");
    h_lt[i]->GetXaxis()->SetTitleOffset(1.1);
    h_lt[i]->GetYaxis()->SetTitle("dN/dc#tau (a.u.)");
    h_lt[i]->SetTitle("");
    if(i == 0) h_lt[i]->Draw("error");
    else h_lt[i]->Draw("error same");  
  }

  l_PRm->Draw();
  l_PRp->Draw();
  l_NPm->Draw();
  l_NPp->Draw();

  TLatex lclts;
  lclts.SetTextSize(0.04);
  lclts.SetTextColor(colpt[0]);
  lclts.DrawLatex(300, 2.5e6, "Data:");
  lclts.DrawLatex(-40, 2e3, "Peak");
  lclts.DrawLatex(275, 2e3, "NP");

  //  c->SaveAs(Form("plots/ANdists/lt_scale.pdf"));
  c->Clear();

  // plot M
  c->SetLogy(0);
  c->SetLeftMargin(0.13);

  for(int i = 0; i < 2; i++) {
    h_m[i]->SetStats(0);
    h_m[i]->Scale(1./h_m[i]->Integral());
    h_m[i]->Scale(h_m[0]->GetMaximum()/h_m[i]->GetMaximum());
    h_m[i]->SetLineColor(colpt[i]);
    h_m[i]->SetMarkerColor(colpt[i]);
    h_m[i]->SetMarkerStyle(20);
    h_m[i]->SetMarkerSize(.5);
    h_m[i]->GetXaxis()->SetTitle("#it{m} (GeV)");
    h_m[i]->GetXaxis()->CenterTitle(true);
    h_m[i]->GetXaxis()->SetTitleOffset(1.1);
    h_m[i]->GetYaxis()->SetTitleOffset(2.);
    h_m[i]->GetXaxis()->SetLabelOffset(0.01);
    h_m[i]->GetYaxis()->SetLabelOffset(0.01);
    h_m[i]->GetYaxis()->SetTitle("dN/d#it{m} (normalized)");
    h_m[i]->SetTitle("");
    if(i > 0){
      //      h_m[i]->SetLineStyle(kDashed);
      h_m[i]->SetMarkerStyle(24);
    }
  }
  
  for(int i = 0; i < 2; i++) {  
    h_m[i]->SetMinimum(0);
    h_m[i]->SetMaximum(0.0425);
    if(i==0) h_m[i]->Draw("error");
    else h_m[i]->Draw("error same");
  }


  TLatex lcm;
  lcm.SetTextSize(0.04);
  lcm.DrawLatex(3.425, 0.037, "#bf{2018 #psi(2S)}");
  lcm.DrawLatex(3.75, 0.0245, "#bf{#it{p}_{T}: 20-100 GeV}");

  TLegend *legcm = new TLegend(0.725, 0.77, 1.025, 0.92);
  legcm->SetTextSize(0.04);
  legcm->SetBorderSize(0);
  legcm->SetFillColorAlpha(kWhite,0);
  legcm->AddEntry(h_m[0], "PR Data", "pl");
  legcm->AddEntry(h_m[1], "MC", "pl");
  legcm->Draw();

  c->SaveAs(Form("plots/ANdists/m_scale.pdf"));
  c->Clear();

  c->SetTopMargin(0.04);

  h_m[0]->SetStats(0);
  h_m[0]->SetLineColor(kBlack);
  h_m[0]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
  h_m[0]->GetXaxis()->SetTitleOffset(1.1);
  h_m[0]->GetYaxis()->SetTitle("dN/dM");
  h_m[0]->GetYaxis()->SetTitleOffset(1.8);
  h_m[0]->SetMaximum(1.2e6);
  h_m[0]->SetTitle("");
  h_m[0]->Draw("error");  

  TLatex lcmf;
  lcmf.SetTextSize(0.04);
  lcmf.SetTextColor(colpt[0]);
  lcmf.DrawLatex(3.15, 1e6, "2018 Data");
  
  //  c->SaveAs(Form("plots/ANdists/m_full.pdf"));
  c->Clear();

  // plot costh
  c->SetLogy(0);
  c->SetTopMargin(0.015);

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
    if(i == 0) h_cos[i]->Draw("error");
    else h_cos[i]->Draw("error same"); 
  }

  TLatex lccf;
  lccf.SetTextSize(0.04);
  lccf.SetTextColor(colc[0]);
  lccf.DrawLatex(0.1, 0.04, "Data PRS");
  lccf.SetTextColor(colc[2]);
  lccf.DrawLatex(0.1, 0.033, "MC");
  lccf.SetTextColor(colc[1]);
  lccf.DrawLatex(0.1, 0.026, "Data NP");

  lccf.SetTextColor(colc[2]);
  lccf.SetTextSize(0.03);
  lccf.DrawLatex(0.08, 0.01, "20<p_{T}<100 GeV");

  //  c->SaveAs(Form("plots/ANdists/cos_full.pdf"));
  c->Clear();
  
  c->Destructor();
  fin->Close();
}
