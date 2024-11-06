// macro to take the stored histos and plot nicely

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


void plot_ANdists()
{
  // PART 1 : reading the histograms

  // four pT dists: PRSR data + the 3 MC
  TH1D **h_pT = new TH1D*[5]; 
  // 6 y dists: PRSR data + MC over 3 pT regions
  TH1D **h_y = new TH1D*[8]; 
  // 7 M dists: PR data + MC over 3 pT regions + full pT data
  TH1D **h_m = new TH1D*[9]; 
  // 4 lifetime dists: data in 4 pT regions
  TH1D **h_lt = new TH1D*[5];
  // 12 costh dists: PR data + NP data + MC over 4 pT regions
  TH1D **h_cos = new TH1D*[15]; 
  
  TFile *fin = new TFile("files/store_ANdists.root");
 
  // read the pT dists
  string lbl_pt[] = {"Data", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC"};
  for(int i = 0; i < 5; i++) {
    h_pT[i] = (TH1D*)fin->Get(Form("h_pT_%s", lbl_pt[i].c_str()));
  }

  // store the y dists
  string lbl_y[] = {"lowPtData", "midPtData", "highPtData", "highestPtData", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC"};
  for(int i = 0; i < 8; i++) {
    h_y[i] = (TH1D*)fin->Get(Form("h_y_%s", lbl_y[i].c_str()));
  }

  // store the M dists
  string lbl_m[] = {"lowPtData", "midPtData", "highPtData", "highestPtData", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC", "fullData"};
  for(int i = 0; i < 9; i++) {
    h_m[i] = (TH1D*)fin->Get(Form("h_m_%s", lbl_m[i].c_str()));
  }

  // store the lifetime dists
  string lbl_lt[] = {"lowPt", "midPt", "highPt", "highestPt", "full"};
  for(int i = 0; i < 5; i++) {
    h_lt[i] = (TH1D*)fin->Get(Form("h_lt_%s", lbl_lt[i].c_str()));
  }

  // store the cos dists
  string lbl_cos[] = {"lowPtPR", "midPtPR", "highPtPR", "highestPtPR", "fullPR", "lowPtNP", "midPtNP", "highPtNP", "highestPtNP", "fullNP", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC", "fullMC"};
  for(int i = 0; i < 15; i++) {
    h_cos[i] = (TH1D*)fin->Get(Form("h_cos_%s", lbl_cos[i].c_str()));
  }

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  
  // plot pT
  c->SetLogy();
  c->SetLeftMargin(0.11);
  int colpt[] = {kBlack, kRed+1, kViolet+2, kGreen+3, kBlue};

  for(int i = 0; i < 5; i++) {

    h_pT[i]->SetStats(0);
    h_pT[i]->SetLineColor(colpt[i]);
    h_pT[i]->SetMarkerColor(colpt[i]);
    h_pT[i]->SetMarkerStyle(20);
    h_pT[i]->SetMarkerSize(.5);
    h_pT[i]->SetMinimum(1e2);
    h_pT[i]->SetMaximum(2e7);
    h_pT[i]->GetXaxis()->SetRangeUser(0, 199);
    h_pT[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV)");
    h_pT[i]->GetXaxis()->SetTitleOffset(1.2);
    h_pT[i]->GetXaxis()->CenterTitle(true);
    h_pT[i]->GetXaxis()->SetLabelOffset(0.01);
    h_pT[i]->GetYaxis()->SetLabelOffset(0.01);
    h_pT[i]->GetYaxis()->SetTitle("dN/d#it{p}_{T}");
    h_pT[i]->SetTitle("");
    if(i == 0) {
      //    h_pT[i]->SetLineWidth(2);
      // h_pT[i]->SetFillColorAlpha(kYellow,0.2);
      h_pT[i]->Draw("error");
    }
    else {
      h_pT[i]->Draw("error same");
    }
  }
  
  TLatex lcpt;
  double yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.91, 1);
  lcpt.SetTextSize(0.04);
  lcpt.DrawLatex(45, yl, "#bf{2018 J/#psi}");
  lcpt.SetTextColor(colpt[0]);
  lcpt.DrawLatex(105, yl, "#bf{PRS data}");
  yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.83, 1);
  lcpt.DrawLatex(105, yl, "#bf{Prompt MC samples}");
  lcpt.SetTextColor(colpt[1]);
  yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.76, 1);
  lcpt.DrawLatex(120, yl, "#bf{[25,46] GeV}");
  lcpt.SetTextColor(colpt[2]);
  yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.69, 1);
  lcpt.DrawLatex(120, yl, "#bf{[40,52] GeV}");
  lcpt.SetTextColor(colpt[3]);
  yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.62, 1);
  lcpt.DrawLatex(120, yl, "#bf{> 46 GeV}");
  lcpt.SetTextColor(colpt[4]);
  yl = getPos(h_pT[0]->GetMinimum(), h_pT[0]->GetMaximum(), 0.55, 1);
  lcpt.DrawLatex(120, yl, "#bf{> 66 GeV}");

  //TLine *ptL = new TLine(120, 0, 120,  exp(0.5*(log(h_pT[0]->GetMaximum())+log(h_pT[0]->GetMinimum()))));
  TLine *ptL = new TLine(120, 0, 120,  2e4);
  ptL->SetLineStyle(kDashed);
  ptL->Draw();

  c->SaveAs(Form("plots/ANdists/pt_all.pdf"));
  c->Clear();

  // plot y
  c->SetLogy();
  int coly[] = {kRed+1, kViolet+2, kGreen+3, kBlue};
  for(int i = 0; i < 8; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(coly[i%4]);
    h_y[i]->SetMinimum(1e3);
    h_y[i]->SetMaximum(6e5);
    h_y[i]->GetXaxis()->SetTitle("#it{y}");
    h_y[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y[i]->GetXaxis()->CenterTitle(true);
    h_y[i]->GetXaxis()->SetLabelOffset(0.01);
    h_y[i]->GetYaxis()->SetLabelOffset(0.01);
    h_y[i]->GetYaxis()->SetTitle("dN/d#it{y}");
    h_y[i]->SetMarkerColor(coly[i%4]);
    h_y[i]->SetMarkerStyle(20);
    h_y[i]->SetMarkerSize(.5);
    h_y[i]->SetTitle("");
    if(i > 3) {
      h_y[i]->SetLineStyle(kDashed);
      h_y[i]->SetMarkerStyle(24);
    }
    if(i == 0) h_y[i]->Draw("error");
    else h_y[i]->Draw("error same");  
  }

  // just formatting h_m for the y legend
  h_m[1]->SetLineColor(kBlack);
  h_m[1]->SetMarkerColor(kBlack);
  h_m[1]->SetLineStyle(kSolid);
  h_m[1]->SetMarkerStyle(20);
  h_m[5]->SetLineColor(kBlack);
  h_m[5]->SetMarkerColor(kBlack);
  h_m[5]->SetLineStyle(kDashed);
  h_m[5]->SetMarkerStyle(24);
  
  TLegend *legcy = new TLegend(0.6, 0.65, 0.9, 0.775);
  legcy->SetTextSize(0.04);
  legcy->SetBorderSize(0);
  legcy->SetFillColorAlpha(kWhite,0);
  legcy->AddEntry(h_m[1], "PRS Data", "pl");
  legcy->AddEntry(h_m[5], "MC", "pl");
  legcy->Draw();

  TLatex lcy;
  lcy.SetTextSize(0.04);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.68, 1);
  lcy.DrawLatex(-1, yl, "#bf{2018 J/#psi}");
  
  //  c->SaveAs(Form("plots/ANdists/y_all.pdf"));
  c->Clear();

  // y scaling MC to data
  for(int i = 0; i < 8; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(coly[i%4]);
    h_y[i]->SetMinimum(1e3);
    h_y[i]->SetMaximum(6e5);
    h_y[i]->GetXaxis()->SetTitle("#it{y}");
    h_y[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y[i]->GetYaxis()->SetTitle("dN/d#it{y}");
    h_y[i]->SetTitle("");
    if(i > 3) {
      h_y[i]->Scale(h_y[i-4]->Integral() / h_y[i]->Integral());
      h_y[i]->SetLineStyle(kDashed);
    }
    if(i == 0) h_y[i]->Draw("error");
    else h_y[i]->Draw("error same");  
  }

  TLegend *legcys = new TLegend(0.575, 0.6, 0.875, 0.725);
  legcys->SetTextSize(0.04);
  legcys->SetBorderSize(0);
  legcys->SetFillColorAlpha(kWhite,0);
  legcys->AddEntry(h_m[1], "PRS Data", "pl");
  legcys->AddEntry(h_m[5], "Scaled MC", "pl");
  legcys->Draw();


  TLatex lcys;
  lcys.SetTextSize(0.04);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.9, 1);
  lcys.DrawLatex(-1, yl, "#bf{2018 J/#psi}");
  lcys.SetTextColor(coly[0]);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.7, 1);
  lcys.DrawLatex(-1, yl, "#bf{[25,45] GeV}");
  lcys.SetTextColor(coly[1]);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.62, 1);
  lcys.DrawLatex(-1, yl, "#bf{[45,50] GeV}");
  lcys.SetTextColor(coly[2]);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.54, 1);
  lcys.DrawLatex(-1, yl, "#bf{[50,70] GeV}");
  lcys.SetTextColor(coly[3]);
  yl = getPos(h_y[0]->GetMinimum(), h_y[0]->GetMaximum(), 0.46, 1);
  lcys.DrawLatex(-1, yl, "#bf{[70,120] GeV}");

  c->SaveAs(Form("plots/ANdists/y_scale.pdf"));
  c->Clear();

  // plot ratio data/MC
  c->SetLogy(0);
  TH1D **h_yr = new TH1D*[4]; 
  for(int i = 0; i < 4; i++) {
    h_yr[i] = (TH1D*)h_y[i]->Clone(Form("rH_%d",i));
    h_yr[i]->Sumw2();
    h_yr[i]->Divide(h_y[i+4]);

    h_yr[i]->SetMinimum(0.51);
    h_yr[i]->SetMaximum(1.49);
    h_yr[i]->SetStats(0);
    h_yr[i]->SetLineColor(coly[i]);
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
  lcyr.DrawLatex(-1.1, yl, "#bf{2018 J/#psi}");
  lcyr.SetTextColor(coly[0]);
  yl = getPos(h_yr[0]->GetMinimum(), h_yr[0]->GetMaximum(), 0.88, 0);
  lcyr.DrawLatex(0.35, yl, "#bf{[25,45] GeV}");
  lcyr.SetTextColor(coly[1]);
  yl = getPos(h_yr[0]->GetMinimum(), h_yr[0]->GetMaximum(), 0.8, 0);
  lcyr.DrawLatex(0.35, yl, "#bf{[45,50] GeV}");
  lcyr.SetTextColor(coly[2]);
  yl = getPos(h_yr[0]->GetMinimum(), h_yr[0]->GetMaximum(), 0.72, 0);
  lcyr.DrawLatex(0.35, yl, "#bf{[50,70] GeV}");
  lcyr.SetTextColor(coly[3]);
  yl = getPos(h_yr[0]->GetMinimum(), h_yr[0]->GetMaximum(), 0.64, 0);
  lcyr.DrawLatex(0.35, yl, "#bf{[70,120] GeV}");

  //  c->SaveAs(Form("plots/ANdists/y_ratio.pdf"));
  c->Clear();

  // plot lifetime
  c->SetLogy();

  for(int i = 0; i < 4; i++) {
    h_lt[i]->SetStats(0);
    h_lt[i]->SetLineColor(coly[i]);
    h_lt[i]->SetMarkerColor(coly[i]);
    h_lt[i]->SetMinimum(4e2);
    h_lt[i]->SetMaximum(3e6);
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
  lclt.DrawLatex(150, yl, "#bf{2018 J/#psi}");
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.95, 1);
  lclt.DrawLatex(400, yl, "#bf{Signal Region Data}");
  lclt.SetTextColor(colpt[1]);
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.88, 1);
  lclt.DrawLatex(500, yl, "#bf{[25,45] GeV}");
  lclt.SetTextColor(colpt[2]);
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.81, 1);
  lclt.DrawLatex(500, yl, "#bf{[45,50] GeV}");
  lclt.SetTextColor(colpt[3]);
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.74, 1);
  lclt.DrawLatex(500, yl, "#bf{[50,70] GeV}");
  lclt.SetTextColor(colpt[4]);
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.67, 1);
  lclt.DrawLatex(500, yl, "#bf{[70,120] GeV}");

  c->SaveAs(Form("plots/ANdists/lt_all.pdf"));
  c->Clear();

  // plot just full lt
  h_lt[4]->SetStats(0);
  h_lt[4]->SetLineColor(colpt[0]);
  h_lt[4]->SetMinimum(1e3);
  h_lt[4]->SetMaximum(4e6);
  h_lt[4]->GetXaxis()->SetTitle("c#tau (#mum)");
  h_lt[4]->GetXaxis()->SetTitleOffset(1.1);
  h_lt[4]->GetYaxis()->SetTitle("dN/dc#tau");
  h_lt[4]->SetTitle("");
  h_lt[4]->Draw("error");
  
  l_PRm->Draw();
  l_PRp->Draw();
  l_NPm->Draw();
  l_NPp->Draw();

  TLatex lcltf;
  lcltf.SetTextSize(0.04);
  lcltf.SetTextColor(colpt[0]);
  lcltf.DrawLatex(300, 1.5e6, "Run 2 Data");

  //  c->SaveAs(Form("plots/ANdists/lt_full.pdf"));
  c->Clear();

  // lifetime scaling to low pT
  for(int i = 0; i < 4; i++) {
    h_lt[i]->SetStats(0);
    h_lt[i]->SetLineColor(coly[i]);
    h_lt[i]->Scale(h_lt[0]->Integral()/h_lt[i]->Integral());
    h_lt[i]->SetMinimum(4e2);
    h_lt[i]->SetMaximum(3e6);
    h_lt[i]->GetXaxis()->SetTitle("c#tau (#mum)");
    h_lt[i]->GetXaxis()->SetTitleOffset(1.1);
    h_lt[i]->GetYaxis()->SetTitle("dN/d c#tau (a.u.)");
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
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.9, 1);
  lclts.DrawLatex(150, yl, "#bf{2018 J/#psi}");
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.5, 1);
  lclts.DrawLatex(125, yl, "#bf{Signal Region Data}");
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.1, 1);
  double xl = getPos(-50, 50, 0.5, 0);
  lclts.SetTextAlign(21);
  lclts.DrawLatex(xl, yl, "#bf{PR}");
  xl = getPos(100, 800, 0.5, 0);
  lclts.DrawLatex(xl, yl, "#bf{NP}");
  lclts.SetTextColor(colpt[1]);
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.43, 1);
  lclts.SetTextAlign(11);
  lclts.DrawLatex(225, yl, "#bf{[25,45] GeV}");
  lclts.SetTextColor(colpt[2]);
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.36, 1);
  lclts.DrawLatex(225, yl, "#bf{[45,50] GeV}");
  lclts.SetTextColor(colpt[3]);
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.29, 1);
  lclts.DrawLatex(225, yl, "#bf{[50,70] GeV}");
  lclts.SetTextColor(colpt[4]);
  yl = getPos(h_lt[0]->GetMinimum(), h_lt[0]->GetMaximum(), 0.22, 1);
  lclts.DrawLatex(225, yl, "#bf{[70,120] GeV}");

  c->SaveAs(Form("plots/ANdists/lt_scale.pdf"));
  c->Clear();

  // plot M
  c->SetLogy(0);
  c->SetLeftMargin(0.12);

  for(int i = 0; i < 8; i++) {
    h_m[i]->SetStats(0);
    h_m[i]->Scale(1./h_m[i]->Integral());
    h_m[i]->Scale(h_m[i%4]->GetMaximum()/h_m[i]->GetMaximum());
    h_m[i]->SetLineColor(coly[i%4]);
    h_m[i]->SetMarkerColor(coly[i%4]);
    h_m[i]->SetMarkerStyle(20);
    h_m[i]->SetMarkerSize(.5);
    h_m[i]->GetXaxis()->SetTitle("#it{m} (GeV)");
    h_m[i]->GetXaxis()->CenterTitle(true);
    h_m[i]->GetXaxis()->SetTitleOffset(1.1);
    h_m[i]->GetXaxis()->SetLabelOffset(0.01);
    h_m[i]->GetYaxis()->SetLabelOffset(0.01);
    h_m[i]->GetYaxis()->SetTitle("dN/d#it{m} (normalized)");
    h_m[i]->SetTitle("");
    if(i > 3){
      h_m[i]->SetLineStyle(kDashed);
      h_m[i]->SetMarkerStyle(24);
    }
  }
  
  for(int i = 0; i < 8; i++) {  
    if(i%4==1 || i%4 == 2) continue;
    
    h_m[i]->SetMinimum(0);
    h_m[i]->SetMaximum(0.085);
    if(i==0) h_m[i]->Draw("error");
    else h_m[i]->Draw("error same");
  }

  TLatex lcm;
  lcm.SetTextSize(0.04);
  lcm.DrawLatex(2.95, 0.074, "#bf{2018 J/#psi}");
  lcm.SetTextColor(coly[0]);
  lcm.DrawLatex(3.16, 0.049, "#bf{#it{p}_{T}: 25-45 GeV}");
  lcm.SetTextColor(coly[3]);
  lcm.DrawLatex(3.16, 0.041, "#bf{#it{p}_{T}: 70-120 GeV}");

  // set a line black just for the legend
  h_m[1]->SetLineColor(kBlack);
  h_m[1]->SetMarkerColor(kBlack);
  h_m[5]->SetLineColor(kBlack);
  h_m[5]->SetMarkerColor(kBlack);
   
  TLegend *legcm = new TLegend(0.725, 0.77, 1.025, 0.92);
  legcm->SetTextSize(0.04);
  legcm->SetBorderSize(0);
  legcm->SetFillColorAlpha(kWhite,0);
  legcm->AddEntry(h_m[1], "PR Data", "pl");
  legcm->AddEntry(h_m[5], "MC", "pl");
  legcm->Draw();
  
  c->SaveAs(Form("plots/ANdists/m_scale.pdf"));
  c->Clear();

  c->SetTopMargin(0.04);

  h_m[8]->SetStats(0);
  h_m[8]->SetLineColor(kBlack);
  h_m[8]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
  h_m[8]->GetXaxis()->SetTitleOffset(1.1);
  h_m[8]->GetYaxis()->SetTitle("dN/dM");
  h_m[8]->GetYaxis()->SetTitleOffset(1.8);
  h_m[8]->SetMaximum(1.2e6);
  h_m[8]->SetTitle("");
  h_m[8]->Draw("error");  

  TLatex lcmf;
  lcmf.SetTextSize(0.04);
  lcmf.SetTextColor(colpt[0]);
  lcmf.DrawLatex(3.15, 1e6, "Run 2 Data");
  
  //  c->SaveAs(Form("plots/ANdists/m_full.pdf"));
  c->Clear();

  // plot costh
  c->SetLogy(0);
  c->SetTopMargin(0.04);
  c->SetRightMargin(0.03);

  for(int i = 10; i < 14; i++) {

    h_cos[i]->SetStats(0);
    h_cos[i]->SetLineColor(coly[i-10]);
    h_cos[i]->SetMarkerColor(coly[i-10]);
    h_cos[i]->SetMarkerStyle(20);
    h_cos[i]->SetMarkerSize(.5);
    h_cos[i]->SetMinimum(0);
    h_cos[i]->Scale(1./h_cos[i]->GetBinContent(1));
    h_cos[i]->SetMaximum(h_cos[i]->GetBinContent(1)*1.5);
    h_cos[i]->GetXaxis()->SetTitle("|cos #theta_{HX}|");
    h_cos[i]->GetYaxis()->SetTitle("Events / 0.05 (normalized)");
    h_cos[i]->SetTitle("");
    h_cos[i]->GetYaxis()->SetMaxDigits(3);
    h_cos[i]->GetXaxis()->SetLabelOffset(0.01);
    h_cos[i]->GetYaxis()->SetLabelOffset(0.01);
    h_cos[i]->GetXaxis()->SetTitleOffset(1.3);
    h_cos[i]->GetXaxis()->CenterTitle(true);
    if(i == 10) h_cos[i]->Draw("error");
    else h_cos[i]->Draw("error same"); 
  }

  TLatex lcc;
  lcc.SetTextSize(0.04);
  // draw the state
  yl = getPos(h_cos[10]->GetMinimum(), h_cos[10]->GetMaximum(), 0.9, 0);
  lcc.DrawLatex(0.1, yl, "#bf{2018 J/#psi}");
  lcc.SetTextColor(coly[0]);
  yl = getPos(h_cos[10]->GetMinimum(), h_cos[10]->GetMaximum(), 0.46, 0);
  lcc.DrawLatex(0.1, yl, "#bf{[25,45] GeV}");
  lcc.SetTextColor(coly[1]);
  yl = getPos(h_cos[10]->GetMinimum(), h_cos[10]->GetMaximum(), 0.38, 0);
  lcc.DrawLatex(0.1, yl, "#bf{[45,50] GeV}");
  lcc.SetTextColor(coly[2]);
  yl = getPos(h_cos[10]->GetMinimum(), h_cos[10]->GetMaximum(), 0.3, 0);
  lcc.DrawLatex(0.1, yl, "#bf{[50,70] GeV}");
  lcc.SetTextColor(coly[3]);
  yl = getPos(h_cos[10]->GetMinimum(), h_cos[10]->GetMaximum(), 0.22, 0);
  lcc.DrawLatex(0.1, yl, "#bf{[70,120] GeV}");
  
  c->SaveAs(Form("plots/ANdists/cos_MC.pdf"));
  c->Clear();

  int colc[] = {kViolet+2, kRed+1, kBlack};

  for(int i = 0; i < 3; i++) {

    int j = (i+1)*5-1;
    
    h_cos[j]->SetStats(0);
    h_cos[j]->Scale(1./h_cos[j]->GetBinContent(1));
    h_cos[j]->SetLineColor(colc[i]);
    h_cos[j]->SetMinimum(0);
    h_cos[j]->SetMaximum(h_cos[j]->GetBinContent(1)*1.19);
    h_cos[j]->GetXaxis()->SetTitle("|cos #theta_{HX}|");
    h_cos[j]->GetYaxis()->SetTitle("Events / 0.05 (normalized)");
    h_cos[j]->SetTitle("");
    h_cos[j]->SetMarkerColor(colc[i]);
    h_cos[j]->SetMarkerStyle(20);
    h_cos[j]->SetMarkerSize(.5);
    h_cos[j]->GetYaxis()->SetMaxDigits(3);
    h_cos[j]->GetXaxis()->SetLabelOffset(0.01);
    h_cos[j]->GetYaxis()->SetLabelOffset(0.01);
    h_cos[j]->GetXaxis()->SetTitleOffset(1.3);
    h_cos[j]->GetXaxis()->CenterTitle(true);
    if(i == 0) h_cos[j]->Draw("error");
    else h_cos[j]->Draw("error same"); 


  }

  TLatex lccf;
  lccf.SetTextSize(0.04);
  // draw the state
  lccf.SetTextAlign(21);
  yl = getPos(h_cos[4]->GetMinimum(), h_cos[4]->GetMaximum(), 0.9, 0);
  lccf.DrawLatex(0.8, yl, "#bf{2018 J/#psi}"); 
  yl = getPos(h_cos[4]->GetMinimum(), h_cos[4]->GetMaximum(), 0.82, 0);
  lccf.DrawLatex(0.8, yl, "#bf{[25,120] GeV}");

  lccf.SetTextAlign(11);
  lccf.SetTextColor(colc[0]);
  yl = getPos(h_cos[4]->GetMinimum(), h_cos[4]->GetMaximum(), 0.46, 0);
  lccf.DrawLatex(0.1, yl, "#bf{Data PRS}");
  lccf.SetTextColor(colc[2]);
  yl = getPos(h_cos[4]->GetMinimum(), h_cos[4]->GetMaximum(), 0.38, 0);
  lccf.DrawLatex(0.1, yl, "#bf{MC}");
  lccf.SetTextColor(colc[1]);
  yl = getPos(h_cos[4]->GetMinimum(), h_cos[4]->GetMaximum(), 0.3, 0);
  lccf.DrawLatex(0.1, yl, "#bf{Data NPS}");


  c->SaveAs(Form("plots/ANdists/cos_full.pdf"));
  c->Clear();
  
  c->Destructor();
  fin->Close();
}
