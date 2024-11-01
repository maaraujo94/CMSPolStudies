// macro to take the stored histos and plot nicely

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}


void plot_ANdists()
{
  string loc = "/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi";
  // PART 1 : reading the histograms

  // four pT dists: PRSR data + the 3 MC
  TH1D **h_pT17 = new TH1D*[5]; 
  TH1D **h_pT18 = new TH1D*[5]; 
  // 6 y dists: PRSR data + MC over 3 pT regions
  TH1D **h_y17 = new TH1D*[8]; 
  TH1D **h_y18 = new TH1D*[8]; 
  // 7 M dists: PR data + MC over 3 pT regions + full pT data
  TH1D **h_m17 = new TH1D*[9]; 
  TH1D **h_m18 = new TH1D*[9]; 
  // 4 lifetime dists: data in 4 pT regions
  TH1D **h_lt17 = new TH1D*[5];
  TH1D **h_lt18 = new TH1D*[5];

  // get 2017 plots
  TFile *fin17 = new TFile(Form("%s/2017/forANPlots/files/store_ANdists.root", loc.c_str()));
  // get 2018 plots
  TFile *fin18 = new TFile(Form("%s/2018/forANPlots/files/store_ANdists.root", loc.c_str()));
 
  // read the pT dists
  string lbl_pt[] = {"Data", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC"};
  for(int i = 0; i < 5; i++) {
    h_pT17[i] = (TH1D*)fin17->Get(Form("h_pT_%s", lbl_pt[i].c_str()));
    h_pT18[i] = (TH1D*)fin18->Get(Form("h_pT_%s", lbl_pt[i].c_str()));
  }

  // store the y dists
  string lbl_y[] = {"lowPtData", "midPtData", "highPtData", "highestPtData", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC"};
  for(int i = 0; i < 8; i++) {
    h_y17[i] = (TH1D*)fin17->Get(Form("h_y_%s", lbl_y[i].c_str()));
    h_y18[i] = (TH1D*)fin18->Get(Form("h_y_%s", lbl_y[i].c_str()));
  }

  // store the M dists
  string lbl_m[] = {"lowPtData", "midPtData", "highPtData", "highestPtData", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC", "fullData"};
  for(int i = 0; i < 9; i++) {
    h_m17[i] = (TH1D*)fin17->Get(Form("h_m_%s", lbl_m[i].c_str()));
    h_m18[i] = (TH1D*)fin18->Get(Form("h_m_%s", lbl_m[i].c_str()));
  }

  // store the lifetime dists
  string lbl_lt[] = {"lowPt", "midPt", "highPt", "highestPt", "full"};
  for(int i = 0; i < 5; i++) {
    h_lt17[i] = (TH1D*)fin17->Get(Form("h_lt_%s", lbl_lt[i].c_str()));
    h_lt18[i] = (TH1D*)fin18->Get(Form("h_lt_%s", lbl_lt[i].c_str()));
  }

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetTopMargin(0.015);
  c->SetRightMargin(0.03);
  
  // plot pT
  c->SetLogy();
  c->SetLeftMargin(0.11);
  int colpt[] = {kBlack, kRed+1, kViolet+2, kGreen+3, kBlue};

  for(int i = 0; i < 5; i++) {

    h_pT17[i]->Scale(h_pT18[i]->Integral() / h_pT17[i]->Integral());
    h_pT17[i]->SetStats(0);
    h_pT17[i]->SetLineColor(colpt[i]);
    h_pT17[i]->SetMarkerColor(colpt[i]);
    h_pT17[i]->SetMarkerStyle(24);
    h_pT17[i]->SetMarkerSize(.5);
    h_pT17[i]->SetLineStyle(kDashed);
    h_pT17[i]->SetMinimum(1e2);
    h_pT17[i]->SetMaximum(2e7);
    h_pT17[i]->GetXaxis()->SetRangeUser(0, 199);
    h_pT17[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV)");
    h_pT17[i]->GetXaxis()->SetTitleOffset(1.2);
    h_pT17[i]->GetXaxis()->CenterTitle(true);
    h_pT17[i]->GetXaxis()->SetLabelOffset(0.01);
    h_pT17[i]->GetYaxis()->SetLabelOffset(0.01);
    h_pT17[i]->GetYaxis()->SetTitle("dN/d#it{p}_{T}");
    h_pT17[i]->SetTitle("");
    if(i == 0) {
      h_pT17[i]->Draw("error");
    }
    else {
      h_pT17[i]->Draw("error same");
    }

    h_pT18[i]->SetLineColor(colpt[i]);
    h_pT18[i]->SetMarkerColor(colpt[i]);
    h_pT18[i]->SetMarkerStyle(20);
    h_pT18[i]->SetMarkerSize(.5);
    h_pT18[i]->Draw("error same");
      
  }

  TLegend *legcpt = new TLegend(0.175, 0.2, 0.475, 0.325);
  legcpt->SetTextSize(0.04);
  legcpt->SetBorderSize(0);
  legcpt->SetFillColorAlpha(kWhite,0);
  legcpt->AddEntry(h_pT17[0], "2017", "pl");
  legcpt->AddEntry(h_pT18[0], "2018", "pl");
  legcpt->Draw();

  TLatex lcpt;
  double yl = getPos(h_pT17[0]->GetMinimum(), h_pT17[0]->GetMaximum(), 0.91, 1);
  lcpt.SetTextSize(0.04);
  lcpt.DrawLatex(45, yl, "#bf{J/#psi}");
  lcpt.SetTextColor(colpt[0]);
  lcpt.DrawLatex(105, yl, "#bf{PRS data}");
  yl = getPos(h_pT17[0]->GetMinimum(), h_pT17[0]->GetMaximum(), 0.83, 1);
  lcpt.DrawLatex(105, yl, "#bf{Prompt MC samples}");
  lcpt.SetTextColor(colpt[1]);
  yl = getPos(h_pT17[0]->GetMinimum(), h_pT17[0]->GetMaximum(), 0.76, 1);
  lcpt.DrawLatex(120, yl, "#bf{[25,46] GeV}");
  lcpt.SetTextColor(colpt[2]);
  yl = getPos(h_pT17[0]->GetMinimum(), h_pT17[0]->GetMaximum(), 0.69, 1);
  lcpt.DrawLatex(120, yl, "#bf{[40,52] GeV}");
  lcpt.SetTextColor(colpt[3]);
  yl = getPos(h_pT17[0]->GetMinimum(), h_pT17[0]->GetMaximum(), 0.62, 1);
  lcpt.DrawLatex(120, yl, "#bf{> 46 GeV}");
  lcpt.SetTextColor(colpt[4]);
  yl = getPos(h_pT17[0]->GetMinimum(), h_pT17[0]->GetMaximum(), 0.55, 1);
  lcpt.DrawLatex(120, yl, "#bf{> 66 GeV}");

  TLine *ptL = new TLine(120, 0, 120,  2e4);
  ptL->SetLineStyle(kDashed);
  ptL->Draw();

  c->SaveAs(Form("plots/pt_all.pdf"));
  c->Clear();

  // plot y
  c->SetLogy();
  int coly[] = {kRed+1, kViolet+2, kGreen+3, kBlue};
  for(int i = 0; i < 8; i++) {

    h_y17[i]->SetStats(0);
    h_y17[i]->SetLineColor(coly[i%4]);
    h_y17[i]->SetMinimum(1e3);
    h_y17[i]->SetMaximum(6e5);
    h_y17[i]->GetXaxis()->SetTitle("#it{y}");
    h_y17[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y17[i]->GetXaxis()->CenterTitle(true);
    h_y17[i]->GetXaxis()->SetLabelOffset(0.01);
    h_y17[i]->GetYaxis()->SetLabelOffset(0.01);
    h_y17[i]->GetYaxis()->SetTitle("dN/d#it{y}");
    h_y17[i]->SetMarkerColor(coly[i%4]);
    h_y17[i]->SetMarkerStyle(20);
    h_y17[i]->SetMarkerSize(.5);
    h_y17[i]->SetTitle("");
    if(i > 3) {
      h_y17[i]->SetLineStyle(kDashed);
      h_y17[i]->SetMarkerStyle(24);
    }
    if(i == 0) h_y17[i]->Draw("error");
    else h_y17[i]->Draw("error same");  
  }

  // just formatting h_m for the y legend
  h_m17[1]->SetLineColor(kBlack);
  h_m17[1]->SetMarkerColor(kBlack);
  h_m17[1]->SetLineStyle(kSolid);
  h_m17[1]->SetMarkerStyle(20);
  h_m17[5]->SetLineColor(kBlack);
  h_m17[5]->SetMarkerColor(kBlack);
  h_m17[5]->SetLineStyle(kDashed);
  h_m17[5]->SetMarkerStyle(24);
  
  TLegend *legcy = new TLegend(0.6, 0.65, 0.9, 0.775);
  legcy->SetTextSize(0.04);
  legcy->SetBorderSize(0);
  legcy->SetFillColorAlpha(kWhite,0);
  legcy->AddEntry(h_m17[1], "PRS Data", "pl");
  legcy->AddEntry(h_m17[5], "MC", "pl");
  legcy->Draw();

  TLatex lcy;
  lcy.SetTextSize(0.04);
  yl = getPos(h_y17[0]->GetMinimum(), h_y17[0]->GetMaximum(), 0.68, 1);
  lcy.DrawLatex(-1, yl, "#bf{2018 J/#psi}");
  
  //  c->SaveAs(Form("plots/y_all.pdf"));
  c->Clear();

  // y scaling MC to data
  for(int i = 0; i < 8; i++) {

    h_y17[i]->SetStats(0);
    h_y17[i]->SetLineColor(coly[i%4]);
    h_y17[i]->SetMinimum(1e3);
    h_y17[i]->SetMaximum(6e5);
    h_y17[i]->GetXaxis()->SetTitle("#it{y}");
    h_y17[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y17[i]->GetYaxis()->SetTitle("dN/d#it{y}");
    h_y17[i]->SetTitle("");
    if(i > 3) {
      h_y17[i]->Scale(h_y17[i-4]->Integral() / h_y17[i]->Integral());
      h_y17[i]->SetLineStyle(kDashed);
    }
    if(i == 0) h_y17[i]->Draw("error");
    else h_y17[i]->Draw("error same");  
  }

  TLegend *legcys = new TLegend(0.575, 0.6, 0.875, 0.725);
  legcys->SetTextSize(0.04);
  legcys->SetBorderSize(0);
  legcys->SetFillColorAlpha(kWhite,0);
  legcys->AddEntry(h_m17[1], "PRS Data", "pl");
  legcys->AddEntry(h_m17[5], "Scaled MC", "pl");
  legcys->Draw();


  TLatex lcys;
  lcys.SetTextSize(0.04);
  yl = getPos(h_y17[0]->GetMinimum(), h_y17[0]->GetMaximum(), 0.9, 1);
  lcys.DrawLatex(-1, yl, "#bf{2018 J/#psi}");
  lcys.SetTextColor(coly[0]);
  yl = getPos(h_y17[0]->GetMinimum(), h_y17[0]->GetMaximum(), 0.7, 1);
  lcys.DrawLatex(-1, yl, "#bf{[25,45] GeV}");
  lcys.SetTextColor(coly[1]);
  yl = getPos(h_y17[0]->GetMinimum(), h_y17[0]->GetMaximum(), 0.62, 1);
  lcys.DrawLatex(-1, yl, "#bf{[45,50] GeV}");
  lcys.SetTextColor(coly[2]);
  yl = getPos(h_y17[0]->GetMinimum(), h_y17[0]->GetMaximum(), 0.54, 1);
  lcys.DrawLatex(-1, yl, "#bf{[50,70] GeV}");
  lcys.SetTextColor(coly[3]);
  yl = getPos(h_y17[0]->GetMinimum(), h_y17[0]->GetMaximum(), 0.46, 1);
  lcys.DrawLatex(-1, yl, "#bf{[70,120] GeV}");

  c->SaveAs(Form("plots/y_scale.pdf"));
  c->Clear();

  // plot ratio data/MC
  c->SetLogy(0);
  TH1D **h_y17r = new TH1D*[4]; 
  for(int i = 0; i < 4; i++) {
    h_y17r[i] = (TH1D*)h_y17[i]->Clone(Form("rH_%d",i));
    h_y17r[i]->Sumw2();
    h_y17r[i]->Divide(h_y17[i+4]);

    h_y17r[i]->SetMinimum(0.51);
    h_y17r[i]->SetMaximum(1.49);
    h_y17r[i]->SetStats(0);
    h_y17r[i]->SetLineColor(coly[i]);
    h_y17r[i]->GetXaxis()->SetTitle("#it{y}");
    h_y17r[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y17r[i]->GetXaxis()->CenterTitle(true);
    h_y17r[i]->GetXaxis()->SetLabelOffset(0.01);
    h_y17r[i]->GetYaxis()->SetLabelOffset(0.01);
    h_y17r[i]->GetYaxis()->SetTitle("Data/MC");
    h_y17r[i]->SetTitle("");

    if(i == 0) h_y17r[i]->Draw("error");
    else h_y17r[i]->Draw("error same");  
  }

  TLatex lcyr;
  lcyr.SetTextSize(0.04);
  yl = getPos(h_y17r[0]->GetMinimum(), h_y17r[0]->GetMaximum(), 0.9, 0);
  lcyr.DrawLatex(-1.1, yl, "#bf{2018 J/#psi}");
  lcyr.SetTextColor(coly[0]);
  yl = getPos(h_y17r[0]->GetMinimum(), h_y17r[0]->GetMaximum(), 0.88, 0);
  lcyr.DrawLatex(0.35, yl, "#bf{[25,45] GeV}");
  lcyr.SetTextColor(coly[1]);
  yl = getPos(h_y17r[0]->GetMinimum(), h_y17r[0]->GetMaximum(), 0.8, 0);
  lcyr.DrawLatex(0.35, yl, "#bf{[45,50] GeV}");
  lcyr.SetTextColor(coly[2]);
  yl = getPos(h_y17r[0]->GetMinimum(), h_y17r[0]->GetMaximum(), 0.72, 0);
  lcyr.DrawLatex(0.35, yl, "#bf{[50,70] GeV}");
  lcyr.SetTextColor(coly[3]);
  yl = getPos(h_y17r[0]->GetMinimum(), h_y17r[0]->GetMaximum(), 0.64, 0);
  lcyr.DrawLatex(0.35, yl, "#bf{[70,120] GeV}");

  //  c->SaveAs(Form("plots/y_ratio.pdf"));
  c->Clear();

  // plot lifetime
  c->SetLogy();

  for(int i = 0; i < 4; i++) {
    h_lt17[i]->Scale(h_lt18[i]->Integral() / h_lt17[i]->Integral());
    h_lt17[i]->SetStats(0);
    h_lt17[i]->SetLineColor(coly[i]);
    h_lt17[i]->SetMarkerColor(coly[i]);
    h_lt17[i]->SetMarkerStyle(24);
    h_lt17[i]->SetMarkerSize(.5);
    h_lt17[i]->SetLineStyle(kDashed);
    h_lt17[i]->SetMinimum(4e2);
    h_lt17[i]->SetMaximum(3e6);
    h_lt17[i]->GetXaxis()->SetTitle("c#tau (#mum)");
    h_lt17[i]->GetXaxis()->SetTitleOffset(1.1);
    h_lt17[i]->GetXaxis()->CenterTitle(true);
    h_lt17[i]->GetXaxis()->SetLabelOffset(0.01);
    h_lt17[i]->GetYaxis()->SetLabelOffset(0.01);
    h_lt17[i]->GetYaxis()->SetTitle("dN/d c#tau");
    h_lt17[i]->SetTitle("");
    if(i == 0) h_lt17[i]->Draw("error");
    else h_lt17[i]->Draw("error same");  

    h_lt18[i]->SetLineColor(coly[i]);
    h_lt18[i]->SetMarkerColor(coly[i]);
    h_lt18[i]->SetMarkerStyle(20);
    h_lt18[i]->SetMarkerSize(.5);
    h_lt18[i]->Draw("error same");
  }

  TLegend *legclt = new TLegend(0.35, 0.725, 0.65, 0.85);
  legclt->SetTextSize(0.04);
  legclt->SetBorderSize(0);
  legclt->SetFillColorAlpha(kWhite,0);
  legclt->AddEntry(h_pT17[0], "2017", "pl");
  legclt->AddEntry(h_pT18[0], "2018", "pl");
  legclt->Draw();
  
  TLine *l_PRm = new TLine(-50, h_lt17[0]->GetMinimum(), -50, h_lt17[0]->GetMaximum());
  l_PRm->SetLineStyle(kDashed);
  l_PRm->SetLineColor(kBlack);
  l_PRm->Draw();
  TLine *l_PRp = new TLine(50, h_lt17[0]->GetMinimum(), 50, h_lt17[0]->GetMaximum());
  l_PRp->SetLineStyle(kDashed);
  l_PRp->SetLineColor(kBlack);
  l_PRp->Draw();

  TLine *l_NPm = new TLine(100, h_lt17[0]->GetMinimum(), 100, h_lt17[0]->GetMaximum());
  l_NPm->SetLineStyle(kDotted);
  l_NPm->SetLineColor(kBlack);
  l_NPm->Draw();
  TLine *l_NPp = new TLine(800, h_lt17[0]->GetMinimum(), 800, h_lt17[0]->GetMaximum());
  l_NPp->SetLineStyle(kDotted);
  l_NPp->SetLineColor(kBlack);
  l_NPp->Draw();

  TLatex lclt;
  lclt.SetTextSize(0.04);
  lclt.SetTextColor(colpt[0]);
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.9, 1);
  lclt.DrawLatex(150, yl, "#bf{J/#psi}");
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.95, 1);
  lclt.DrawLatex(400, yl, "#bf{Signal Region Data}");
  lclt.SetTextColor(colpt[1]);
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.88, 1);
  lclt.DrawLatex(500, yl, "#bf{[25,45] GeV}");
  lclt.SetTextColor(colpt[2]);
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.81, 1);
  lclt.DrawLatex(500, yl, "#bf{[45,50] GeV}");
  lclt.SetTextColor(colpt[3]);
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.74, 1);
  lclt.DrawLatex(500, yl, "#bf{[50,70] GeV}");
  lclt.SetTextColor(colpt[4]);
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.67, 1);
  lclt.DrawLatex(500, yl, "#bf{[70,120] GeV}");

  c->SaveAs(Form("plots/lt_all.pdf"));
  c->Clear();

  // lifetime scaling to low pT
  for(int i = 0; i < 4; i++) {
    h_lt17[i]->SetStats(0);
    h_lt17[i]->SetLineColor(coly[i]);
    h_lt17[i]->Scale(h_lt17[0]->Integral()/h_lt17[i]->Integral());
    h_lt17[i]->SetMinimum(4e2);
    h_lt17[i]->SetMaximum(3e6);
    h_lt17[i]->GetXaxis()->SetTitle("c#tau (#mum)");
    h_lt17[i]->GetXaxis()->SetTitleOffset(1.1);
    h_lt17[i]->GetYaxis()->SetTitle("dN/d c#tau (a.u.)");
    h_lt17[i]->SetTitle("");
    if(i == 0) h_lt17[i]->Draw("error");
    else h_lt17[i]->Draw("error same");  
  }

  l_PRm->Draw();
  l_PRp->Draw();
  l_NPm->Draw();
  l_NPp->Draw();

  TLatex lclts;
  lclts.SetTextSize(0.04);
  lclts.SetTextColor(colpt[0]);
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.9, 1);
  lclts.DrawLatex(150, yl, "#bf{2018 J/#psi}");
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.5, 1);
  lclts.DrawLatex(125, yl, "#bf{Signal Region Data}");
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.1, 1);
  double xl = getPos(-50, 50, 0.5, 0);
  lclts.SetTextAlign(21);
  lclts.DrawLatex(xl, yl, "#bf{PR}");
  xl = getPos(100, 800, 0.5, 0);
  lclts.DrawLatex(xl, yl, "#bf{NP}");
  lclts.SetTextColor(colpt[1]);
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.43, 1);
  lclts.SetTextAlign(11);
  lclts.DrawLatex(225, yl, "#bf{[25,45] GeV}");
  lclts.SetTextColor(colpt[2]);
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.36, 1);
  lclts.DrawLatex(225, yl, "#bf{[45,50] GeV}");
  lclts.SetTextColor(colpt[3]);
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.29, 1);
  lclts.DrawLatex(225, yl, "#bf{[50,70] GeV}");
  lclts.SetTextColor(colpt[4]);
  yl = getPos(h_lt17[0]->GetMinimum(), h_lt17[0]->GetMaximum(), 0.22, 1);
  lclts.DrawLatex(225, yl, "#bf{[70,120] GeV}");

  c->SaveAs(Form("plots/lt_scale.pdf"));
  c->Clear();

  // plot M
  c->SetLogy(0);
  c->SetLeftMargin(0.12);

  for(int i = 0; i < 8; i++) {
    h_m17[i]->SetStats(0);
    h_m17[i]->Scale(1./h_m17[i]->Integral());
    h_m17[i]->Scale(h_m17[i%4]->GetMaximum()/h_m17[i]->GetMaximum());
    h_m17[i]->SetLineColor(coly[i%4]);
    h_m17[i]->SetMarkerColor(coly[i%4]);
    h_m17[i]->SetMarkerStyle(20);
    h_m17[i]->SetMarkerSize(.5);
    h_m17[i]->GetXaxis()->SetTitle("#it{m} (GeV)");
    h_m17[i]->GetXaxis()->CenterTitle(true);
    h_m17[i]->GetXaxis()->SetTitleOffset(1.1);
    h_m17[i]->GetXaxis()->SetLabelOffset(0.01);
    h_m17[i]->GetYaxis()->SetLabelOffset(0.01);
    h_m17[i]->GetYaxis()->SetTitle("dN/d#it{m} (normalized)");
    h_m17[i]->SetTitle("");
    if(i > 3){
      h_m17[i]->SetLineStyle(kDashed);
      h_m17[i]->SetMarkerStyle(24);
    }
  }
  
  for(int i = 0; i < 8; i++) {  
    if(i%4==1 || i%4 == 2) continue;
    
    h_m17[i]->SetMinimum(0);
    h_m17[i]->SetMaximum(0.085);
    if(i==0) h_m17[i]->Draw("error");
    else h_m17[i]->Draw("error same");
  }

  TLatex lcm;
  lcm.SetTextSize(0.04);
  lcm.DrawLatex(2.95, 0.074, "#bf{2018 J/#psi}");
  lcm.SetTextColor(coly[0]);
  lcm.DrawLatex(3.16, 0.049, "#bf{#it{p}_{T}: 25-45 GeV}");
  lcm.SetTextColor(coly[3]);
  lcm.DrawLatex(3.16, 0.041, "#bf{#it{p}_{T}: 70-120 GeV}");

  // set a line black just for the legend
  h_m17[1]->SetLineColor(kBlack);
  h_m17[1]->SetMarkerColor(kBlack);
  h_m17[5]->SetLineColor(kBlack);
  h_m17[5]->SetMarkerColor(kBlack);
   
  TLegend *legcm = new TLegend(0.725, 0.77, 1.025, 0.92);
  legcm->SetTextSize(0.04);
  legcm->SetBorderSize(0);
  legcm->SetFillColorAlpha(kWhite,0);
  legcm->AddEntry(h_m17[1], "PR Data", "pl");
  legcm->AddEntry(h_m17[5], "MC", "pl");
  legcm->Draw();
  
  c->SaveAs(Form("plots/m_scale.pdf"));
  c->Clear();

  c->SetTopMargin(0.04);

  h_m17[8]->SetStats(0);
  h_m17[8]->SetLineColor(kBlack);
  h_m17[8]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
  h_m17[8]->GetXaxis()->SetTitleOffset(1.1);
  h_m17[8]->GetYaxis()->SetTitle("dN/dM");
  h_m17[8]->GetYaxis()->SetTitleOffset(1.8);
  h_m17[8]->SetMaximum(1.2e6);
  h_m17[8]->SetTitle("");
  h_m17[8]->Draw("error");  

  TLatex lcmf;
  lcmf.SetTextSize(0.04);
  lcmf.SetTextColor(colpt[0]);
  lcmf.DrawLatex(3.15, 1e6, "Run 2 Data");
  
  //  c->SaveAs(Form("plots/m_full.pdf"));
  c->Clear();

  c->Destructor();
  fin17->Close();
}
