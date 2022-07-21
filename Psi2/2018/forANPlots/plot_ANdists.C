// macro to take the stored histos and plot nicely
void plot_ANdists()
{
  // PART 1 : reading the histograms

  // 3 pT dists: PRSR data + the 2 MC
  TH1D **h_pT = new TH1D*[3]; 
  // 4 y dists: PRSR data + MC over 2 pT regions
  TH1D **h_y = new TH1D*[4]; 
  // 5 M dists: PR data + MC over 2 pT regions + full pT data
  TH1D **h_m = new TH1D*[5]; 
  // 3 lifetime dists: data in 3 pT regions
  TH1D **h_lt = new TH1D*[3];
  // 9 costh dists: PR data + NP data + MC over 3 pT regions
  TH1D **h_cos = new TH1D*[9]; 
  
  TFile *fin = new TFile("files/store_ANdists.root");
 
  // read the pT dists
  string lbl_pt[] = {"Data", "lowPtMC", "highPtMC"};
  for(int i = 0; i < 3; i++) {
    h_pT[i] = (TH1D*)fin->Get(Form("h_pT_%s", lbl_pt[i].c_str()));
  }

  // store the y dists
  string lbl_y[] = {"lowPtData", "highPtData", "lowPtMC", "highPtMC"};
  for(int i = 0; i < 4; i++) {
    h_y[i] = (TH1D*)fin->Get(Form("h_y_%s", lbl_y[i].c_str()));
  }

  // store the M dists
  string lbl_m[] = {"lowPtData", "highPtData", "lowPtMC", "highPtMC", "fullData"};
  for(int i = 0; i < 5; i++) {
    h_m[i] = (TH1D*)fin->Get(Form("h_m_%s", lbl_m[i].c_str()));
  }

  // store the lifetime dists
  string lbl_lt[] = {"lowPt", "highPt", "full"};
  for(int i = 0; i < 3; i++) {
    h_lt[i] = (TH1D*)fin->Get(Form("h_lt_%s", lbl_lt[i].c_str()));
  }

  // store the cos dists
  string lbl_cos[] = {"fullPR", "lowPtPR", "highPtPR", "fullNP", "lowPtNP", "highPtNP", "fullMC", "lowPtMC", "highPtMC"};
  for(int i = 0; i < 9; i++) {
    h_cos[i] = (TH1D*)fin->Get(Form("h_cos_%s", lbl_cos[i].c_str()));
  }

  TCanvas *c = new TCanvas("", "", 900, 900);

  // plot pT
  c->SetLogy();
  c->SetLeftMargin(0.11);
  int colpt[] = {kBlack, kRed+1, kBlue};
  double ptlim[] = {46, 300};

  for(int i = 0; i < 3; i++) {

    h_pT[i]->SetStats(0);
    h_pT[i]->SetLineColor(colpt[i]);
    h_pT[i]->SetMinimum(1e2);
    h_pT[i]->SetMaximum(1e7);
    h_pT[i]->GetXaxis()->SetTitle("p_{T}(#mu#mu) (GeV)");
    h_pT[i]->GetXaxis()->SetTitleOffset(1.1);
    h_pT[i]->GetYaxis()->SetTitle("dN/dp_{T}");
    h_pT[i]->SetTitle("");
    if(i == 0) h_pT[i]->Draw("histo");
    else {
      // set bins above pt_max to zero
      int nb = h_pT[i]->GetNbinsX();
      for(int ib = 0; ib < nb; ib++) {
	double pt = h_pT[i]->GetXaxis()->GetBinUpEdge(ib+1);
	if (pt > ptlim[i-1])
	  h_pT[i]->SetBinContent(ib+1, 0);
      }
    
    h_pT[i]->Draw("histo same");
    }
  }
  
  TLatex lcpt;
  lcpt.SetTextSize(0.04);
  lcpt.SetTextColor(colpt[0]);
  lcpt.DrawLatex(120, 3.5e6, "Peak data");
  lcpt.SetTextColor(colpt[1]);
  lcpt.DrawLatex(120, 1.5e6, "MC low p_{T}");
  lcpt.SetTextColor(colpt[2]);
  lcpt.DrawLatex(120, 6e5, "MC high p_{T}");

  TLine *ptL = new TLine(120, 0, 120,  exp(0.5*(log(h_pT[0]->GetMaximum())+log(h_pT[0]->GetMinimum()))));
  ptL->SetLineStyle(kDashed);
  ptL->Draw();

  c->SaveAs(Form("plots/ANdists/pt_all.pdf"));
  c->Clear();
  
  // plot y
  c->SetLogy();
  int coly[] = {kRed+1, kBlue};
  for(int i = 0; i < 4; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(coly[i%2]);
    h_y[i]->SetMinimum(5e2);
    h_y[i]->SetMaximum(7e4);
    h_y[i]->GetXaxis()->SetTitle("y(#mu#mu)");
    h_y[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y[i]->GetYaxis()->SetTitle("dN/dy");
    h_y[i]->SetTitle("");
    if(i > 1) h_y[i]->SetLineStyle(kDashed);
    if(i == 0) h_y[i]->Draw("histo");
    else h_y[i]->Draw("histo same");  
  }

  TLatex lcy;
  lcy.SetTextSize(0.025);
  lcy.SetTextColor(coly[0]);
  lcy.DrawLatex(0.875, 1.5e4, "Low-p_{T}");
  lcy.SetTextColor(coly[1]);
  lcy.DrawLatex(0.875, 2e3, "High-p_{T}");
  
  c->SaveAs(Form("plots/ANdists/y_all.pdf"));
  c->Clear();
  
  // y scaling MC to data
  for(int i = 0; i < 4; i++) {

    h_y[i]->SetStats(0);
    h_y[i]->SetLineColor(coly[i%2]);
    h_y[i]->SetMinimum(5e2);
    h_y[i]->SetMaximum(7e4);
    h_y[i]->GetXaxis()->SetTitle("y(#mu#mu)");
    h_y[i]->GetXaxis()->SetTitleOffset(1.1);
    h_y[i]->GetYaxis()->SetTitle("dN/dy");
    h_y[i]->SetTitle("");
    if(i > 1) {
      h_y[i]->Scale(h_y[i-2]->Integral() / h_y[i]->Integral());
      h_y[i]->SetLineStyle(kDashed);
    }
    if(i == 0) h_y[i]->Draw("histo");
    else h_y[i]->Draw("histo same");  
  }

  TLatex lcys;
  lcys.SetTextSize(0.025);
  lcys.SetTextColor(coly[0]);
  lcys.DrawLatex(0.875, 1.5e4, "Low-p_{T}");
  lcys.SetTextColor(coly[1]);
  lcys.DrawLatex(0.875, 2e3, "High-p_{T}");

  c->SaveAs(Form("plots/ANdists/y_scale.pdf"));
  c->Clear();
  
  // plot lifetime
  c->SetLogy();

  for(int i = 0; i < 2; i++) {

    h_lt[i]->SetStats(0);
    h_lt[i]->SetLineColor(coly[i]);
    h_lt[i]->SetMinimum(5e2);
    h_lt[i]->SetMaximum(3e5);
    h_lt[i]->GetXaxis()->SetTitle("c#tau (#mum)");
    h_lt[i]->GetXaxis()->SetTitleOffset(1.1);
    h_lt[i]->GetYaxis()->SetTitle("dN/dc#tau");
    h_lt[i]->SetTitle("");
    if(i == 0) h_lt[i]->Draw("histo");
    else h_lt[i]->Draw("histo same");  
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
  TLine *l_NPp = new TLine(500, h_lt[0]->GetMinimum(), 500, h_lt[0]->GetMaximum());
  l_NPp->SetLineStyle(kDotted);
  l_NPp->SetLineColor(kBlack);
  l_NPp->Draw();

  TLatex lclt;
  lclt.SetTextSize(0.04);
  lclt.SetTextColor(colpt[0]);
  lclt.DrawLatex(300, 2e5, "Data:");
  lclt.SetTextColor(colpt[1]);
  lclt.DrawLatex(300, 1.2e5, "low p_{T}");
  lclt.SetTextColor(colpt[2]);
  lclt.DrawLatex(300, 7e4, "high p_{T}");

  c->SaveAs(Form("plots/ANdists/lt_all.pdf"));
  c->Clear();
  
  // plot just full lt
  h_lt[2]->SetStats(0);
  h_lt[2]->SetLineColor(colpt[0]);
  h_lt[2]->SetMinimum(5e2);
  h_lt[2]->SetMaximum(3e5);
  h_lt[2]->GetXaxis()->SetTitle("c#tau (#mum)");
  h_lt[2]->GetXaxis()->SetTitleOffset(1.1);
  h_lt[2]->GetYaxis()->SetTitle("dN/dc#tau");
  h_lt[2]->SetTitle("");
  h_lt[2]->Draw("histo");
  
  l_PRm->Draw();
  l_PRp->Draw();
  l_NPm->Draw();
  l_NPp->Draw();

  TLatex lcltf;
  lcltf.SetTextSize(0.04);
  lcltf.SetTextColor(colpt[0]);
  lcltf.DrawLatex(300, 2e5, "2018 Data");

  c->SaveAs(Form("plots/ANdists/lt_full.pdf"));
  c->Clear();
  
  // lifetime scaling to low pT
  for(int i = 0; i < 2; i++) {
    h_lt[i]->SetStats(0);
    h_lt[i]->SetLineColor(coly[i]);
    h_lt[i]->Scale(h_lt[0]->Integral()/h_lt[i]->Integral());
    h_lt[i]->SetMinimum(5e2);
    h_lt[i]->SetMaximum(3e5);
    h_lt[i]->GetXaxis()->SetTitle("c#tau (#mum)");
    h_lt[i]->GetXaxis()->SetTitleOffset(1.1);
    h_lt[i]->GetYaxis()->SetTitle("dN/dc#tau (a.u.)");
    h_lt[i]->SetTitle("");
    if(i == 0) h_lt[i]->Draw("histo");
    else h_lt[i]->Draw("histo same");  
  }

  l_PRm->Draw();
  l_PRp->Draw();
  l_NPm->Draw();
  l_NPp->Draw();

  TLatex lclts;
  lclts.SetTextSize(0.04);
  lclts.SetTextColor(colpt[0]);
  lclts.DrawLatex(300, 2e5, "Data:");
  lclts.DrawLatex(-40, 1e3, "Peak");
  lclts.DrawLatex(275, 1e3, "NP");
  
  lclts.SetTextColor(colpt[1]);
  lclts.DrawLatex(300, 1.2e5, "low p_{T}");
  lclts.SetTextColor(colpt[2]);
  lclts.DrawLatex(300, 7e4, "high p_{T}");
  
  c->SaveAs(Form("plots/ANdists/lt_scale.pdf"));
  c->Clear();
  
  // plot M
  c->SetLogy(0);
  c->SetLeftMargin(0.12);

  for(int i = 0; i < 4; i++) {
    h_m[i]->SetStats(0);
    h_m[i]->Scale(1./h_m[i]->Integral());
    h_m[i]->Scale(h_m[i%2]->GetMaximum()/h_m[i]->GetMaximum());
    h_m[i]->SetLineColor(coly[i%2]);
    h_m[i]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
    h_m[i]->GetXaxis()->SetTitleOffset(1.1);
    h_m[i]->GetYaxis()->SetTitle("dN/dM (normalized)");
    h_m[i]->SetTitle("");
    if(i > 1) h_m[i]->SetLineStyle(kDashed);
  }

  for(int i = 0; i < 4; i++) {   
    h_m[i]->SetMinimum(0);
    h_m[i]->SetMaximum(0.035);
    if(i==0) h_m[i]->Draw("histo");
    else h_m[i]->Draw("histo same");
  }

  TLatex lcm;
  lcm.SetTextSize(0.04);
  lcm.SetTextColor(coly[0]);
  lcm.DrawLatex(3.8, 0.03, "low-p_{T}");
  lcm.SetTextColor(coly[1]);
  lcm.DrawLatex(3.8, 0.027, "high-p_{T}");
  lcm.SetTextColor(kBlack);
  
  c->SaveAs(Form("plots/ANdists/m_scale.pdf"));
  c->Clear();
  
  h_m[4]->SetStats(0);
  h_m[4]->SetLineColor(kBlack);
  h_m[4]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
  h_m[4]->GetXaxis()->SetTitleOffset(1.1);
  h_m[4]->GetYaxis()->SetTitle("dN/dM");
  h_m[4]->SetMaximum(3.3e4);
  h_m[4]->SetTitle("");
  h_m[4]->Draw("histo");  

  TLatex lcmf;
  lcmf.SetTextSize(0.04);
  lcmf.SetTextColor(colpt[0]);
  lcmf.DrawLatex(3.8, 3e4, "2018 Data");

  c->SaveAs(Form("plots/ANdists/m_full.pdf"));
  c->Clear();
  
  // plot costh
  c->SetLogy(0);

  for(int i = 7; i < 9; i++) {

    h_cos[i]->SetStats(0);
    h_cos[i]->Scale(0.06/h_cos[i]->GetBinContent(1));
    h_cos[i]->SetLineColor(coly[i-7]);
    h_cos[i]->SetMinimum(0);
    h_cos[i]->SetMaximum(0.09);
    h_cos[i]->GetXaxis()->SetTitle("|cos#theta|");
    h_cos[i]->GetXaxis()->SetTitleOffset(1.);
    h_cos[i]->GetYaxis()->SetTitle("dN/d|cos#theta| (a.u.)");
    if(i == 7) h_cos[i]->Draw("histo");
    else h_cos[i]->Draw("histo same"); 
  }

  TLatex lcc;
  lcc.SetTextSize(0.04);
  lcc.SetTextColor(coly[0]);
  lcc.DrawLatex(0.05, 0.03, "MC low p_{T}");
  lcc.SetTextColor(coly[1]);
  lcc.DrawLatex(0.05, 0.023, "MC high p_{T}");
  
  c->SaveAs(Form("plots/ANdists/cos_MC.pdf"));
  c->Clear();
  
  int colc[] = {kViolet-1, kRed, kBlack};

  for(int i = 0; i < 3; i++) {

    int j = i*3;
    
    h_cos[j]->SetStats(0);
    h_cos[j]->Scale(0.075/h_cos[j]->GetBinContent(1));
    h_cos[j]->SetLineColor(colc[i]);
    h_cos[j]->SetMinimum(0);
    h_cos[j]->SetMaximum(0.09);
    h_cos[j]->GetXaxis()->SetTitle("|cos#theta|");
    h_cos[j]->GetXaxis()->SetTitleOffset(1.);
    h_cos[j]->GetYaxis()->SetTitle("dN/d|cos#theta| (a.u.)");
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
  
  fin->Close();
}
