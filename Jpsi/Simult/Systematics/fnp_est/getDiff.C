// get the ratio between old and new fNP_psi for psi(2S) 
void getDiff()
{
  // get the two fNP_psi for psi(2S)
  TFile *fin = new TFile("files/NPFrac_comp.root");
  TH1D *h_oldF = (TH1D*)fin->Get("h_fnpc_b");
  h_oldF->SetDirectory(0);
  TH1D *h_newF = (TH1D*)fin->Get("h_fnp_psi");
  h_newF->SetDirectory(0);
  fin->Close();

  const int n_pt = h_oldF->GetNbinsX();

  // calculate ratio (factor to apply to J/psi fNP_psi to get estimate of value wih "new" method)
  double r_fNP[n_pt], pt[n_pt];
  double er_fNP[n_pt], ept[n_pt];
  for(int i = 0; i < n_pt; i++) {
    pt[i] = 0.5*(h_newF->GetXaxis()->GetBinUpEdge(i+1)+h_newF->GetXaxis()->GetBinLowEdge(i+1));
    ept[i] = 0.5*(h_newF->GetXaxis()->GetBinUpEdge(i+1)-h_newF->GetXaxis()->GetBinLowEdge(i+1));
    r_fNP[i] = h_newF->GetBinContent(i+1)/h_oldF->GetBinContent(i+1);
    er_fNP[i] = 0;

    cout << i << " " << r_fNP[i] << endl;
  }
  TGraphErrors *g_rfNP = new TGraphErrors(n_pt, pt, r_fNP, ept, er_fNP);
  double r_f_max =TMath::MaxElement(n_pt, r_fNP);
  cout << r_f_max << endl;
  
  // get the J/psi fNP_psi from "old" method and estimate "new" value
  TFile *fin_J = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/PR_fit/files/NPFrac.root");
  TH2D *h_fnp_J = (TH2D*)fin_J->Get("h_fNPc");
  h_fnp_J->SetDirectory(0);
  fin_J->Close();

  // new value
  TH1D *h_fnp_J_1d = h_fnp_J->ProjectionY("h_fNP_J", 1, 1);
  const int n_ptJ = h_fnp_J_1d->GetNbinsX();
  h_fnp_J_1d->Scale(100);
  double fn_v[n_ptJ], ptJ[n_ptJ];
  double efn_v[n_ptJ], eptJ[n_ptJ];
  for(int i = 0; i < n_ptJ; i++) {
    ptJ[i] = 0.5*(h_fnp_J_1d->GetXaxis()->GetBinUpEdge(i+1)+h_fnp_J_1d->GetXaxis()->GetBinLowEdge(i+1));
    eptJ[i] = 0.5*(h_fnp_J_1d->GetXaxis()->GetBinUpEdge(i+1)-h_fnp_J_1d->GetXaxis()->GetBinLowEdge(i+1));

    fn_v[i] = h_fnp_J_1d->GetBinContent(i+1)*r_f_max;
    efn_v[i] = h_fnp_J_1d->GetBinError(i+1)*r_f_max;
  }
  TGraphErrors *g_fnpJ = new TGraphErrors(n_ptJ, ptJ, fn_v, eptJ, efn_v);
  
  // plot the ratio
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.11);
  c->SetTopMargin(0.015);

  double ptmin_P = h_newF->GetXaxis()->GetBinLowEdge(1);
  double ptmax_P = h_newF->GetXaxis()->GetBinUpEdge(n_pt);

  TH1F *fr1 = c->DrawFrame(ptmin_P-5, 0, ptmax_P+5, 1);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f_{NP} ratio");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("");

  g_rfNP->SetMarkerStyle(20);
  g_rfNP->SetMarkerColor(kBlack);
  g_rfNP->SetLineColor(kBlack);
  g_rfNP->Draw("psame");

  c->SaveAs("plots/ratioF.pdf");
  c->Clear();

  // plot the new fn
  double ptmin_J = h_fnp_J_1d->GetXaxis()->GetBinLowEdge(1);
  double ptmax_J = h_fnp_J_1d->GetXaxis()->GetBinUpEdge(n_ptJ);

  TH1F *ff1 = c->DrawFrame(ptmin_J-5, 0, ptmax_J+5, 40);
  ff1->SetXTitle("p_{T} (GeV)");
  ff1->SetYTitle("f_{NP}^{#psi} (%)");
  ff1->GetYaxis()->SetTitleOffset(1.3);
  ff1->GetYaxis()->SetLabelOffset(0.01);
  ff1->SetTitle("");

  g_fnpJ->SetMarkerStyle(20);
  g_fnpJ->SetMarkerSize(.5);
  g_fnpJ->SetMarkerColor(kBlack);
  g_fnpJ->SetLineColor(kBlack);
  g_fnpJ->Draw("psame");

  h_fnp_J_1d->SetMarkerStyle(25);
  h_fnp_J_1d->SetMarkerSize(.5);
  h_fnp_J_1d->SetMarkerColor(kBlack);
  h_fnp_J_1d->SetLineColor(kBlack);  
  h_fnp_J_1d->SetLineStyle(kDashed);
  h_fnp_J_1d->Draw("same");

  TLegend *leg2 = new TLegend(0.65, 0.85, 0.95, 0.95);
  leg2->SetTextSize(0.03);
  leg2->SetBorderSize(0);
  leg2->SetFillColorAlpha(kWhite,0);
  leg2->AddEntry(h_fnp_J_1d, "previous f_{NP}^{#psi}", "pl");
  leg2->AddEntry(g_fnpJ, "new f_{NP}^{#psi} (estimate)", "pl");
  leg2->Draw();

  c->SaveAs("plots/compF.pdf");
  c->Clear();
  c->Destructor();

  // store the new fNP to run the framework with it
  // has to be a th2d
  h_fnp_J->SetName("h_fNPc_old");
  TH2D *h_fNPc = (TH2D*)h_fnp_J->Clone("h_fNPc");
  h_fNPc->Scale(r_f_max);

  TFile *fout = new TFile("files/NPFrac.root", "recreate");
 h_fNPc->Write();
  fout->Close();
}
