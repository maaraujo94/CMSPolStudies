// macro to plot the bkg fracs
// f_bkg, f_bkg^NP, f_NP
void compBkgFrac()
{
  // get fit fBG
  TH1D **fSB_b = new TH1D*[3];
  TH1D **h_fSB = new TH1D*[3];

  string basel = "/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/";
  string loc[3] = {"Simult", "Simult_dR2", "Simult_dR1"};
  for(int i = 0; i < 3; i++) {
    TFile *fin1 = new TFile(Form("%s%s/bkgFits/files/bkgFrac.root", basel.c_str(), loc[i].c_str()));
    fSB_b[i] = (TH1D*)fin1->Get("fbkg_unc");
    fSB_b[i]->SetDirectory(0);
    fin1->Close();

    fSB_b[i]->Scale(100.);
  }

  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1F *fr1 = c->DrawFrame(15, 0.0, 105, 100);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f (%)");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("Run 2 f_{bkg} comparison");

  string lbl[3] = {"no cut", "#DeltaR>0.15", "#DeltaR>0.17"};
  
  for(int i = 0; i < 3; i++) {
    fSB_b[i]->SetMarkerColor(i+1);
    fSB_b[i]->SetLineColor(i+1);
    fSB_b[i]->SetMarkerStyle(20);
    fSB_b[i]->SetMarkerSize(.5);
    fSB_b[i]->Draw("error same");
  }
  
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  for(int i = 0; i < 3; i++) 
    leg->AddEntry(fSB_b[i], Form("%s", lbl[i].c_str()), "pl");
  leg->Draw();
  
  c->SaveAs("plots/fBG_comp.pdf");
  c->Clear();

  // get fit fBG^NP
  TH1D **fSBN_b = new TH1D*[3];

  for(int i = 0; i < 3; i++) {
    TFile *fin1 = new TFile(Form("%s%s/bkgFits/files/bkgFrac_NP.root", basel.c_str(), loc[i].c_str()));
    fSBN_b[i] = (TH1D*)fin1->Get("fbkg_unc"); // wide-pT f_bkg
    fSBN_b[i]->SetDirectory(0);
    fin1->Close();

    fSBN_b[i]->Scale(100.);
  }

  TH1F *fr1n = c->DrawFrame(15, 0.0, 105, 50);
  fr1n->SetXTitle("p_{T} (GeV)");
  fr1n->SetYTitle("f (%)");
  fr1n->GetYaxis()->SetTitleOffset(1.3);
  fr1n->GetYaxis()->SetLabelOffset(0.01);
  fr1n->SetTitle("Run 2 f_{bkg}^{NP} comparison");

  for(int i = 0; i < 3; i++) {
    fSBN_b[i]->SetMarkerColor(i+1);
    fSBN_b[i]->SetLineColor(i+1);
    fSBN_b[i]->SetMarkerStyle(20);
    fSBN_b[i]->SetMarkerSize(.5);
    fSBN_b[i]->Draw("error same");
  }
  
  leg->Draw();
  
  c->SaveAs("plots/fBGN_comp.pdf");
  c->Clear();

  // get fit fNP
  TH1D **h_fNP = new TH1D*[3];

  for(int i = 0; i < 3; i++) {
    TFile *fin1 = new TFile(Form("%s%s/PR_fit/files/NPFrac.root", basel.c_str(), loc[i].c_str()));
    TH2D* fNP2d = (TH2D*)fin1->Get("h_fnp"); 
    fNP2d->SetDirectory(0);
    fin1->Close();

    h_fNP[i] = fNP2d->ProjectionY(Form("fnp_1d_%d", i), 1, 1);
    h_fNP[i]->Scale(100.);
  }

  TH1F *fr2 = c->DrawFrame(15, 0.0, 105, 40);
  fr2->SetXTitle("p_{T} (GeV)");
  fr2->SetYTitle("f (%)");
  fr2->GetYaxis()->SetTitleOffset(1.3);
  fr2->GetYaxis()->SetLabelOffset(0.01);
  fr2->SetTitle("Run 2 f_{NP} comparison");

  for(int i = 0; i < 3; i++) {
    h_fNP[i]->SetLineColor(i+1);
    h_fNP[i]->SetMarkerColor(i+1);
    h_fNP[i]->SetMarkerStyle(20);
    h_fNP[i]->SetMarkerSize(.5);
    h_fNP[i]->Draw("error same");
  }
  
  leg->Draw();
  
  c->SaveAs("plots/fNP_comp.pdf");
  c->Clear();
  c->Destructor();
}
