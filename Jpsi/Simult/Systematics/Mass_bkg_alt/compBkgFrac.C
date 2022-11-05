// macro to plot the bkg fracs
// f_bkg, f_bkg^NP, f_NP, f_NP^c
void compBkgFrac()
{
  int cols[] = {kBlack, kRed};
  
  // get fit fBG
  TH1D **fSB_b = new TH1D*[2];

  string basel = "/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult";
  string loc[2] = {"/bkgFits", "/Systematics/Mass_bkg_alt"};
  for(int i = 0; i < 2; i++) {
    TFile *fin1 = new TFile(Form("%s%s/files/bkgFrac.root", basel.c_str(), loc[i].c_str()));
    fSB_b[i] = (TH1D*)fin1->Get("fbkg_unc");
    fSB_b[i]->SetDirectory(0);
    fin1->Close();

    fSB_b[i]->Scale(100.);
  }

  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1F *fr1 = c->DrawFrame(20, 0.0, 125, 14);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f (%)");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("f_{bkg} comparison");

  string lbl[2] = {"Base", "Alt"};
  
  for(int i = 0; i < 2; i++) {
    fSB_b[i]->SetMarkerColor(cols[i]);
    fSB_b[i]->SetLineColor(cols[i]);
    fSB_b[i]->SetMarkerStyle(20);
    fSB_b[i]->SetMarkerSize(.5);
    fSB_b[i]->Draw("error same");
  }
  
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  for(int i = 0; i < 2; i++) 
    leg->AddEntry(fSB_b[i], Form("%s", lbl[i].c_str()), "pl");
  leg->Draw();
  
  c->SaveAs("plots/fBG_comp.pdf");
  c->Clear();

  // get fit fBG^NP
  TH1D **fSBN_b = new TH1D*[2];

  for(int i = 0; i < 2; i++) {
    TFile *fin1 = new TFile(Form("%s%s/files/bkgFrac_NP.root", basel.c_str(), loc[i].c_str()));
    fSBN_b[i] = (TH1D*)fin1->Get("fbkg_unc");
    fSBN_b[i]->SetDirectory(0);
    fin1->Close();

    //    cout << "dying here?" << endl;
    fSBN_b[i]->Scale(100.);
    //cout << "check" << endl;
  }

  TH1F *fr1n = c->DrawFrame(20, 0.0, 125, 14);
  fr1n->SetXTitle("p_{T} (GeV)");
  fr1n->SetYTitle("f (%)");
  fr1n->GetYaxis()->SetTitleOffset(1.3);
  fr1n->GetYaxis()->SetLabelOffset(0.01);
  fr1n->SetTitle("f_{bkg}^{NP} comparison");

  for(int i = 0; i < 2; i++) {
    fSBN_b[i]->SetMarkerColor(cols[i]);
    fSBN_b[i]->SetLineColor(cols[i]);
    fSBN_b[i]->SetMarkerStyle(20);
    fSBN_b[i]->SetMarkerSize(.5);
    fSBN_b[i]->Draw("error same");
  }
  
  leg->Draw();
  
  c->SaveAs("plots/fBGN_comp.pdf");
  c->Clear();
  c->Destructor();
}
