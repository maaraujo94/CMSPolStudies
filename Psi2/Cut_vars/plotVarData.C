// code to plot the cut variables for the Data Psi(2S) sample
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y, ct/cterr
- vProb
*/

void plotVarData()
{
  TFile *infile = new TFile("../Store_data_codes/Psi2_data_cuts.root");
  TH1D *h_muPpT = (TH1D*)infile->Get("h_muPpT");
  TH1D *h_muNpT = (TH1D*)infile->Get("h_muNpT");
  TH1D *h_muPEta = (TH1D*)infile->Get("h_muPEta");
  TH1D *h_muNEta = (TH1D*)infile->Get("h_muNEta");
  TH1D *h_JMass = (TH1D*)infile->Get("h_JMass");
  TH1D *h_JPt = (TH1D*)infile->Get("h_JPt");
  TH1D *h_Jy = (TH1D*)infile->Get("h_Jy");
  TH1D *h_Jlts = (TH1D*)infile->Get("h_Jlts");
  TH1D *h_vP = (TH1D*)infile->Get("h_vP");
  
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLogy();

  h_muPpT->GetYaxis()->SetRangeUser(1e1, 2e6);
  h_muPpT->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
  h_muPpT->SetLineColor(kRed);
  h_muPpT->SetMarkerColor(kRed);
  h_muPpT->SetTitle("Inclusive data muon p_{T}");
  h_muPpT->Draw("error");
  h_muNpT->SetLineColor(kBlue);
  h_muNpT->SetMarkerColor(kBlue);
  h_muNpT->Draw("error same");
  TLine *mupT = new TLine(5.6, 1e1, 5.6, 2e6);
  mupT->SetLineStyle(kDashed);
  mupT->SetLineColor(kBlack);
  mupT->Draw("lsame");
  c->SaveAs("plots/data_muon_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  h_muPEta->GetYaxis()->SetRangeUser(0, 1.5e5);
  h_muPEta->GetXaxis()->SetTitle("#eta(#mu)");
  h_muPEta->SetLineColor(kRed);
  h_muPEta->SetMarkerColor(kRed);
  h_muPEta->SetTitle("Inclusive data muon #eta");
  h_muPEta->Draw("error");
  h_muNEta->SetLineColor(kBlue);
  h_muNEta->SetMarkerColor(kBlue);
  h_muNEta->Draw("error same");
  TLine *muEta_1 = new TLine(-1.4, 0, -1.4, 1.5e5);
  muEta_1->SetLineStyle(kDashed);
  muEta_1->SetLineColor(kBlack);
  muEta_1->Draw("lsame");
  TLine *muEta_2 = new TLine(1.4, 0, 1.4, 1.5e5);
  muEta_2->SetLineStyle(kDashed);
  muEta_2->SetLineColor(kBlack);
  muEta_2->Draw("lsame");
  c->SaveAs("plots/data_muon_eta.pdf");
  c->Clear();

  c->SetLogy(0);
  h_JMass->GetYaxis()->SetRangeUser(0, 3.5e5);
  h_JMass->GetXaxis()->SetTitle("M(#psi(2S)) (GeV)");
  h_JMass->SetTitle("Inclusive data #psi(2S) mass");
  h_JMass->Draw("error");
  c->SaveAs("plots/data_psi2_mass.pdf");
  c->Clear();

  c->SetLogy();
  //h_JPt->SetStats(0);
  h_JPt->GetYaxis()->SetRangeUser(2e1, 5e6);
  h_JPt->GetXaxis()->SetTitle("p_{T}(#psi(2S)) (GeV)");
  //h_JPt->GetYaxis()->SetTitle("dN/dp_{T}");
  h_JPt->SetTitle("Inclusive data #psi(2S) p_{T}");
  h_JPt->Draw("error");
  TLine *JPt_1 = new TLine(25, 2e1, 25, 5e6);
  JPt_1->SetLineStyle(kDashed);
  JPt_1->SetLineColor(kBlack);
  JPt_1->Draw("lsame");
  TLine *JPt_2 = new TLine(100, 2e1, 100, 5e6);
  JPt_2->SetLineStyle(kDashed);
  JPt_2->SetLineColor(kBlack);
  JPt_2->Draw("lsame");
  c->SaveAs("plots/data_psi2_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  //h_Jy->SetStats(0);
  h_Jy->GetYaxis()->SetRangeUser(0, 1.5e5);
  h_Jy->GetXaxis()->SetTitle("y(#psi(2S))");
  //h_Jy->GetYaxis()->SetTitle("dN/dy");
  h_Jy->SetTitle("Inclusive data #psi(2S) y");
  h_Jy->Draw("error");
  TLine *JRap_1 = new TLine(-1.2, 0, -1.2, 1.5e5);
  JRap_1->SetLineStyle(kDashed);
  JRap_1->SetLineColor(kBlack);
  JRap_1->Draw("lsame");
  TLine *JRap_2 = new TLine(1.2, 0, 1.2, 1.5e5);
  JRap_2->SetLineStyle(kDashed);
  JRap_2->SetLineColor(kBlack);
  JRap_2->Draw("lsame");
  c->SaveAs("plots/data_psi2_rap.pdf");
  c->Clear();

  c->SetLogy();
  h_Jlts->GetYaxis()->SetRangeUser(1e4, 1e6);
  h_Jlts->GetXaxis()->SetTitle("|c#tau|/#sigma_{c#tau}");
  h_Jlts->SetTitle("Inclusive data #psi(2S) lts");
  h_Jlts->Draw("error");
  TLine *Jlts_0 = new TLine(-2.5, 1e4, -2.5, 1e6);
  Jlts_0->SetLineStyle(kDashed);
  Jlts_0->SetLineColor(kBlack);
  Jlts_0->Draw("lsame");
  TLine *Jlts_1 = new TLine(2.5, 1e4, 2.5, 1e6);
  Jlts_1->SetLineStyle(kDashed);
  Jlts_1->SetLineColor(kBlack);
  Jlts_1->Draw("lsame");
  TLine *Jlts_2 = new TLine(4, 1e4, 4, 1e6);
  Jlts_2->SetLineStyle(kDashed);
  Jlts_2->SetLineColor(kBlack);
  Jlts_2->Draw("lsame");
  c->SaveAs("plots/data_psi2_lts.pdf");
  c->Clear();

  c->SetLogy(0);
  h_vP->GetYaxis()->SetRangeUser(1e4, 2e5);
  h_vP->GetXaxis()->SetTitle("vProb");
  h_vP->SetTitle("Inclusive data #psi(2S) vProb");
  h_vP->Draw("error");
  TLine *vP_0 = new TLine(0.01, 1e4, 0.01, 2e5);
  vP_0->SetLineStyle(kDashed);
  vP_0->SetLineColor(kBlack);
  vP_0->Draw("lsame");
  c->SaveAs("plots/data_psi2_vP.pdf");

  infile->Close();
  
}
