// code to plot the cut variables for the Data Jpsi sample
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y, ct/cterr
*/

void plotVarData()
{
  TFile *infile = new TFile("../Store_data_codes/Jpsi_data_cuts.root");
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
  double norm;

  norm = h_muPpT->Integral();
  h_muPpT->Scale(1./norm);
  h_muNpT->Scale(1./norm);
  h_muPpT->GetYaxis()->SetRangeUser(1e-6, 1e-1);
  h_muPpT->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
  h_muPpT->SetLineColor(kRed);
  h_muPpT->SetTitle("Inclusive data muon p_{T}");
  h_muPpT->Draw("hist");
  h_muNpT->SetLineColor(kBlue);
  h_muNpT->Draw("hist same");
  TLine *mupT = new TLine(5.6, 1e-6, 5.6, 1e-1);
  mupT->SetLineStyle(kDashed);
  mupT->SetLineColor(kBlack);
  mupT->Draw("lsame");
  c->SaveAs("plots/data_muon_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_muPEta->Integral();
  h_muPEta->Scale(1./norm);
  h_muNEta->Scale(1./norm);
  h_muPEta->GetYaxis()->SetRangeUser(0, 0.035);
  h_muPEta->GetXaxis()->SetTitle("#eta(#mu)");
  h_muPEta->SetLineColor(kRed);
  h_muPEta->SetTitle("Inclusive data muon #eta");
  h_muPEta->Draw("hist");
  h_muNEta->SetLineColor(kBlue);
  h_muNEta->Draw("hist same");
  TLine *muEta_1 = new TLine(-1.4, 0, -1.4, 0.035);
  muEta_1->SetLineStyle(kDashed);
  muEta_1->SetLineColor(kBlack);
  muEta_1->Draw("lsame");
  TLine *muEta_2 = new TLine(1.4, 0, 1.4, 0.035);
  muEta_2->SetLineStyle(kDashed);
  muEta_2->SetLineColor(kBlack);
  muEta_2->Draw("lsame");
  c->SaveAs("plots/data_muon_eta.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_JMass->Integral();
  h_JMass->Scale(1./norm);
  h_JMass->GetYaxis()->SetRangeUser(0, 0.07);
  h_JMass->GetXaxis()->SetTitle("M(J/#psi) (GeV)");
  h_JMass->SetTitle("Inclusive data J/#psi mass");
  h_JMass->Draw("hist");
  /*TLine *JMass_1 = new TLine(3, 0, 3, 0.07);
  JMass_1->SetLineStyle(kDashed);
  JMass_1->SetLineColor(kBlack);
  JMass_1->Draw("lsame");
  TLine *JMass_2 = new TLine(3.2, 0, 3.2, 0.07);
  JMass_2->SetLineStyle(kDashed);
  JMass_2->SetLineColor(kBlack);
  JMass_2->Draw("lsame");*/
  c->SaveAs("plots/data_jpsi_mass.pdf");
  c->Clear();

  c->SetLogy();
  norm = h_JPt->Integral();
  //h_JPt->Scale(1./norm);
  h_JPt->SetStats(0);
  h_JPt->GetYaxis()->SetRangeUser(1e2, 1e7);
  h_JPt->GetXaxis()->SetTitle("p_{T}(J/#psi) (GeV)");
  h_JPt->GetYaxis()->SetTitle("dN/dp_{T}");
  h_JPt->SetTitle("Inclusive data J/#psi p_{T}");
  h_JPt->Draw("hist");
  TLine *JPt_1 = new TLine(25, 1e2, 25, 1e7);
  JPt_1->SetLineStyle(kDashed);
  JPt_1->SetLineColor(kBlack);
  JPt_1->Draw("lsame");
  TLine *JPt_2 = new TLine(100, 1e2, 100, 1e7);
  JPt_2->SetLineStyle(kDashed);
  JPt_2->SetLineColor(kBlack);
  JPt_2->Draw("lsame");
  c->SaveAs("plots/data_jpsi_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_Jy->Integral();
  //h_Jy->Scale(1./norm);
  h_Jy->SetStats(0);
  h_Jy->GetYaxis()->SetRangeUser(0, 6e5);
  h_Jy->GetXaxis()->SetTitle("y(J/#psi)");
  h_Jy->GetYaxis()->SetTitle("dN/dy");
  h_Jy->SetTitle("Inclusive data J/#psi y");
  h_Jy->Draw("hist");
  TLine *JRap_1 = new TLine(-1.2, 0, -1.2, 6e5);
  JRap_1->SetLineStyle(kDashed);
  JRap_1->SetLineColor(kBlack);
  JRap_1->Draw("lsame");
  TLine *JRap_2 = new TLine(1.2, 0, 1.2, 6e5);
  JRap_2->SetLineStyle(kDashed);
  JRap_2->SetLineColor(kBlack);
  JRap_2->Draw("lsame");
  c->SaveAs("plots/data_jpsi_rap.pdf");
  c->Clear();

  c->SetLogy();
  norm = h_Jlts->Integral();
  h_Jlts->Scale(1./norm);
  h_Jlts->GetYaxis()->SetRangeUser(1e-3, 3e-1);
  h_Jlts->GetXaxis()->SetTitle("|c#tau|/#sigma_{c#tau}");
  h_Jlts->SetTitle("Inclusive data J/#psi lts");
  h_Jlts->Draw("hist");
  TLine *Jlts_0 = new TLine(-2.5, 1e-3, -2.5, 3e-1);
  Jlts_0->SetLineStyle(kDashed);
  Jlts_0->SetLineColor(kBlack);
  Jlts_0->Draw("lsame");
  TLine *Jlts_1 = new TLine(2.5, 1e-3, 2.5, 3e-1);
  Jlts_1->SetLineStyle(kDashed);
  Jlts_1->SetLineColor(kBlack);
  Jlts_1->Draw("lsame");
  TLine *Jlts_2 = new TLine(4, 1e-3, 4, 3e-1);
  Jlts_2->SetLineStyle(kDashed);
  Jlts_2->SetLineColor(kBlack);
  Jlts_2->Draw("lsame");
  c->SaveAs("plots/data_jpsi_lts.pdf");
  c->Clear();

  c->SetLogy(0);
  h_vP->GetYaxis()->SetRangeUser(1.5e5, 4e5);
  h_vP->GetXaxis()->SetTitle("vProb");
  h_vP->SetTitle("Inclusive data J/#psi vProb");
  h_vP->Draw("error");
  TLine *vP_0 = new TLine(0.01, 1.5e5, 0.01, 4e5);
  vP_0->SetLineStyle(kDashed);
  vP_0->SetLineColor(kBlack);
  vP_0->Draw("lsame");
  c->SaveAs("plots/data_jpsi_vP.pdf");

  infile->Close();
  
}
