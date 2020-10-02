// code to plot the cut variables for the B->J/psi K sample
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y, ct/cterr
*/

void plotVarBJK()
{
  TFile *infile = new TFile("../../Store_data_codes/2018/BtoJK_data_cuts.root");
  TH1D *h_muPpT = (TH1D*)infile->Get("h_muPpT");
  TH1D *h_muNpT = (TH1D*)infile->Get("h_muNpT");
  TH1D *h_muPEta = (TH1D*)infile->Get("h_muPEta");
  TH1D *h_muNEta = (TH1D*)infile->Get("h_muNEta");
  TH1D *h_JMass = (TH1D*)infile->Get("h_JMass");
  TH1D *h_JPt = (TH1D*)infile->Get("h_JPt");
  TH1D *h_Jy = (TH1D*)infile->Get("h_Jy");
  TH1D *h_Jlts = (TH1D*)infile->Get("h_Jlts");
  TH1D *h_BMass = (TH1D*)infile->Get("h_BMass");
  
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLogy();
  double norm;

  norm = h_muPpT->Integral();
  h_muPpT->Scale(1./norm);
  h_muNpT->Scale(1./norm);
  h_muPpT->GetYaxis()->SetRangeUser(1e-6, 1e-1);
  h_muPpT->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
  h_muPpT->SetLineColor(kRed);
  h_muPpT->SetTitle("B #rightarrow J/#psi K muon p_{T}");
  h_muPpT->Draw("hist");
  h_muNpT->SetLineColor(kBlue);
  h_muNpT->Draw("hist same");
  TLine *mupT = new TLine(5.6, 1e-6, 5.6, 1e-1);
  mupT->SetLineStyle(kDashed);
  mupT->SetLineColor(kBlack);
  mupT->Draw("lsame");
  c->SaveAs("plots/exc_muon_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_muPEta->Integral();
  h_muPEta->Scale(1./norm);
  h_muNEta->Scale(1./norm);
  h_muPEta->GetYaxis()->SetRangeUser(0, 0.035);
  h_muPEta->GetXaxis()->SetTitle("#eta(#mu)");
  h_muPEta->SetLineColor(kRed);
  h_muPEta->SetTitle("B #rightarrow J/#psi K muon #eta");
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
  c->SaveAs("plots/exc_muon_eta.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_JMass->Integral();
  h_JMass->Scale(1./norm);
  h_JMass->GetYaxis()->SetRangeUser(0, 0.07);
  h_JMass->GetXaxis()->SetTitle("M(J/#psi) (GeV)");
  h_JMass->SetTitle("B #rightarrow J/#psi K J/#psi mass");
  h_JMass->Draw("hist");
  TLine *JMass_1 = new TLine(3, 0, 3, 0.07);
  JMass_1->SetLineStyle(kDashed);
  JMass_1->SetLineColor(kBlack);
  JMass_1->Draw("lsame");
  TLine *JMass_2 = new TLine(3.2, 0, 3.2, 0.07);
  JMass_2->SetLineStyle(kDashed);
  JMass_2->SetLineColor(kBlack);
  JMass_2->Draw("lsame");
  c->SaveAs("plots/exc_jpsi_mass.pdf");
  c->Clear();

  c->SetLogy();
  norm = h_JPt->Integral();
  h_JPt->Scale(1./norm);
  h_JPt->GetYaxis()->SetRangeUser(4e-6, 5e-1);
  h_JPt->GetXaxis()->SetTitle("p_{T}(J/#psi) (GeV)");
  h_JPt->SetTitle("B #rightarrow J/#psi K J/#psi p_{T}");
  h_JPt->Draw("hist");
  TLine *JPt_1 = new TLine(25, 4e-6, 25, 5e-1);
  JPt_1->SetLineStyle(kDashed);
  JPt_1->SetLineColor(kBlack);
  JPt_1->Draw("lsame");
  TLine *JPt_2 = new TLine(100, 4e-6, 100, 5e-1);
  JPt_2->SetLineStyle(kDashed);
  JPt_2->SetLineColor(kBlack);
  JPt_2->Draw("lsame");
  c->SaveAs("plots/exc_jpsi_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_Jy->Integral();
  h_Jy->Scale(1./norm);
  h_Jy->GetYaxis()->SetRangeUser(0, 0.03);
  h_Jy->GetXaxis()->SetTitle("y(J/#psi)");
  h_Jy->SetTitle("B #rightarrow J/#psi K J/#psi y");
  h_Jy->Draw("hist");
  TLine *JRap_1 = new TLine(-1.2, 0, -1.2, 0.03);
  JRap_1->SetLineStyle(kDashed);
  JRap_1->SetLineColor(kBlack);
  JRap_1->Draw("lsame");
  TLine *JRap_2 = new TLine(1.2, 0, 1.2, 0.03);
  JRap_2->SetLineStyle(kDashed);
  JRap_2->SetLineColor(kBlack);
  JRap_2->Draw("lsame");
  c->SaveAs("plots/exc_jpsi_rap.pdf");
  c->Clear();

  c->SetLogy();
  norm = h_Jlts->Integral();
  h_Jlts->Scale(1./norm);
  h_Jlts->GetYaxis()->SetRangeUser(1e-3, 3e-1);
  h_Jlts->GetXaxis()->SetTitle("|c#tau|/#sigma_{c#tau}");
  h_Jlts->SetTitle("B #rightarrow J/#psi K J/#psi lts");
  h_Jlts->Draw("hist");
  TLine *Jlts = new TLine(2.5, 1e-3, 2.5, 3e-1);
  Jlts->SetLineStyle(kDashed);
  Jlts->SetLineColor(kBlack);
  Jlts->Draw("lsame");
  c->SaveAs("plots/exc_jpsi_lts.pdf");

  c->SetLogy(0);
  norm = h_BMass->Integral();
  h_BMass->Scale(1./norm);
  h_BMass->GetYaxis()->SetRangeUser(0, 0.05);
  h_BMass->GetXaxis()->SetTitle("M(B) (GeV)");
  h_BMass->SetTitle("B #rightarrow J/#psi K B mass");
  h_BMass->Draw("hist");
  TLine *BMass_1 = new TLine(5.24, 0, 5.24, 0.05);
  BMass_1->SetLineStyle(kDashed);
  BMass_1->SetLineColor(kBlack);
  BMass_1->Draw("lsame");
  TLine *BMass_2 = new TLine(5.32, 0, 5.32, 0.05);
  BMass_2->SetLineStyle(kDashed);
  BMass_2->SetLineColor(kBlack);
  BMass_2->Draw("lsame");
  c->SaveAs("plots/exc_B_mass.pdf");

  infile->Close();
  
}
