// code to plot the cut variables for the MC Jpsi sample
/* variables to plot
- single muon pT, eta
- dimuon mass, pT, y, ct/cterr
*/

void compVarMC()
{
  TFile *infile = new TFile("../../Store_data_codes/2018/Jpsi_MC_comp.root");
  TH1D *h_muPpT_0 = (TH1D*)infile->Get("h_muPpT_0");
  TH1D *h_muNpT_0 = (TH1D*)infile->Get("h_muNpT_0");
  TH1D *h_muPEta_0 = (TH1D*)infile->Get("h_muPEta_0");
  TH1D *h_muNEta_0 = (TH1D*)infile->Get("h_muNEta_0");
  TH1D *h_JMass_0 = (TH1D*)infile->Get("h_JMass_0");
  TH1D *h_JPt_0 = (TH1D*)infile->Get("h_JPt_0");
  TH1D *h_Jy_0 = (TH1D*)infile->Get("h_Jy_0");
  
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLogy();
  double norm;

  norm = h_muPpT_0->Integral();
  h_muPpT_0->Scale(1./norm);
  h_muNpT_0->Scale(1./norm);
  h_muPpT_0->GetYaxis()->SetRangeUser(1e-3, 1e-1);
  h_muPpT_0->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
  h_muPpT_0->SetLineColor(kRed);
  h_muPpT_0->SetTitle("Inclusive MC (low p_{T}) muon p_{T}");
  h_muPpT_0->Draw("hist");
  h_muNpT_0->SetLineColor(kBlue);
  h_muNpT_0->Draw("hist same");
  c->SaveAs("plots/MC_comp_muon_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_muPEta_0->Integral();
  h_muPEta_0->Scale(1./norm);
  h_muNEta_0->Scale(1./norm);
  h_muPEta_0->GetYaxis()->SetRangeUser(0, 0.02);
  h_muPEta_0->GetXaxis()->SetTitle("#eta(#mu)");
  h_muPEta_0->SetLineColor(kRed);
  h_muPEta_0->SetTitle("Inclusive MC (low p_{T}) muon #eta");
  h_muPEta_0->Draw("hist");
  h_muNEta_0->SetLineColor(kBlue);
  h_muNEta_0->Draw("hist same");
  c->SaveAs("plots/MC_comp_muon_eta.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_JMass_0->Integral();
  h_JMass_0->Scale(1./norm);
  h_JMass_0->GetYaxis()->SetRangeUser(0, 0.07);
  h_JMass_0->GetXaxis()->SetTitle("M(J/#psi) (GeV)");
  h_JMass_0->SetTitle("Inclusive MC (low p_{T}) J/#psi mass");
  h_JMass_0->Draw("hist");
  c->SaveAs("plots/MC_comp_jpsi_mass.pdf");
  c->Clear();

  c->SetLogy();
  norm = h_JPt_0->Integral();
  h_JPt_0->Scale(1./norm);
  h_JPt_0->GetYaxis()->SetRangeUser(1e-3, 1e-1);
  h_JPt_0->GetXaxis()->SetTitle("p_{T}(J/#psi) (GeV)");
  h_JPt_0->SetTitle("Inclusive MC (low p_{T}) J/#psi p_{T}");
  h_JPt_0->Draw("hist");
  c->SaveAs("plots/MC_comp_jpsi_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_Jy_0->Integral();
  h_Jy_0->Scale(1./norm);
  h_Jy_0->GetYaxis()->SetRangeUser(0, 0.02);
  h_Jy_0->GetXaxis()->SetTitle("y(J/#psi)");
  h_Jy_0->SetTitle("Inclusive MC (low p_{T}) J/#psi y");
  h_Jy_0->Draw("hist");
  c->SaveAs("plots/MC_comp_jpsi_rap.pdf");
  c->Clear();

  TH1D *h_muPpT_1 = (TH1D*)infile->Get("h_muPpT_1");
  TH1D *h_muNpT_1 = (TH1D*)infile->Get("h_muNpT_1");
  TH1D *h_muPEta_1 = (TH1D*)infile->Get("h_muPEta_1");
  TH1D *h_muNEta_1 = (TH1D*)infile->Get("h_muNEta_1");
  TH1D *h_JMass_1 = (TH1D*)infile->Get("h_JMass_1");
  TH1D *h_JPt_1 = (TH1D*)infile->Get("h_JPt_1");
  TH1D *h_Jy_1 = (TH1D*)infile->Get("h_Jy_1");

  c->SetLogy();
  norm = h_muPpT_1->Integral();
  h_muPpT_1->Scale(1./norm);
  h_muNpT_1->Scale(1./norm);
  h_muPpT_1->GetYaxis()->SetRangeUser(1e-3, 1e-1);
  h_muPpT_1->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
  h_muPpT_1->SetLineColor(kRed);
  h_muPpT_1->SetTitle("Inclusive MC (high p_{T}) muon p_{T}");
  h_muPpT_1->Draw("hist");
  h_muNpT_1->SetLineColor(kBlue);
  h_muNpT_1->Draw("hist same");
  c->SaveAs("plots/MC_hpt_comp_muon_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_muPEta_1->Integral();
  h_muPEta_1->Scale(1./norm);
  h_muNEta_1->Scale(1./norm);
  h_muPEta_1->GetYaxis()->SetRangeUser(0, 0.02);
  h_muPEta_1->GetXaxis()->SetTitle("#eta(#mu)");
  h_muPEta_1->SetLineColor(kRed);
  h_muPEta_1->SetTitle("Inclusive MC (high p_{T}) muon #eta");
  h_muPEta_1->Draw("hist");
  h_muNEta_1->SetLineColor(kBlue);
  h_muNEta_1->Draw("hist same");
  c->SaveAs("plots/MC_hpt_comp_muon_eta.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_JMass_1->Integral();
  h_JMass_1->Scale(1./norm);
  h_JMass_1->GetYaxis()->SetRangeUser(0, 0.07);
  h_JMass_1->GetXaxis()->SetTitle("M(J/#psi) (GeV)");
  h_JMass_1->SetTitle("Inclusive MC (high p_{T}) J/#psi mass");
  h_JMass_1->Draw("hist");
  c->SaveAs("plots/MC_hpt_comp_jpsi_mass.pdf");
  c->Clear();

  c->SetLogy();
  norm = h_JPt_1->Integral();
  h_JPt_1->Scale(1./norm);
  h_JPt_1->GetYaxis()->SetRangeUser(1e-3, 1e-1);
  h_JPt_1->GetXaxis()->SetTitle("p_{T}(J/#psi) (GeV)");
  h_JPt_1->SetTitle("Inclusive MC (high p_{T}) J/#psi p_{T}");
  h_JPt_1->Draw("hist");
  c->SaveAs("plots/MC_hpt_comp_jpsi_pt.pdf");
  c->Clear();

  c->SetLogy(0);
  norm = h_Jy_1->Integral();
  h_Jy_1->Scale(1./norm);
  h_Jy_1->GetYaxis()->SetRangeUser(0, 0.02);
  h_Jy_1->GetXaxis()->SetTitle("y(J/#psi)");
  h_Jy_1->SetTitle("Inclusive MC (high p_{T}) J/#psi y");
  h_Jy_1->Draw("hist");
  c->SaveAs("plots/MC_hpt_comp_jpsi_rap.pdf");
  c->Clear();

  infile->Close();
  
}
