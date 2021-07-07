// code to plot the cut variables for the MC Jpsi low-pT samples
/* variables to plot
   - single muon pT, eta
   - dimuon mass, pT, y, ct/cterr
*/

void compHptMC()
{
  string lbl[] = {"hpt", "vhpt"};
  double norm;

  TH1D **h_muPpT = new TH1D*[2];
  TH1D **h_muNpT = new TH1D*[2];
  TH1D **h_muPEta = new TH1D*[2];
  TH1D **h_muNEta = new TH1D*[2];
  TH1D **h_JMass = new TH1D*[2];
  TH1D **h_JPt = new TH1D*[2];
  TH1D **h_Jy = new TH1D*[2];
  TH1D **h_Jlt = new TH1D*[2];

  TFile *infile = new TFile("../Store_data_codes/Jpsi_MCh_comp.root");
  TCanvas *c = new TCanvas("", "", 700, 700);
  // cycle over 2017, 2018
  for(int i_y = 0; i_y < 2; i_y++) {
    for(int i_s = 0; i_s < 2; i_s++) {
      h_muPpT[i_s]  = (TH1D*)infile->Get(Form("h%d_muPpT_%s",  i_y+7, lbl[i_s].c_str()));
      h_muNpT[i_s]  = (TH1D*)infile->Get(Form("h%d_muNpT_%s",  i_y+7, lbl[i_s].c_str()));
      h_muPEta[i_s] = (TH1D*)infile->Get(Form("h%d_muPEta_%s", i_y+7, lbl[i_s].c_str()));
      h_muNEta[i_s] = (TH1D*)infile->Get(Form("h%d_muNEta_%s", i_y+7, lbl[i_s].c_str()));
      h_JMass[i_s]  = (TH1D*)infile->Get(Form("h%d_JMass_%s",  i_y+7, lbl[i_s].c_str()));
      h_JPt[i_s]    = (TH1D*)infile->Get(Form("h%d_JPt_%s",    i_y+7, lbl[i_s].c_str()));
      h_Jy[i_s]     = (TH1D*)infile->Get(Form("h%d_Jy_%s",     i_y+7, lbl[i_s].c_str()));
      h_Jlt[i_s]    = (TH1D*)infile->Get(Form("h%d_Jlt_%s",    i_y+7, lbl[i_s].c_str()));
    }
    
    c->SetLogy();
    for(int i_s = 0; i_s < 2; i_s++) {
      norm = h_muPpT[i_s]->Integral();
      h_muPpT[i_s]->Scale(1./norm);
      h_muNpT[i_s]->Scale(1./norm);
      h_muPpT[i_s]->GetYaxis()->SetRangeUser(1e-5, 3e-2);
      h_muPpT[i_s]->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
      h_muPpT[i_s]->SetLineColor(kRed);
      h_muPpT[i_s]->SetStats(0);
      h_muPpT[i_s]->SetTitle(Form("201%d MC (high p_{T}) muon p_{T}", i_y+7));
      if(i_s==0) h_muPpT[i_s]->Draw("hist");
      else {
	h_muPpT[i_s]->SetLineStyle(kDashed);
	h_muNpT[i_s]->SetLineStyle(kDashed);
	h_muPpT[i_s]->Draw("hist same");
      }
      h_muNpT[i_s]->SetLineColor(kBlue);
      h_muNpT[i_s]->Draw("hist same");
    }
    c->SaveAs(Form("plots_comp/MC%d_hpt_muon_pt.pdf", i_y+7));
    c->Clear();

    c->SetLogy(0);
    for(int i_s = 0; i_s < 2; i_s++) {
      norm = h_muPEta[i_s]->Integral();
      h_muPEta[i_s]->Scale(1./norm);
      h_muNEta[i_s]->Scale(1./norm);
      h_muPEta[i_s]->GetYaxis()->SetRangeUser(0, 0.016);
      h_muPEta[i_s]->GetXaxis()->SetTitle("#eta(#mu)");
      h_muPEta[i_s]->SetLineColor(kRed);
      h_muPEta[i_s]->SetStats(0);
      h_muPEta[i_s]->SetTitle(Form("201%d MC (high p_{T}) muon #eta", i_y+7));
      if(i_s == 0) h_muPEta[i_s]->Draw("hist");
      else {
	h_muPEta[i_s]->SetLineStyle(kDashed);
	h_muNEta[i_s]->SetLineStyle(kDashed);
	h_muPEta[i_s]->Draw("hist same");
      }
      h_muNEta[i_s]->SetLineColor(kBlue);
      h_muNEta[i_s]->Draw("hist same");
    }
    c->SaveAs(Form("plots_comp/MC%d_hpt_muon_eta.pdf", i_y+7));
    c->Clear();

    c->SetLogy(0);
    for(int i_s = 0; i_s < 2; i_s++) {
      norm = h_JMass[i_s]->Integral();
      h_JMass[i_s]->Scale(1./norm);
      h_JMass[i_s]->GetYaxis()->SetRangeUser(0, 0.04);
      h_JMass[i_s]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
      h_JMass[i_s]->SetStats(0);
      h_JMass[i_s]->SetTitle(Form("201%d MC (high p_{T}) dimuon mass", i_y+7));
      if(i_s == 0) h_JMass[i_s]->Draw("hist");
      else {
	//h_JMass[i_s]->SetLineStyle(kDashed);
	h_JMass[i_s]->SetLineColor(kRed);
	h_JMass[i_s]->Draw("hist same");
      }
    }
    c->SaveAs(Form("plots_comp/MC%d_hpt_dimuon_mass.pdf", i_y+7));
    c->Clear();

    c->SetLogy();
    for(int i_s = 0; i_s < 2; i_s++) {
      norm = h_JPt[i_s]->Integral();
      h_JPt[i_s]->Scale(1./norm);
      h_JPt[i_s]->GetYaxis()->SetRangeUser(5e-3, 2e-2);
      h_JPt[i_s]->GetXaxis()->SetTitle("p_{T}(#mu#mu) (GeV)");
      h_JPt[i_s]->SetStats(0);
      h_JPt[i_s]->SetTitle(Form("201%d MC (high p_{T}) dimuon p_{T}", i_y+7));
      if(i_s == 0) h_JPt[i_s]->Draw("hist");
      else {
	//h_JPt[i_s]->SetLineStyle(kDashed);
	h_JPt[i_s]->SetLineColor(kRed);
	h_JPt[i_s]->Draw("hist same");
      }
    }
    c->SaveAs(Form("plots_comp/MC%d_hpt_dimuon_pt.pdf", i_y+7));
    c->Clear();
    
    c->SetLogy(0);
    for(int i_s = 0; i_s < 2; i_s++) {
      norm = h_Jy[i_s]->Integral();
      h_Jy[i_s]->Scale(1./norm);
      h_Jy[i_s]->GetYaxis()->SetRangeUser(0, 0.016);
      h_Jy[i_s]->GetXaxis()->SetTitle("y(#mu#mu)");
      h_Jy[i_s]->SetStats(0);
      h_Jy[i_s]->SetTitle(Form("201%d MC (high p_{T}) dimuon y", i_y+7));
      if(i_s == 0) h_Jy[i_s]->Draw("hist");
      else {
	//h_Jy[i_s]->SetLineStyle(kDashed);
	h_Jy[i_s]->SetLineColor(kRed);
	h_Jy[i_s]->Draw("hist same");
      }
    }
    c->SaveAs(Form("plots_comp/MC%d_hpt_dimuon_rap.pdf", i_y+7));
    c->Clear();

    c->SetLogy();
    for(int i_s = 0; i_s < 2; i_s++) {
      norm = h_Jlt[i_s]->Integral();
      h_Jlt[i_s]->Scale(1./norm);
      h_Jlt[i_s]->GetYaxis()->SetRangeUser(1e-5, 2e-1);
      h_Jlt[i_s]->GetXaxis()->SetTitle("c#tau(#mu#mu) (cm)");
      h_Jlt[i_s]->SetStats(0);
      h_Jlt[i_s]->SetTitle(Form("201%d MC (high p_{T}) dimuon lifetime", i_y+7));
      if(i_s == 0) h_Jlt[i_s]->Draw("hist");
      else {
	//h_Jlt[i_s]->SetLineStyle(kDashed);
	h_Jlt[i_s]->SetLineColor(kRed);
	h_Jlt[i_s]->Draw("hist same");
      }
    }
    c->SaveAs(Form("plots_comp/MC%d_hpt_dimuon_lt.pdf", i_y+7));
    c->Clear();

  }
  infile->Close();
  c->Destructor();
     
}
