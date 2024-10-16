// code to plot the fit results - HX vs CS

#import "../ptbins.C"

void plotRes_ax()
{
  // get the fit results
  string lbl[3] = {"theta", "phi", "thph"};

  // HX PR
  TFile *fInd_HXPR = new TFile("../PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lHXPR = new TGraphErrors*[3];
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lHXPR[i_t] = (TGraphErrors*)fInd_HXPR->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd_HXPR->Close();

  // HX NP
  TFile *fInd_HXNP = new TFile("../NP_fit/files/finalFitRes.root");
  TGraphErrors **graph_lHXNP = new TGraphErrors*[3];
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lHXNP[i_t] = (TGraphErrors*)fInd_HXNP->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd_HXNP->Close();

  // CS PR
  TFile *fInd_CSPR = new TFile("../../Simult2d_CS/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lCSPR = new TGraphErrors*[3];
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lCSPR[i_t] = (TGraphErrors*)fInd_CSPR->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd_CSPR->Close();

  // CS NP
  TFile *fInd_CSNP = new TFile("../../Simult2d_CS/NP_fit/files/finalFitRes.root");
  TGraphErrors **graph_lCSNP = new TGraphErrors*[3];
  for(int i_t = 0; i_t < 3; i_t++) {
    graph_lCSNP[i_t] = (TGraphErrors*)fInd_CSNP->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd_CSNP->Close();

  // draw the fit results
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.02);

  // PR lambdas
  string lbl_id[] = {"#theta", "#phi", "#theta#phi"};
  string lbl_sv[] = {"th", "ph", "tp"};
  for(int i = 0; i < 3; i++) {
    TH1F *flbd_PR = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
    flbd_PR->SetXTitle("p_{T} (GeV)");
    flbd_PR->SetYTitle(Form("#lambda_{%s}^{PR}", lbl_id[i].c_str()));
    flbd_PR->GetYaxis()->SetTitleOffset(1.3);
    flbd_PR->GetYaxis()->SetLabelOffset(0.01);
    flbd_PR->SetTitle("");
    
    graph_lHXPR[i]->SetLineColor(kBlack);
    graph_lHXPR[i]->SetMarkerStyle(20);
    graph_lHXPR[i]->SetMarkerSize(.5);
    graph_lHXPR[i]->SetMarkerColor(kBlack);
    graph_lHXPR[i]->Draw("p same");

    graph_lCSPR[i]->SetLineColor(kBlue);
    graph_lCSPR[i]->SetMarkerStyle(25);
    graph_lCSPR[i]->SetMarkerSize(.5);
    graph_lCSPR[i]->SetMarkerColor(kBlue);
    graph_lCSPR[i]->Draw("p same");
  
    TLine *zero = new TLine(ptBins[0]-5, 0, ptBins[nPtBins], 0);
    zero->SetLineColor(kBlack);
    zero->SetLineStyle(kDashed);
    zero->Draw();

    TLegend *legPR = new TLegend(0.2, 0.2, 0.5, 0.3);
    legPR->SetTextSize(0.03);
    legPR->SetBorderSize(0);
    legPR->SetFillColorAlpha(kWhite, 0);
    legPR->AddEntry(graph_lHXPR[i], Form("#lambda_{%s}^{HX}", lbl_id[i].c_str()), "pl");
    legPR->AddEntry(graph_lCSPR[i], Form("#lambda_{%s}^{CS}", lbl_id[i].c_str()), "pl");
    legPR->Draw();
  
    c->SaveAs(Form("HX_vs_CS/l%sPR_comp.pdf", lbl_sv[i].c_str()));
    c->Clear();
  }

  // NP lambdas
  for(int i = 0; i < 3; i++) {
    TH1F *flbd_NP = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
    flbd_NP->SetXTitle("p_{T} (GeV)");
    flbd_NP->SetYTitle(Form("#lambda_{%s}^{NP}", lbl_id[i].c_str()));
    flbd_NP->GetYaxis()->SetTitleOffset(1.3);
    flbd_NP->GetYaxis()->SetLabelOffset(0.01);
    flbd_NP->SetTitle("");
    
    graph_lHXNP[i]->SetLineColor(kBlack);
    graph_lHXNP[i]->SetMarkerStyle(20);
    graph_lHXNP[i]->SetMarkerSize(.5);
    graph_lHXNP[i]->SetMarkerColor(kBlack);
    graph_lHXNP[i]->Draw("p same");

    graph_lCSNP[i]->SetLineColor(kBlue);
    graph_lCSNP[i]->SetMarkerStyle(25);
    graph_lCSNP[i]->SetMarkerSize(.5);
    graph_lCSNP[i]->SetMarkerColor(kBlue);
    graph_lCSNP[i]->Draw("p same");
  
    TLine *zero = new TLine(ptBins[0]-5, 0, ptBins[nPtBins], 0);
    zero->SetLineColor(kBlack);
    zero->SetLineStyle(kDashed);
    zero->Draw();

    TLegend *legNP = new TLegend(0.2, 0.2, 0.5, 0.3);
    legNP->SetTextSize(0.03);
    legNP->SetBorderSize(0);
    legNP->SetFillColorAlpha(kWhite, 0);
    legNP->AddEntry(graph_lHXNP[i], Form("#lambda_{%s}^{HX}", lbl_id[i].c_str()), "pl");
    legNP->AddEntry(graph_lCSNP[i], Form("#lambda_{%s}^{CS}", lbl_id[i].c_str()), "pl");
    legNP->Draw();
  
    c->SaveAs(Form("HX_vs_CS/l%sNP_comp.pdf", lbl_sv[i].c_str()));
    c->Clear();
  }

  c->Destructor();
 
}
