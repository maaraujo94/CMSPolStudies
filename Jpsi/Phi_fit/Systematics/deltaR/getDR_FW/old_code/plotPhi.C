void plotPhi()
{
  TFile *fin = new TFile("Phistore.root");
  TH2D* r_PhiPM = (TH2D*)fin->Get("r_PhiPM");
  TH1D* r_PhiP  = (TH1D*)fin->Get("r_PhiP");
  TH1D* r_PhiM  = (TH1D*)fin->Get("r_PhiM");
  TH2D* t_PhiPM = (TH2D*)fin->Get("t_PhiPM"); 
  TH1D* t_PhiP  = (TH1D*)fin->Get("t_PhiP");
  TH1D* t_PhiM  = (TH1D*)fin->Get("t_PhiM");

  TCanvas *c = new TCanvas("", "", 600, 600);

  r_PhiPM->SetStats(0);
  r_PhiPM->SetTitle("reco #phi, p_ {T} > 50");
  r_PhiPM->Draw("colz");

  c->SaveAs("phiPM_r.pdf");
  c->Clear();

  t_PhiPM->SetStats(0);
  t_PhiPM->SetTitle("trig #phi, p_ {T} > 50");
  t_PhiPM->Draw("colz");
  
  c->SaveAs("phiPM_t.pdf");
  c->Clear();

  r_PhiP->SetStats(0);
  r_PhiP->SetTitle("reco #phi, p_ {T} > 50");
  r_PhiP->SetLineColor(kRed);
  r_PhiP->GetXaxis()->SetTitle("#phi");
  r_PhiP->Draw();

  r_PhiM->SetLineColor(kBlue);
  r_PhiM->Draw("same");

  c->SaveAs("phi_r.pdf");
  c->Clear();

  t_PhiP->SetStats(0);
  t_PhiP->SetTitle("trig #phi, p_ {T} > 50");
  t_PhiP->SetLineColor(kRed);
  t_PhiP->GetXaxis()->SetTitle("#phi");
  t_PhiP->Draw();

  t_PhiM->SetLineColor(kBlue);
  t_PhiM->Draw("same");

  c->SaveAs("phi_t.pdf");
  c->Clear();

  c->Destructor();
  fin->Close();

}
