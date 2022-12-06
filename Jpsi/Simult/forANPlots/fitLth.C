// code to fit the lth results

#import "../ptbins.C"

void fitLth()
{
  // get the fit results
  // get lambda values for each bin
  TFile *fIndS = new TFile("../Systematics/mainDiffs/files/finalUnc.root");
  TGraphAsymmErrors* graph_lthS = (TGraphAsymmErrors*)fIndS->Get("lth_fPR");
  fIndS->Close();

  // remove the 2017-2018 factor from the total systematic uncertainty
  double sigY = 0.012;
  double new_sigP[nPtBins], new_sigM[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    double sigP = graph_lthS->GetEYhigh()[i];
    double sigM = graph_lthS->GetEYlow()[i];

    new_sigP[i] = sqrt(pow(sigP,2) - pow(sigY,2));
    new_sigM[i] = sqrt(pow(sigM,2) - pow(sigY,2));

  }
  TGraphAsymmErrors* g_lth = new TGraphAsymmErrors(graph_lthS->GetN(), graph_lthS->GetX(), graph_lthS->GetY(), graph_lthS->GetEXlow(), graph_lthS->GetEXhigh(), new_sigM, new_sigP);
  
  // plot both fits
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);
  
  TF1 *fcon = new TF1("fcon", "[0]", ptBins[0], ptBins[nPtBins]);
  fcon->SetParameter(0, 0.2);
  fcon->SetLineColor(kBlue);
  fcon->SetLineStyle(kDashed);
  g_lth->Fit(fcon, "0");

  TF1 *flin = new TF1("flin", "[0]+x*[1]", ptBins[0], ptBins[nPtBins]);
  flin->SetParameters(0.1, 0.1);
  flin->SetLineColor(kViolet);
  flin->SetLineStyle(kDashed);
  g_lth->Fit(flin, "0");
  
  // draw the fit results

  // draw lambda_th(pT)
  TH1F *fc = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("#lambda_{#theta}");
  fc->GetYaxis()->SetTitleOffset(1.3);
  fc->GetYaxis()->SetLabelOffset(0.01);
    
  g_lth->SetLineColor(kBlack);
  g_lth->SetMarkerStyle(20);
  g_lth->SetMarkerSize(.5);
  g_lth->SetMarkerColor(kBlack);
  g_lth->Draw("p same");

  fcon->Draw("lsame");
 
  TLine *zero = new TLine(ptBins[0]-5, 0, ptBins[nPtBins], 0);
  zero->SetLineColor(kBlack);
  zero->SetLineStyle(kDashed);
  zero->Draw();
  
  TLatex lc;
  lc.SetTextSize(0.03);
  lc.SetTextColor(kBlue);
  lc.DrawLatex(30, -0.3, Form("b = %.3f #pm %.3f", fcon->GetParameter(0), fcon->GetParError(0)));
  lc.DrawLatex(30, -0.45, Form("#chi^{2}/ndf = %.1f/%d", fcon->GetChisquare(), fcon->GetNDF()));
  lc.DrawLatex(30, -0.6, Form("P(#chi^{2},ndf) = %.1f %%", TMath::Prob(fcon->GetChisquare(), fcon->GetNDF())*100.));
  lc.DrawLatex(35, 0.05, "#lambda_{#theta} = b");
    
  c->SaveAs("plots/fit_lth_con.pdf");
  c->Clear();

  TH1F *fl = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#lambda_{#theta}");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
    
  g_lth->SetLineColor(kBlack);
  g_lth->SetMarkerStyle(20);
  g_lth->SetMarkerSize(.5);
  g_lth->SetMarkerColor(kBlack);
  g_lth->Draw("p same");

  flin->Draw("lsame");
 
  zero->Draw();
  
  TLatex ll;
  ll.SetTextSize(0.03);
  ll.SetTextColor(kViolet);
  ll.DrawLatex(30, -0.15, Form("m = ( %.2f #pm %.2f ) #times 10^{-3}", flin->GetParameter(1)*1e3, flin->GetParError(1)*1e3));
  ll.DrawLatex(30, -0.3, Form("b = %.3f #pm %.3f", flin->GetParameter(0), flin->GetParError(0)));
  ll.DrawLatex(30, -0.45, Form("#chi^{2}/ndf = %.1f/%d", flin->GetChisquare(), flin->GetNDF()));
  ll.DrawLatex(30, -0.6, Form("P(#chi^{2},ndf) = %.1f %%", TMath::Prob(flin->GetChisquare(), flin->GetNDF())*100.));
  ll.DrawLatex(35, 0.05, "#lambda_{#theta} = m*p_{T} + b");
  
  c->SaveAs("plots/fit_lth_lin.pdf");
  c->Clear();  
  c->Destructor();
}
