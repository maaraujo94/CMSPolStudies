// code to plot the fit results

#import "../Simult/ptbins.C"

void combIndSim()
{
  // get the fit results
  // get lambda values for each bin
  const int n_lbl = 2;
  string lbl[] = {"J", "NP"};
  string sv_lbl[] = {"", "NP"};
  TFile *fInd7 = new TFile("../2017/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth7 = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lth7[i_t] = (TGraphErrors*)fInd7->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd7->Close();
  TFile *fInd8 = new TFile("../2018/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lth8 = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lth8[i_t] = (TGraphErrors*)fInd8->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fInd8->Close();
  TFile *fIndS = new TFile("../Simult/PR_fit/files/finalFitRes.root");
  TGraphErrors **graph_lthS = new TGraphErrors*[n_lbl];
  for(int i_t = 0; i_t < n_lbl; i_t++) {
    graph_lthS[i_t] = (TGraphErrors*)fIndS->Get(Form("graph_lambda_%s", lbl[i_t].c_str()));
  }    
  fIndS->Close();

  // run everything for both labels
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  
  for(int i_l = 0; i_l < n_lbl; i_l++) {
    // get the systematic uncertainty envelope
    int nS = graph_lthS[i_l]->GetN();
    double unc_y[nS];
    for(int i = 0; i < nS; i++) {
      // TODO: find psi' corresponding value
      unc_y[i] = 0.0111;
    }
    TGraphErrors *g_unc = new TGraphErrors(nS, graph_lthS[i_l]->GetX(), graph_lthS[i_l]->GetY(), graph_lthS[i_l]->GetEX(), unc_y);

    // combine both sets of data for fitting
    int nC = graph_lth7[i_l]->GetN() * 2;
    double xvC[nC], xeC[nC], yvC[nC], yeC[nC];
    for(int i = 0; i < nC/2; i++) {
      xvC[i] = graph_lth7[i_l]->GetX()[i];
      yvC[i] = graph_lth7[i_l]->GetY()[i];
      xeC[i] = graph_lth7[i_l]->GetEX()[i];
      yeC[i] = graph_lth7[i_l]->GetEY()[i];
    
      xvC[i+nC/2] = graph_lth8[i_l]->GetX()[i];
      yvC[i+nC/2] = graph_lth8[i_l]->GetY()[i];
      xeC[i+nC/2] = graph_lth8[i_l]->GetEX()[i];
      yeC[i+nC/2] = graph_lth8[i_l]->GetEY()[i];    
    }
    TGraphErrors *g_lth = new TGraphErrors(nC, xvC, yvC, xeC, yeC);
    TGraphErrors *g_lthSF = new TGraphErrors(graph_lthS[i_l]->GetN(), graph_lthS[i_l]->GetX(), graph_lthS[i_l]->GetY(), graph_lthS[i_l]->GetEX(), graph_lthS[i_l]->GetEY());

    double n_st = 0.2;
    if(i_l == 1) n_st = -0.2;
    
    TF1 *flin = new TF1("flin", "[0]*x+[1]", ptBins[0], ptBins[nPtBins]);
    flin->SetParameters(0.1, n_st);
    flin->SetLineColor(kViolet+2);
    flin->SetLineStyle(kDashed);
    g_lth->Fit(flin);

    TF1 *flinS = new TF1("flinS", "[0]*x+[1]", ptBins[0], ptBins[nPtBins]);
    flinS->SetParameters(0.1, n_st);
    flinS->SetLineColor(kBlack);
    flinS->SetLineStyle(kDashed);
    g_lthSF->Fit(flinS);
   
    // draw the fit results

    // draw lambda_th(pT)
    TH1F *fl = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
    fl->SetXTitle("p_{T} (GeV)");
    fl->SetYTitle("#lambda_{#theta}");
    fl->GetYaxis()->SetTitleOffset(1.3);
    fl->GetYaxis()->SetLabelOffset(0.01);
    if(i_l==0) fl->SetTitle("prompt #psi(2S) #lambda_{#theta}");
    else fl->SetTitle("non-prompt #psi(2S) #lambda_{#theta}");

    graph_lth7[i_l]->SetLineColor(kBlue);
    graph_lth7[i_l]->SetMarkerStyle(20);
    graph_lth7[i_l]->SetMarkerSize(.75);
    graph_lth7[i_l]->SetMarkerColor(kBlue);
    graph_lth7[i_l]->Draw("px same");
    graph_lth8[i_l]->SetLineColor(kRed);
    graph_lth8[i_l]->SetMarkerStyle(20);
    graph_lth8[i_l]->SetMarkerSize(.75);
    graph_lth8[i_l]->SetMarkerColor(kRed);
    graph_lth8[i_l]->Draw("px same");
    graph_lthS[i_l]->SetLineColor(kBlack);
    graph_lthS[i_l]->SetMarkerStyle(20);
    graph_lthS[i_l]->SetMarkerSize(.5);
    graph_lthS[i_l]->SetMarkerColor(kBlack);
    graph_lthS[i_l]->Draw("p same");
 
    TLine *zero = new TLine(ptBins[0]-5, 0, ptBins[nPtBins], 0);
    zero->SetLineColor(kBlack);
    zero->SetLineStyle(kDashed);
    zero->Draw();

    TLegend *leg = new TLegend(0.82, 0.75, 0.97, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(graph_lth7[i_l], "2017", "pl");
    leg->AddEntry(graph_lth8[i_l], "2018", "pl");
    leg->AddEntry(graph_lthS[i_l], "Run2", "pl");
    leg->Draw();
  
    c->SaveAs(Form("plots/par_lth%s_all.pdf", sv_lbl[i_l].c_str()));

    g_unc->SetLineColor(kBlack);
    g_unc->SetFillColorAlpha(kBlack, 0.5);
    g_unc->Draw("3");

    c->SaveAs(Form("plots/par_lth%s_all_band.pdf", sv_lbl[i_l].c_str()));

    TH1F *fl1 = c->DrawFrame(ptBins[0]-5, -1, ptBins[nPtBins], 1);
    fl1->SetXTitle("p_{T} (GeV)");
    fl1->SetYTitle("#lambda_{#theta}");
    fl1->GetYaxis()->SetTitleOffset(1.3);
    fl1->GetYaxis()->SetLabelOffset(0.01);
    if(i_l==0) fl1->SetTitle("prompt #psi(2S) #lambda_{#theta}");
    else fl1->SetTitle("non-prompt #psi(2S) #lambda_{#theta}");

    graph_lth7[i_l]->Draw("px same");
    graph_lth8[i_l]->Draw("px same");
    graph_lthS[i_l]->Draw("p same");
 
    zero->Draw();

    leg->Draw();
  
    flin->Draw("lsame");
    flinS->Draw("lsame");

    TLatex lc;
    lc.SetTextSize(0.03);
    if(i_l == 0) {
      lc.DrawLatex(30, -0.15, Form("m = ( %.2f #pm %.2f ) #times 10^{-3}", flinS->GetParameter(0)*1e3, flinS->GetParError(0)*1e3));
      lc.DrawLatex(30, -0.3, Form("b = %.3f #pm %.3f", flinS->GetParameter(1), flinS->GetParError(1)));
      lc.DrawLatex(30, -0.45, Form("#chi^{2}/ndf = %.0f/%d", flinS->GetChisquare(), flinS->GetNDF()));
      lc.DrawLatex(60, -0.45, Form("P(#chi^{2},ndf) = %.1f %%", TMath::Prob(flinS->GetChisquare(), flinS->GetNDF())*100.));
    }
    else {
      lc.DrawLatex(30, 0.35, Form("m = ( %.2f #pm %.2f ) #times 10^{-3}", flinS->GetParameter(0)*1e3, flinS->GetParError(0)*1e3));
      lc.DrawLatex(30, 0.2, Form("b = %.3f #pm %.3f", flinS->GetParameter(1), flinS->GetParError(1)));
      lc.DrawLatex(30, 0.05, Form("#chi^{2}/ndf = %.0f/%d", flinS->GetChisquare(), flinS->GetNDF()));
      lc.DrawLatex(60, 0.05, Form("P(#chi^{2},ndf) = %.1f %%", TMath::Prob(flinS->GetChisquare(), flinS->GetNDF())*100.));
    }
    
    lc.SetTextColor(kViolet+2);
    lc.DrawLatex(30, 0.85, Form("m = ( %.2f #pm %.2f ) #times 10^{-3}", flin->GetParameter(0)*1e3, flin->GetParError(0)*1e3));
    lc.DrawLatex(30, 0.7, Form("b = %.3f #pm %.3f", flin->GetParameter(1), flin->GetParError(1)));
    lc.DrawLatex(30, 0.55, Form("#chi^{2}/ndf = %.0f/%d", flin->GetChisquare(), flin->GetNDF()));
    lc.DrawLatex(60, 0.55, Form("P(#chi^{2},ndf) = %.1f %%", TMath::Prob(flin->GetChisquare(), flin->GetNDF())*100.));

    c->SaveAs(Form("plots/fit_lth%s_all.pdf", sv_lbl[i_l].c_str()));
    c->Clear();
  }
  c->Destructor();
}
