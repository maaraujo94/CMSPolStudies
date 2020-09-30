// code to plot and fit the costh dists in each bin of pt, and store the fit results

void fit_1d()
{
  // open files and read TTrees
  TFile *fin = new TFile("../../Store_data_codes/data_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/MC_cos.root");
  TTree *treeM = (TTree*)fin2->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int mEvt = treeM->GetEntries();

  cout << dEvt << " data events after cuts and " << mEvt << " MC events after cuts" << endl;

  // prepare binning and histograms for plots
  const int nbins = 9;
  double pTbins[10] = {12, 14, 15.5, 17.5, 19, 21, 22.5, 25, 29, 70};
  double cos_lim[9] = {0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.75, 1};
  
  TH1D **datacos = new TH1D*[nbins];
  TH1D **mccos = new TH1D*[nbins];
  for(int i = 0; i < nbins; i++) {
    datacos[i] = new TH1D(Form("c_data_%d",i), Form("data costh bin %d",i), 100, -cos_lim[i], cos_lim[i]);
    mccos[i] = new TH1D(Form("c_mc_%d",i), Form("mc costh bin %d",i), 100, -cos_lim[i], cos_lim[i]);
  }
  
  // definitions to store data and MC events
  Double_t data_cos, data_pt;
  Double_t mc_cos, mc_pt;

  treeD->SetBranchAddress("costh", &data_cos);
  treeD->SetBranchAddress("JpsiPt", &data_pt);

  treeM->SetBranchAddress("costh", &mc_cos);
  treeM->SetBranchAddress("JpsiPt", &mc_pt);

  // cycle over data and MC, fill the costh histogram acc to binning
  double pT_mean[9];
  int n_pts[9];
  for(int i = 0; i < nbins; i++)
    {
      pT_mean[i] = 0;
      n_pts[i] = 0;
    }
      
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      for (int j = 0; j < nbins; j++)
	if (data_pt > pTbins[j] && data_pt < pTbins[j+1])
	  {
	    datacos[j]->Fill(data_cos);
	    pT_mean[j] += data_pt;
	    n_pts[j]++;
	  }
    }

  // getting the mean pT of each bin
  for(int i = 0; i < nbins; i++)
    pT_mean[i] /= n_pts[i];
  
  for(int i = 0; i < mEvt; i++)
    {
      treeM->GetEntry(i);
      for (int j = 0; j < nbins; j++)
	if (mc_pt > pTbins[j] && mc_pt < pTbins[j+1])
	  mccos[j]->Fill(mc_cos);
    }
  
  // plot all costh histograms
  TCanvas *c = new TCanvas("", "", 700, 700);
  for(int i = 0; i < nbins; i++) {
    c->Clear();
    datacos[i]->Draw("e");
    c->SaveAs(Form("plots/data_cosa_bin%d.pdf", i));
    c->Clear();
    mccos[i]->Draw("e");
    c->SaveAs(Form("plots/mc_cosa_bin%d.pdf", i));
  }

  // get the ratio histogram for each bin
  TH1D **ratiocos = new TH1D*[nbins];

  for(int i = 0; i < nbins; i++)
    {
      ratiocos[i] = (TH1D*)datacos[i]->Clone(Form("ratio %d",i));
      ratiocos[i]->Sumw2();
      ratiocos[i]->Divide(mccos[i]);
    }

  // fit lambda function to histogram
  TF1 *Wdist = new TF1("wdist", "[0]*(1+x*x*[1])", -1, 1);
  Wdist->SetParNames("A", "lambda_th");
  Wdist->SetParameters(10, 0.1);

  double fit_lim[9] = {0.3, 0.35, 0.4, 0.45, 0.5, 0.5, 0.55, 0.55, 0.7};
  double l_th[nbins], el_th[nbins], chisq[nbins];
  int ndf[nbins];
  
  for(int i = 0; i < nbins; i++)
    {
      ratiocos[i]->Fit("wdist", "", "", -fit_lim[i], fit_lim[i]);
      l_th[i] = Wdist->GetParameter(1);
      el_th[i] = Wdist->GetParError(1);
      chisq[i] = Wdist->GetChisquare();
      ndf[i] = Wdist->GetNDF();
      
      c->Clear();
      ratiocos[i]->SetTitle(Form("ratio bin %d",i));
      ratiocos[i]->Draw("ep");
      c->SaveAs(Form("plots/ratio_cosa_bin%d.pdf", i));
    }

  // plotting the lambda_th obtained for each bin
  double pt_lo[nbins], pt_hi[nbins];
  for(int i = 0; i < nbins; i++)
    {
      pt_lo[i] = pT_mean[i]-pTbins[i];
      pt_hi[i] = pTbins[i+1]-pT_mean[i];
    }
  
  TGraphAsymmErrors *lambda_th = new TGraphAsymmErrors(nbins, pT_mean, l_th, pt_lo, pt_hi, el_th, el_th);
  
  c->Clear();

  TH1F *fc = c->DrawFrame(9, -1, 71, 1);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("#lambda_{#theta}");
  fc->GetYaxis()->SetTitleOffset(1.3);
  c->Modified();
  c->SetTitle("");
  
  lambda_th->Draw("p");
  
  TF1 *zero = new TF1("zero", "0", 0, 100);
  zero->SetLineStyle(kDashed);
  zero->SetLineColor(kBlack);
  zero->Draw("lsame");
  c->SaveAs("plots/lambda_th.pdf");  

  //writing a tex file with information for each bin of the fit
  ofstream tex;
  tex.open("fit_data.tex");

  tex << "\\begin{table}[h!]" << endl;
  tex << "\\centering" << endl;
  tex << "\\begin{tabular}{c|c|c|c|c}" << endl;
  tex << "$\\pt$ bin & mean $\\pt$ & $\\cos\\theta_{HX}$ limit & $\\lambda_\\theta$ & $\\chi^2/$ndf \\\\" << endl;
  tex << "\\hline" << endl;
  for(int i = 0; i < nbins; i++)
    {
      int l_norm = ceil(-log10(el_th[i]))+1;
      
      tex << "$[" << setprecision(1) << fixed << pTbins[i] << "," << pTbins[i+1] << "]$ & " << setprecision(1) << fixed << pT_mean[i] << " & [" << setprecision(2) << fixed << -fit_lim[i] << "," << fit_lim[i] << "] & $" << setprecision(l_norm) << fixed << l_th[i] << "\\pm" << el_th[i] << "$ & " << setprecision(1) << fixed << chisq[i] << "/" << ndf[i] << " \\\\" << endl;
    }
  tex << "\\end{tabular}" << endl;
  tex << "\\caption{Fit parameters for each of the 9 $\\pt$ bins}" << endl;
  tex << "\\label{t:fit}" << endl;
  tex << "\\end{table}" << endl;
  tex.close();


  fin->Close();
  fin2->Close();
  
}
