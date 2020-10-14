// function that makes sure all the fit sections are continuous
double connS(double *xvar, double *pvar)
{
  // get variables
  double y = xvar[0];
  double c1 = pvar[0], c2 = pvar[1], m = pvar[2];

  // define fit limits
  double y1 = 0.25, y2 = 0.5, y3 = 0.7;
  
  // convert to easy system
  double m1 = (c2-c1)/(y2-y1);
  double b1 = c1-m1*y1;
  double m2 = m;
  double b2 = c2-m2*y3;

  // return the correct function depending on input y
  if(y < 0 || y > 1.2) {
    cout << "incorrect y input" << endl;
    return 0;
  }
  
  else if(y <= y1) return c1;
  else if(y <= y2) return b1+m1*y;
  else if(y <= y3) return c2;
  else return b2 + m2*y; 
}

// fit mass distributions in several pt and |y| bins,
// get sigma and mu distribution, plot
// fit sigma distribution, save results
void fitMassSigma()
{
  // define pT and |y| binning
  int pt_bins = 7;
  double pt_min[pt_bins], pt_max[pt_bins];
  for(int i = 0; i < pt_bins; i++) {
    pt_min[i] = 25.+10.*i;
    pt_max[i] = 25+10.*(i+1.);
  }
  pt_max[1] = 46.;
  pt_min[2] = 46.;
  pt_max[3] = 66.;
  pt_min[4] = 66.;
  pt_max[6] = 100.;
  
  int y_bins = 12;
  double y_min[] = {0, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.7, 0.8, 0.9, 1.0, 1.1};
  double y_max[] = {0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2};
  
  // define histograms for mass distributions
  double mMin = 3.05, mMax = 3.15;
  TH1D*** h_m = new TH1D**[pt_bins];
  for(int j = 0; j < pt_bins; j++) {
    h_m[j] = new TH1D*[y_bins];
    for(int i = 0; i < y_bins; i++) {
      h_m[j][i] = new TH1D(Form("h_m%d_%d", j+1, i+1), Form("%.2f<|y|<%.2f, %.0f < p_{T} < %.0f GeV", y_min[i], y_max[i], pt_min[j], pt_max[j]), 100, mMin, mMax);
    }
  }

  // open file, get tree
  TFile *fin1 = new TFile("../Store_data_codes/MC18_cos.root");
  TTree *tree1 = (TTree*)fin1->Get("MC_cos");
  TFile *fin2 = new TFile("../Store_data_codes/MC18_hpt_cos.root");
  TTree *tree2 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../Store_data_codes/MC18_vhpt_cos.root");
  TTree *tree3 = (TTree*)fin3->Get("MC_cos");

  int nEvt, perc;
  double y, mass, pt;

  nEvt = tree1->GetEntries();
  perc = nEvt/100;
  tree1->SetBranchAddress("JpsiMass", &mass);
  tree1->SetBranchAddress("JpsiRap", &y);
  tree1->SetBranchAddress("JpsiPt", &pt);

  // fill histograms
  for(int i_evt = 0; i_evt < nEvt; i_evt++) {
    tree1->GetEntry(i_evt);
    if(mass > mMin && mass < mMax && pt > pt_min[0] && pt < 46) {
      for(int i_pt = 0; i_pt < pt_bins; i_pt++) {
	if(pt > pt_min[i_pt] && pt < pt_max[i_pt]) {
	  for(int i_y = 0; i_y < y_bins; i_y++) { 
	    if(y > y_min[i_y] && y < y_max[i_y]) {
	      h_m[i_pt][i_y]->Fill(mass);
	    }
	  }
	}
      }
    }
    if((i_evt+1)%perc == 0) cout << (i_evt+1)/perc << "% done with MC 1" << endl; 
  }
  fin1->Close();

  nEvt = tree2->GetEntries();
  perc = nEvt/100;
  tree2->SetBranchAddress("JpsiMass", &mass);
  tree2->SetBranchAddress("JpsiRap", &y);
  tree2->SetBranchAddress("JpsiPt", &pt);

  // fill histograms
  for(int i_evt = 0; i_evt < nEvt; i_evt++) {
    tree2->GetEntry(i_evt);
    if(mass > mMin && mass < mMax && pt > 46 && pt < 66) {
      for(int i_pt = 0; i_pt < pt_bins; i_pt++) {
	if(pt > pt_min[i_pt] && pt < pt_max[i_pt]) {
	  for(int i_y = 0; i_y < y_bins; i_y++) { 
	    if(y > y_min[i_y] && y < y_max[i_y]) {
	      h_m[i_pt][i_y]->Fill(mass);
	    }
	  }
	}
      }
    }
    if((i_evt+1)%perc == 0) cout << (i_evt+1)/perc << "% done with MC 2" << endl; 
  }
  fin2->Close();

  nEvt = tree3->GetEntries();
  perc = nEvt/100;
  tree3->SetBranchAddress("JpsiMass", &mass);
  tree3->SetBranchAddress("JpsiRap", &y);
  tree3->SetBranchAddress("JpsiPt", &pt);

  // fill histograms
  for(int i_evt = 0; i_evt < nEvt; i_evt++) {
    tree3->GetEntry(i_evt);
    if(mass > mMin && mass < mMax && pt > 66 && pt < pt_max[pt_bins-1]) {
      for(int i_pt = 0; i_pt < pt_bins; i_pt++) {
	if(pt > pt_min[i_pt] && pt < pt_max[i_pt]) {
	  for(int i_y = 0; i_y < y_bins; i_y++) { 
	    if(y > y_min[i_y] && y < y_max[i_y]) {
	      h_m[i_pt][i_y]->Fill(mass);
	    }
	  }
	}
      }
    }
    if((i_evt+1)%perc == 0) cout << (i_evt+1)/perc << "% done with MC 3" << endl; 
  }
  fin3->Close();

  // plot and fit mass histos, get sigma and mu histos
  TCanvas *c = new TCanvas("", "", 700, 700);

  // binning for sigma / mu histos
  double pt_h[pt_bins+1], y_h[y_bins+1];
  for(int i = 0; i < pt_bins; i++)
    pt_h[i] = pt_min[i];
  pt_h[pt_bins] = pt_max[pt_bins-1];
  for(int i = 0; i < y_bins; i++)
    y_h[i] = y_min[i];
  y_h[y_bins] = y_max[y_bins-1];       

  TH1D **h_spt = new TH1D*[pt_bins];
  TH1D **h_mpt = new TH1D*[pt_bins];

  // cycle in pt + y
  for(int i_pt = 0; i_pt < pt_bins; i_pt++) {
    // initialize sigma and mu histos
    h_spt[i_pt] = new TH1D(Form("h_sigp%d", i_pt+1), "#sigma(|y|) in p_{T} bins", y_bins, y_h);
    h_mpt[i_pt] = new TH1D(Form("h_mup%d", i_pt+1), "#mu(|y|) in p_{T} bins", y_bins, y_h);
    
    for(int i_y = 0; i_y < y_bins; i_y++) {

      // plotting mass histos
      h_m[i_pt][i_y]->SetLineColor(kBlack);
      h_m[i_pt][i_y]->SetMarkerColor(kBlack);
      h_m[i_pt][i_y]->SetStats(0);
      h_m[i_pt][i_y]->GetXaxis()->SetTitle("M^{J/#psi} (GeV)");
      h_m[i_pt][i_y]->Draw("error");

      // fit mass with gaussian
      TF1 *f3 = new TF1("f3", "[0]*TMath::Gaus(x, [1],[2])", mMin ,mMax);
      f3->SetParameters(100, 3.097, 0.02);
      f3->SetLineColor(kViolet);
      h_m[i_pt][i_y]->Fit("f3", "Q");

      TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg->SetTextSize(0.03);
      leg->AddEntry(h_m[i_pt][i_y], "Data", "pl");
      leg->AddEntry(f3, "Gaussian", "l");
      leg->Draw();

      // get sigma and mu from fit
      h_spt[i_pt]->SetBinContent(i_y+1, f3->GetParameter(2));
      h_spt[i_pt]->SetBinError(i_y+1, f3->GetParError(2));
      h_mpt[i_pt]->SetBinContent(i_y+1, f3->GetParameter(1));
      h_mpt[i_pt]->SetBinError(i_y+1, f3->GetParError(1));

      c->SaveAs(Form("plots_fit/mass_pT%d_y%d.pdf", i_pt+1, i_y+1));
      c->Clear();
    }
    h_spt[i_pt]->SetLineColor(i_pt+1);
    h_spt[i_pt]->SetMarkerColor(i_pt+1);
    h_spt[i_pt]->SetStats(0);
    h_mpt[i_pt]->SetLineColor(i_pt+1);
    h_mpt[i_pt]->SetMarkerColor(i_pt+1);
    h_mpt[i_pt]->SetStats(0);
  }

  // plotting sigma histos
  h_spt[0]->GetYaxis()->SetRangeUser(0.015, 0.045);
  h_spt[0]->GetXaxis()->SetTitle("|y|");
  h_spt[0]->GetYaxis()->SetTitle("#sigma (GeV)");
  h_spt[0]->Draw("error");
  for(int i = 1; i < pt_bins; i++) {
    h_spt[i]->Draw("same");
  }

  TLegend *leg = new TLegend(0.3, 0.55, 0.55, 0.9);
  leg->SetTextSize(0.03);
  for(int i = 0; i < pt_bins; i++)
    leg->AddEntry(h_spt[i], Form("[%.0f, %.0f] GeV", pt_min[i], pt_max[i]), "pl");
  leg->Draw();
  
  c->SaveAs("plots_fit/sigma_fit.pdf");
  c->Clear();

  // define continuous fit function, fit sigma
  TF1 **cont = new TF1*[pt_bins];
  cont[0] = new TF1("conti0", "connS", 0., 1.2, 3, 1);
  cont[0]->SetParameters(0.22, 0.25, 0.1);
  cont[0]->SetParNames("C1", "C2", "M2");
  cont[0]->SetLineColor(1);
  h_spt[0]->Draw("error");
  h_spt[0]->Fit("conti0");
  for(int i = 1; i < pt_bins; i++) {
    h_spt[i]->Draw("same");
    cont[i] = new TF1(Form("conti%d", i), "connS", 0., 1.2, 3, 1);
    cont[i]->SetParameters(0.22, 0.25, 0.1);
    cont[i]->SetParNames("C1", "C2", "M2");
    cont[i]->SetLineColor(i+1);
    h_spt[i]->Fit(Form("conti%d", i));
  }
  for(int i = 0; i < pt_bins; i++) cont[i]->Draw("lsame");

  leg->Draw();
  
  c->SaveAs("plots_fit/sigma_cont.pdf");
  c->Clear();

  // plot mu histos
  h_mpt[0]->GetYaxis()->SetRangeUser(3.09, 3.11);
  h_mpt[0]->GetXaxis()->SetTitle("|y|");
  h_mpt[0]->GetYaxis()->SetTitle("#mu (GeV)");
  h_mpt[0]->Draw("error");
  for(int i = 1; i < pt_bins; i++) {
    h_mpt[i]->Draw("same");
  }

  leg->Draw();
  
  c->SaveAs("plots_fit/mu_fit.pdf");
  c->Clear();
  
  // output tex and txt file with fit results
  ofstream tex;
  tex.open("text_output/fit_sig.tex");
  tex << "\\begin{table}[h!]" << endl;
  tex << "\\centering" << endl;
  tex << "\\begin{tabular}{c||c|c|c||c}" << endl;
  tex << "$\\pt$ bin & $C_1$ (MeV) & $C_2$ (MeV) & $m$ (MeV) & $\\chi^2/$ndf \\\\" << endl;
  tex << "\\hline" << endl;
  for(int i = 0; i < pt_bins; i++) {
    tex << Form("$[%.0f, %.0f]$ GeV & ", pt_min[i], pt_max[i]);
    for(int j = 0; j < 3; j++) {
      int l_norm = ceil(-log10(cont[i]->GetParError(j)*1000))+1;

      tex << "$" << setprecision(l_norm) << fixed << cont[i]->GetParameter(j)*1000 << "\\pm" << cont[i]->GetParError(j)*1000 << "$ & ";
    }
    tex << setprecision(1) << fixed << cont[i]->GetChisquare() << "/" << cont[i]->GetNDF() << "\\\\" << endl;
  }
  tex << "\\end{tabular}" << endl;
  tex << "\\caption{Fit results for each of the 7 $\\pt$ bins}" << endl;
  tex << "\\label{t:fit}" << endl;
  tex << "\\end{table}" << endl;
  tex.close();

  ofstream txt;
  txt.open("text_output/fit_sig.txt");
  txt << Form("y1=%.2f\ty2=%.1f\ty3=%.1f", 0.25, 0.5, 0.7) << endl;
  txt << "pt_min\tpt_max\tc1\tec1\tc2\tec2\tm\tem" << endl;
  for(int i = 0; i < pt_bins; i++) {
    txt << pt_min[i] << "\t" << pt_max[i];
    for(int j = 0; j < 3; j++)
      txt << "\t" << cont[i]->GetParameter(j) << "\t" << cont[i]->GetParError(j);
    txt << endl;
  }
  txt.close();
}
