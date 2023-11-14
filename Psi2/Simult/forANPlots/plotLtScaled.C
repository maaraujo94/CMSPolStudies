#import "../ptbins.C"

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

const int n_colors = 50;
int st_col;

int do_round(double val)
{
  int valR = (int)val;
  if (val-valR > 0.5) return valR+1;
  else return valR;
}

void init_color()
{
  const int n_stops = 5;

  float stops[n_stops] = {0, 0.25, 0.5, 0.75, 1.0};
  float red[n_stops] = {1,1,0,0,0.9};
  float green[n_stops] = {0,0.6,0.7,0.5,0.1};
  float blue[n_stops] = {0,0,0,1,0};

  // For each defined gradient...
  for (int g = 1; g < n_stops; g++) {
      // create the colors...
      int nColorsGradient = (Int_t) (floor(n_colors*stops[g]) - floor(n_colors*stops[g-1]));
      for (int c = 0; c < nColorsGradient; c++) {
	TColor *col = new TColor( Float_t(red[g-1]   + c * (red[g]   - red[g-1])  / nColorsGradient),
				 Float_t(green[g-1] + c * (green[g] - green[g-1])/ nColorsGradient),
				 Float_t(blue[g-1]  + c * (blue[g]  - blue[g-1]) / nColorsGradient),
				  1);
	if(c == 0 && g == 1) st_col = col->GetNumber();
      }
  }

}

int getCol(double y_avg) {
  double y_min = 0, y_max = 8;
  float cv = (n_colors-1)/y_max * y_avg + st_col;
  return do_round(cv);
}

// macro to plot the lifetime dist + fit 
void plotLtScaled()
{
  // get the lifetime distributions
  // prepare binning and histograms for plots
  TH2D **h_d2d = new TH2D*[3];
  string lbl[] = {"SR", "LSB", "RSB"};
  string lbl_m[] = {"3.57 < M(#mu#mu) < 3.81 GeV",
		    "3.4 < M(#mu#mu) < 3.52 GeV",
		    "3.82 < M(#mu#mu) < 4.0 GeV"};
  TFile *fin = new TFile("files/ltStore.root");
  for(int i = 0; i < 3; i++) {
    fin->GetObject(Form("ltH_%s", lbl[i].c_str()), h_d2d[i]);
    h_d2d[i]->SetDirectory(0);
  }
  fin->Close();

  int tbins = h_d2d[0]->GetNbinsX();
  double lowt = h_d2d[0]->GetXaxis()->GetBinLowEdge(1);
  double hit = h_d2d[0]->GetXaxis()->GetBinUpEdge(tbins);
  double wbin = (hit-lowt)/(double)tbins;

  // define aux vals for plotting
  double pr_lim = 0.05;
  double np_lim = 0.1;
  double lowPlot = 0.1;

  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.13);
  c->SetTopMargin(0.015);
  c->SetLogy();

  // Make 1d histos
  TH1D ***h_d1d = new TH1D**[3];
  for(int j = 0; j < 3; j++) {
    h_d1d[j] = new TH1D*[nPtBins];
    for(int i = 0; i < nPtBins; i++) {
      h_d1d[j][i] = h_d2d[j]->ProjectionX(Form("ltH_%s%.0f", lbl[j].c_str(), ptBins[i]), i+1, i+1);
      h_d1d[j][i]->SetTitle(Form("Run 2 %s data c#tau (%.1f < p_{T} < %.1f GeV)", lbl[j].c_str(), ptBins[i], ptBins[i+1]));
    }
  }

  double min_bin = h_d1d[0][0]->GetXaxis()->FindBin(0.1+1e-6);
  init_color();
  double sc_fac[] = {2e-2, 1e-2, 1e-2};
  
  //plot all pt bins together
  for(int j = 0; j < 3; j++) {
    //scale the histos
    for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
      h_d1d[j][i_pt]->Scale(1./h_d1d[j][i_pt]->Integral(min_bin,tbins));
    }
    double plot_max = h_d1d[j][nPtBins-1]->GetBinContent(16)*10;
    double plot_min = plot_max*sc_fac[j];
    
    TH1F *fh = c->DrawFrame(lowPlot, plot_min, hit, plot_max);
    fh->SetXTitle("c#tau (mm)");
    fh->SetYTitle(Form("Events per %.0f #mum (scaled)", wbin*1000.));
    fh->GetYaxis()->SetTitleOffset(1.8);
    fh->GetYaxis()->SetLabelOffset(0.01);

    TLegend *leg = new TLegend(0.7, 0.55, 1., 0.95);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(kWhite,0);

    // start the pT cycle
    for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
      h_d1d[j][i_pt]->SetMarkerStyle(20+5*(i_pt%2));
      h_d1d[j][i_pt]->SetMarkerColor(getCol(i_pt));
      h_d1d[j][i_pt]->SetLineColor(getCol(i_pt));
      h_d1d[j][i_pt]->SetMarkerSize(0.75);
      h_d1d[j][i_pt]->Draw("error same");  

      leg->AddEntry(h_d1d[j][i_pt], Form("[%.0f, %.0f] GeV", ptBins[i_pt], ptBins[i_pt+1]), "pl");
    }

    // draw the state
    TLatex lc;
    lc.SetTextSize(0.04);
    double xp = getPos(lowPlot, hit, 0.1, 0);
    double yp = getPos(plot_min, plot_max, 0.925, 1);
    lc.DrawLatex(xp, yp, Form("#bf{#psi(2S)}"));
    lc.SetTextSize(0.03);
    yp = getPos(plot_min, plot_max, 0.855, 1);
    lc.DrawLatex(xp, yp, Form("#bf{%s}", lbl_m[j].c_str()));

    leg->Draw();


    c->SaveAs(Form("plots/lifetimeS/plot_%s.pdf", lbl[j].c_str()));
    c->Clear();
  }
  c->Destructor();
}
