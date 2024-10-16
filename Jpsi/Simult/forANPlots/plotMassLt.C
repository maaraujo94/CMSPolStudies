// macro to plot a mass:lifetime 2d histo for all data
void plotMassLt()
{
  TH2D *hist = new TH2D("name", "Full Data", 36, 2.92, 3.28, 180, -100, 800);

  // open files and read TTrees
  TFile *finD = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/dataS_cos.root");
  TTree *treeD = (TTree*)finD->Get("data_cos");

  int dEvt = treeD->GetEntries();

  // definitions to store data and MC events
  Double_t data_pt, data_lt, data_m, data_y;
  
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(data_pt > 25 && data_pt < 120 && abs(data_y) < 1.2)
	hist->Fill(data_m, data_lt*1e4);
    }
  finD->Close();

  double m_min[] = {2.92, 3.0, 3.21};
  double m_max[] = {2.95, 3.2, 3.28};
  double lt_min[] = {-50, 100};
  double lt_max[] = {50, 800};

  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLeftMargin(0.11);
  //c->SetRightMargin(0.03);
  c->SetRightMargin(0.11);
  c->SetTopMargin(0.02);
  c->SetLogz();
 
  hist->SetStats(0);
  hist->GetXaxis()->SetTitle("#it{M} (GeV)");
  hist->GetYaxis()->SetTitle("c#tau (#mum)");
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetLabelOffset(0.01);
  hist->GetXaxis()->SetLabelOffset(0.01);
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetXaxis()->CenterTitle(true);
  hist->SetTitle("");
  hist->Draw("COLZ");

    for(int i = 0; i <3; i++){
    for(int j = 0; j< 2; j++){
      TLine *l1 = new TLine(m_min[i], lt_min[j], m_min[i], lt_max[j]);
      l1->SetLineColor(kRed);
      l1->SetLineWidth(2);
      l1->Draw();
      TLine *l2 = new TLine(m_max[i], lt_min[j], m_max[i], lt_max[j]);
      l2->SetLineColor(kRed);
      l2->SetLineWidth(2);
      l2->Draw();
      TLine *l3 = new TLine(m_min[i], lt_min[j], m_max[i], lt_min[j]);
      l3->SetLineColor(kRed);
      l3->SetLineWidth(2);
      l3->Draw();
      TLine *l4 = new TLine(m_min[i], lt_max[j], m_max[i], lt_max[j]);
      l4->SetLineColor(kRed);
      l4->SetLineWidth(2);
      l4->Draw();
    }
  }

  TLatex dl;
  dl.SetTextSize(0.035);
  dl.SetTextColor(kRed);
  dl.DrawLatex(m_min[1]*1.002, lt_max[1]*0.5, "NPS");
  dl.DrawLatex(m_min[1]*1.002, lt_max[0]*0.3, "PRS");

  c->SaveAs("plots/2dMaps/massLtMap.pdf");
  c->Clear();
  c->Destructor();
}
