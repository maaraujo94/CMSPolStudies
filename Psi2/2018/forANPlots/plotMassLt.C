// macro to plot a mass:lifetime 2d histo for all data
void plotMassLt()
{
  TH2D *hist = new TH2D("name", "2018 Data", 36, 3.4, 4.0, 120, -0.01, 0.05);

  // open files and read TTrees
  TFile *finD = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/data18_cos.root");
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
      if(data_pt > 20 && data_pt < 100 && abs(data_y) < 1.2)
	hist->Fill(data_m, data_lt);
    }
  finD->Close();

  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLeftMargin(0.11);
  c->SetLogz();
 
  hist->SetStats(0);
  hist->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
  hist->GetYaxis()->SetTitle("c#tau (cm)");
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->Draw("COL");
  c->SaveAs("plots/2dMaps/massLtMap.pdf");
  c->Clear();

}
