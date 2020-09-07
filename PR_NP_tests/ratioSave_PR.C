// code to get the 2d fine-binned PR/PR (diff cuts) ratio hist
// plots PR 1, 2 and ratio
// saves ratio (normal and |costh|), number of entries in each histo bin

void ratioSave_PR()
{
  double M_q = 3.097; // using J/psi mass

  // open files and read TTrees
  TFile *fin = new TFile("../Store_data_codes/data_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  
  int dEvt = treeD->GetEntries();

  // prepare binning and histograms for plots
  const int nPtBins = 29;
  double ptBins[nPtBins+1];
  for(int i=0 ; i<9; i++) ptBins[i] = (i+12.)/M_q;
  for(int i=0; i<15; i++) ptBins[i+9] = (22.+2.*i)/M_q;;
  for(int i=0; i<4; i++) ptBins[i+24] = (52.5+2.5*i)/M_q;
  for(int i=0; i<2; i++) ptBins[i+28] = (65. + 5.*i)/M_q;
  for(int i=0; i<30; i++) cout << ptBins[i]*M_q << ",";
  cout << endl;
  
  TH2D *PRHist1 = new TH2D("PR1", "Data (PR-1)", 36, -0.9, 0.9, nPtBins, ptBins);
  TH2D *PRHist2 = new TH2D("PR2", "Data (PR-2)", 36, -0.9, 0.9, nPtBins, ptBins);

  TH2D *PRHist1_ab = new TH2D("PR1_ab", "Data (PR-1)", 18, 0, 0.9, nPtBins, ptBins);
  TH2D *PRHist2_ab = new TH2D("PR2_ab", "Data (PR-2)", 18, 0, 0.9, nPtBins, ptBins);
  
  // definitions to store data and MC events
  Double_t data_cos, data_pt, lts;

  treeD->SetBranchAddress("costh", &data_cos);
  treeD->SetBranchAddress("JpsiPt", &data_pt);
  treeD->SetBranchAddress("lts", &lts);
  
  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(lts < 1) {
	PRHist1->Fill(data_cos, data_pt/M_q);
	PRHist1_ab->Fill(abs(data_cos), data_pt/M_q);
      }
      else {
	PRHist2->Fill(data_cos, data_pt/M_q);
	PRHist2_ab->Fill(abs(data_cos), data_pt/M_q);
      }
    }

  cout << PRHist1->GetEntries() << " entries for PR-1 and " << PRHist2->GetEntries() << " entries for PR-2" << endl;
  
  // plot all costh histograms
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLogz();

  double relMax = 2.5e5;
  double absMax = 2*relMax;
  
  PRHist1->SetStats(0);
  PRHist1->GetXaxis()->SetTitle("cos#theta_{HX}");
  PRHist1->GetYaxis()->SetTitle("p_{T}/M ");
  PRHist1->SetMaximum(relMax);
  PRHist1->Draw("COLZ");
  c->SaveAs("plots/pr_1_2d.pdf");
  c->Clear();

  PRHist2->SetStats(0);
  PRHist2->GetXaxis()->SetTitle("cos#theta_{HX}");
  PRHist2->GetYaxis()->SetTitle("p_{T}/M");
  PRHist2->SetMaximum(relMax);
  PRHist2->Draw("COLZ");
  c->SaveAs("plots/pr_2_2d.pdf");
  c->Clear();
  
  PRHist1_ab->SetStats(0);
  PRHist1_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  PRHist1_ab->GetYaxis()->SetTitle("p_{T}/M");
  PRHist1_ab->SetMaximum(absMax);
  PRHist1_ab->Draw("COLZ");
  c->SaveAs("plots/pr_1_2d_abs.pdf");
  c->Clear();

  PRHist2_ab->SetStats(0);
  PRHist2_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  PRHist2_ab->GetYaxis()->SetTitle("p_{T}/M");
  PRHist2_ab->SetMaximum(absMax);
  PRHist2_ab->Draw("COLZ");
  c->SaveAs("plots/pr_2_2d_abs.pdf");
  c->Clear();

  ofstream f_ent;
  f_ent.open("text_output/nEntries.txt");
 
  f_ent << "[pTmin, pTmax]; [cosmin, cosmax]: PR(1) entries / PR(2) entries" << endl << endl;
  for(int pt = 0; pt < nPtBins; pt++) {
    for(int cos = 0; cos < 36; cos++) {
 	f_ent << "[" << ptBins[pt] << ", " << ptBins[pt+1] << "]; ";
	f_ent << "[" << -0.9+cos*0.05 << ", " << -0.9+(cos+1.)*0.05 << "]: ";
	f_ent << PRHist1->GetBinContent(cos+1, pt+1) << " / " << PRHist2->GetBinContent(cos+1, pt+1) << endl;
   }
    f_ent << endl;
  }
  f_ent.close();
  
  // get the ratio histogram for each bin
  TH2D *ratioHist = new TH2D("ratioH", "ratioH", 36, -0.9, 0.9, nPtBins, ptBins);

  c->SetLogz(0);
  ratioHist = (TH2D*)PRHist1->Clone(Form("ratioHist"));
  ratioHist->SetMaximum();
  ratioHist->Sumw2();
  ratioHist->Divide(PRHist2);
  ratioHist->SetTitle("PR(1)/PR(2)");
  ratioHist->Draw("COLZ");
  c->SaveAs("plots/ratio_2d.pdf");
  c->Clear();

  TH2D *ratioHist_ab = new TH2D("ratioH_ab", "ratioH abs", 18, 0, 0.9, nPtBins, ptBins);

  c->SetLogz(0);
  ratioHist_ab = (TH2D*)PRHist1_ab->Clone(Form("ratioHist_ab"));
  ratioHist_ab->SetMaximum();
  ratioHist_ab->Sumw2();
  ratioHist_ab->Divide(PRHist2_ab);
  ratioHist_ab->SetTitle("PR(1)/PR(2)");
  ratioHist_ab->Draw("COLZ");
  c->SaveAs("plots/ratio_2d_abs.pdf");
  c->Clear();
  
  TFile *outfile = new TFile("files/ratioHist.root", "recreate");
  ratioHist->Write();
  ratioHist_ab->Write();
  outfile->Close();

}
