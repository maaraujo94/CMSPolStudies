// code to get the 2d fine-binned data/mc ratio hist
// plots data, mc and ratio
// saves ratio (normal and |costh|), number of entries in each histo bin

void ratioSave()
{
  double M_q = 1;//3.097; // using J/psi mass

  // open files and read TTrees
  TFile *fin = new TFile("../Store_data_codes/data18_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  
  int dEvt = treeD->GetEntries();

  cout << treeD->GetEntries("lts<2.5") << " PR events after cuts and " << treeD->GetEntries("lts>2.5") << " NP events after cuts" << endl;

  // prepare binning and histograms for plots
  const int nPtBins = 68;
  double ptBins[nPtBins+1];
  for(int i=0 ; i<15; i++) ptBins[i] = (i+25.)/M_q;
  for(int i=0; i<30; i++) ptBins[i+15] = (40.+2.*i)/M_q;;
  for(int i=0; i<14; i++) ptBins[i+45] = (100+2.5*i)/M_q;
  for(int i=0; i<5; i++) ptBins[i+59] = (135. + 5.*i)/M_q;
  for(int i=0; i<5; i++) ptBins[i+64] = (160. + 10.*i)/M_q;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i]*M_q << ",";
  cout << endl;

  TH2D *dataHist = new TH2D("dataH", "Data (PR)", 40, -1., 1., nPtBins, ptBins);
  TH2D *mcHist = new TH2D("mcH", "Data (NP)", 40, -1., 1., nPtBins, ptBins);

  TH2D *dataHist_ab = new TH2D("dataH_ab", "Data (PR)", 20, 0, 1., nPtBins, ptBins);
  TH2D *mcHist_ab = new TH2D("mcH_ab", "Data (NP)", 20, 0, 1., nPtBins, ptBins);
  
  // definitions to store data and MC events
  Double_t data_cos, data_pt, lts;

  treeD->SetBranchAddress("costh", &data_cos);
  treeD->SetBranchAddress("JpsiPt", &data_pt);
  treeD->SetBranchAddress("lts", &lts);
  
  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(lts < 2.5) {
	dataHist->Fill(data_cos, data_pt/M_q);
	dataHist_ab->Fill(abs(data_cos), data_pt/M_q);
      }
      else {
	mcHist->Fill(data_cos, data_pt/M_q);
	mcHist_ab->Fill(abs(data_cos), data_pt/M_q);
      }
    }
  
  // plot all costh histograms
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetLogz();

  double relMax = 2.5e5;
  double absMax = 2*relMax;
  
  dataHist->SetStats(0);
  dataHist->GetXaxis()->SetTitle("cos#theta_{HX}");
  dataHist->GetYaxis()->SetTitle("p_{T} (GeV)");
  //dataHist->SetMaximum(relMax);
  dataHist->Draw("COLZ");
  c->SaveAs("plots/pr_2d.pdf");
  c->Clear();

  mcHist->SetStats(0);
  mcHist->GetXaxis()->SetTitle("cos#theta_{HX}");
  mcHist->GetYaxis()->SetTitle("p_{T} (GeV)");
  //mcHist->SetMaximum(relMax);
  mcHist->Draw("COLZ");
  c->SaveAs("plots/np_2d.pdf");
  c->Clear();
  
  dataHist_ab->SetStats(0);
  dataHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  dataHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");
  //dataHist_ab->SetMaximum(absMax);
  dataHist_ab->Draw("COLZ");
  c->SaveAs("plots/pr_2d_abs.pdf");
  c->Clear();

  mcHist_ab->SetStats(0);
  mcHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  mcHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");
  //mcHist_ab->SetMaximum(absMax);
  mcHist_ab->Draw("COLZ");
  c->SaveAs("plots/np_2d_abs.pdf");
  c->Clear();

  ofstream f_ent;
  f_ent.open("text_output/nEntries.txt");
 
  f_ent << "[pTmin, pTmax]; [cosmin, cosmax]: PR entries / NP entries" << endl << endl;
  for(int pt = 0; pt < nPtBins; pt++) {
    for(int cos = 0; cos < 40; cos++) {
 	f_ent << "[" << ptBins[pt] << ", " << ptBins[pt+1] << "]; ";
	f_ent << "[" << -1.+cos*0.05 << ", " << -1.+(cos+1.)*0.05 << "]: ";
	f_ent << dataHist->GetBinContent(cos+1, pt+1) << " / " << mcHist->GetBinContent(cos+1, pt+1) << endl;
   }
    f_ent << endl;
  }
  f_ent.close();
  
  // get the ratio histogram for each bin
  TH2D *ratioHist = new TH2D("ratioH", "ratioH", 40, -1., 1., nPtBins, ptBins);

  c->SetLogz(0);
  ratioHist = (TH2D*)dataHist->Clone(Form("ratioHist"));
  ratioHist->SetMaximum();
  ratioHist->Sumw2();
  ratioHist->Divide(mcHist);
  ratioHist->SetTitle("PR/NP");
  ratioHist->Draw("COLZ");
  c->SaveAs("plots/ratio_2d.pdf");
  c->Clear();

  TH2D *ratioHist_ab = new TH2D("ratioH_ab", "ratioH abs", 20, 0, 1., nPtBins, ptBins);

  c->SetLogz(0);
  ratioHist_ab = (TH2D*)dataHist_ab->Clone(Form("ratioHist_ab"));
  ratioHist_ab->SetMaximum();
  ratioHist_ab->Sumw2();
  ratioHist_ab->Divide(mcHist_ab);
  ratioHist_ab->SetTitle("PR/NP");
  ratioHist_ab->Draw("COLZ");
  c->SaveAs("plots/ratio_2d_abs.pdf");
  c->Clear();
  
  TFile *outfile = new TFile("files/ratioHist.root", "recreate");
  ratioHist->Write();
  ratioHist_ab->Write();
  outfile->Close();
  
}
