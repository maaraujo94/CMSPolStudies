// code to get the 2d fine-binned data/mc ratio hist
// plots data, mc and ratio
// saves ratio (normal and |costh|), number of entries in each histo bin

void ratioSave()
{
  // open files and read TTrees
  TFile *fin = new TFile("../../Store_data_codes/2018/data18_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/2018/MC18_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../../Store_data_codes/2018/MC18_hpt_cos.root");
  TTree *treeM2 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("../../Store_data_codes/2018/MC18_vhpt_cos.root");
  TTree *treeM3 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM2->GetEntries();

  cout << "Running first sample" << endl;
  cout << treeD->GetEntries("lts<2.5 && JpsiPt > 25 && JpsiPt < 46") << " data events after cuts and " << treeM1->GetEntries("JpsiPt<46") << " MC events after cuts" << endl;

  // prepare binning and histograms for plots 
  const int nPtBins = 18;
  double ptBins[nPtBins+1];
  for(int i=0; i<15; i++) ptBins[i] = i+25.;
  for(int i=0; i<4; i++) ptBins[i+15] = 40.+2.*i;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;
 
  TH2D *dataHist = new TH2D("dataH", "Data (PR)", 40, -1., 1., nPtBins, ptBins);
  TH2D *mcHist = new TH2D("mcH", "MC", 40, -1., 1., nPtBins, ptBins);

  TH2D *dataHist_ab = new TH2D("dataH_ab", "Data (PR)", 20, 0, 1., nPtBins, ptBins);
  TH2D *mcHist_ab = new TH2D("mcH_ab", "MC", 20, 0, 1., nPtBins, ptBins);
  
  // definitions to store data and MC events
  Double_t data_cos, data_pt, lts;
  Double_t mc_cos, mc_pt;
 
  treeD->SetBranchAddress("costh", &data_cos);
  treeD->SetBranchAddress("JpsiPt", &data_pt);
  treeD->SetBranchAddress("lts", &lts);
  
  treeM1->SetBranchAddress("costh", &mc_cos);
  treeM1->SetBranchAddress("JpsiPt", &mc_pt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(lts < 2.5) {
	dataHist->Fill(data_cos, data_pt);
	dataHist_ab->Fill(abs(data_cos), data_pt);
      }
    }

  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      mcHist->Fill(mc_cos, mc_pt);
      mcHist_ab->Fill(abs(mc_cos), mc_pt);
    }

  // plot all costh histograms
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.11);
  c->SetLogz();
  
  dataHist->SetStats(0);
  dataHist->GetXaxis()->SetTitle("cos#theta_{HX}");
  dataHist->GetYaxis()->SetTitle("p_{T} (GeV)");
  dataHist->Draw("COLZ");
  c->SaveAs("plots/data_2d_1.pdf");
  c->Clear();

  mcHist->SetStats(0);
  mcHist->GetXaxis()->SetTitle("cos#theta_{HX}");
  mcHist->GetYaxis()->SetTitle("p_{T} (GeV)");
  mcHist->Draw("COLZ");
  c->SaveAs("plots/mc_2d_1.pdf");
  c->Clear();

  dataHist_ab->SetStats(0);
  dataHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  dataHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");
  dataHist_ab->Draw("COLZ");
  c->SaveAs("plots/data_2d_abs_1.pdf");
  c->Clear();

  mcHist_ab->SetStats(0);
  mcHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  mcHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");
  mcHist_ab->Draw("COLZ");
  c->SaveAs("plots/mc_2d_abs_1.pdf");
  c->Clear();
  
  ofstream f_ent;
  f_ent.open("text_output/nEntries_1.txt");
 
  f_ent << "[pTmin, pTmax]; [cosmin, cosmax]: data entries / mc entries" << endl << endl;
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
  ratioHist->Sumw2();
  ratioHist->Divide(mcHist);
  ratioHist->SetTitle("PR/MC");
  ratioHist->Draw("COLZ");
  c->SaveAs("plots/ratio_2d_1.pdf");
  c->Clear();

  TH2D *ratioHist_ab = new TH2D("ratioH_ab", "ratioH abs", 20, 0, 1., nPtBins, ptBins);

  c->SetLogz(0);
  ratioHist_ab = (TH2D*)dataHist_ab->Clone(Form("ratioHist_ab"));
  ratioHist_ab->Sumw2();
  ratioHist_ab->Divide(mcHist_ab);
  ratioHist_ab->SetTitle("PR/MC");
  ratioHist_ab->Draw("COLZ");
  c->SaveAs("plots/ratio_2d_abs_1.pdf");
  c->Clear();

    cout << "Running second sample" << endl;
  cout << treeD->GetEntries("lts<2.5 && JpsiPt > 46 && JpsiPt < 66") << " data events after cuts and " << treeM2->GetEntries("JpsiPt<66") << " MC events after cuts" << endl;

  // prepare binning and histograms for plots 
  const int nPtBins_hpt = 10;
  double ptBins_hpt[nPtBins_hpt+1];
  for(int i=0; i<11; i++) ptBins_hpt[i] = 46.+2.*i;
  for(int i=0; i<nPtBins_hpt+1; i++) cout << ptBins_hpt[i] << ",";
  cout << endl;
  
  TH2D *dataHist_hpt = new TH2D("dataH_hpt", "Data (PR)", 40, -1.0, 1.0, nPtBins_hpt, ptBins_hpt);
  TH2D *mcHist_hpt = new TH2D("mcH_hpt", "MC", 40, -1.0, 1.0, nPtBins_hpt, ptBins_hpt);

  TH2D *dataHist_ab_hpt = new TH2D("dataH_ab_hpt", "Data (PR)", 20 , 0, 1.0, nPtBins_hpt, ptBins_hpt);
  TH2D *mcHist_ab_hpt = new TH2D("mcH_ab_hpt", "MC", 20 , 0, 1.0, nPtBins_hpt, ptBins_hpt);
  
  // definitions to store data and MC events
  treeM2->SetBranchAddress("costh", &mc_cos);
  treeM2->SetBranchAddress("JpsiPt", &mc_pt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(lts < 2.5) {
	dataHist_hpt->Fill(data_cos, data_pt);
	dataHist_ab_hpt->Fill(abs(data_cos), data_pt);
      }
    }

  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
      mcHist_hpt->Fill(mc_cos, mc_pt);
      mcHist_ab_hpt->Fill(abs(mc_cos), mc_pt);
    }
  
  // plot all costh histograms
  c->SetRightMargin(0.11);
  c->SetLogz();
  
  dataHist_hpt->SetStats(0);
  dataHist_hpt->GetXaxis()->SetTitle("cos#theta_{HX}");
  dataHist_hpt->GetYaxis()->SetTitle("p_{T} (GeV)");
  dataHist_hpt->Draw("COLZ");
  c->SaveAs("plots/data_2d_2.pdf");
  c->Clear();

  mcHist_hpt->SetStats(0);
  mcHist_hpt->GetXaxis()->SetTitle("cos#theta_{HX}");
  mcHist_hpt->GetYaxis()->SetTitle("p_{T} (GeV)");
  mcHist_hpt->Draw("COLZ");
  c->SaveAs("plots/mc_2d_2.pdf");
  c->Clear();

  dataHist_ab_hpt->SetStats(0);
  dataHist_ab_hpt->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  dataHist_ab_hpt->GetYaxis()->SetTitle("p_{T} (GeV)");
  dataHist_ab_hpt->Draw("COLZ");
  c->SaveAs("plots/data_2d_abs_2.pdf");
  c->Clear();

  mcHist_ab_hpt->SetStats(0);
  mcHist_ab_hpt->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  mcHist_ab_hpt->GetYaxis()->SetTitle("p_{T} (GeV)");
  mcHist_ab_hpt->Draw("COLZ");
  c->SaveAs("plots/mc_2d_abs_2.pdf");
  c->Clear();
  
  ofstream f_ent_hpt;
  f_ent_hpt.open("text_output/nEntries_2.txt");
 
  f_ent_hpt << "[pTmin, pTmax]; [cosmin, cosmax]: data entries / mc entries" << endl << endl;
  for(int pt = 0; pt < nPtBins_hpt; pt++) {
    for(int cos = 0; cos < 40; cos++) {
      f_ent_hpt << "[" << ptBins_hpt[pt] << ", " << ptBins_hpt[pt+1] << "]; ";
      f_ent_hpt << "[" << -1.+cos*0.05 << ", " << -1.+(cos+1.)*0.05 << "]: ";
      f_ent_hpt << dataHist_hpt->GetBinContent(cos+1, pt+1) << " / " << mcHist_hpt->GetBinContent(cos+1, pt+1) << endl;
    }
    f_ent_hpt << endl;
  }
  f_ent_hpt.close();
  
  // get the ratio histogram for each bin
  TH2D *ratioHist_hpt = new TH2D("ratioH_hpt", "ratioH", 40, -1.0, 1.0, nPtBins_hpt, ptBins_hpt);

  c->SetLogz(0);
  ratioHist_hpt = (TH2D*)dataHist_hpt->Clone(Form("ratioHist_hpt"));
  ratioHist_hpt->Sumw2();
  ratioHist_hpt->Divide(mcHist_hpt);
  ratioHist_hpt->SetTitle("PR/MC");
  ratioHist_hpt->Draw("COLZ");
  c->SaveAs("plots/ratio_2d_2.pdf");
  c->Clear();

  TH2D *ratioHist_ab_hpt = new TH2D("ratioH_ab_hpt", "ratioH abs", 20 , 0, 1.0, nPtBins_hpt, ptBins_hpt);

  c->SetLogz(0);
  ratioHist_ab_hpt = (TH2D*)dataHist_ab_hpt->Clone(Form("ratioHist_ab_hpt"));
  ratioHist_ab_hpt->Sumw2();
  ratioHist_ab_hpt->Divide(mcHist_ab_hpt);
  ratioHist_ab_hpt->SetTitle("PR/MC");
  ratioHist_ab_hpt->Draw("COLZ");
  c->SaveAs("plots/ratio_2d_abs_2.pdf");
  c->Clear();

  cout << "Running third sample" << endl;
  cout << treeD->GetEntries("lts<2.5 && JpsiPt > 66 && JpsiPt < 100") << " data events after cuts and " << treeM3->GetEntries("JpsiPt<100") << " MC events after cuts" << endl;

  // prepare binning and histograms for plots 
  const int nPtBins_vhpt = 17;
  double ptBins_vhpt[nPtBins_vhpt+1];
  for(int i=0; i<18; i++) ptBins_vhpt[i] = 66.+2.*i;
  for(int i=0; i<nPtBins_vhpt+1; i++) cout << ptBins_vhpt[i] << ",";
  cout << endl;
  
  TH2D *dataHist_vhpt = new TH2D("dataH_vhpt", "Data (PR)", 40, -1.0, 1.0, nPtBins_vhpt, ptBins_vhpt);
  TH2D *mcHist_vhpt = new TH2D("mcH_vhpt", "MC", 40, -1.0, 1.0, nPtBins_vhpt, ptBins_vhpt);

  TH2D *dataHist_ab_vhpt = new TH2D("dataH_ab_vhpt", "Data (PR)", 20 , 0, 1.0, nPtBins_vhpt, ptBins_vhpt);
  TH2D *mcHist_ab_vhpt = new TH2D("mcH_ab_vhpt", "MC", 20 , 0, 1.0, nPtBins_vhpt, ptBins_vhpt);
  
  // definitions to store data and MC events
  treeM3->SetBranchAddress("costh", &mc_cos);
  treeM3->SetBranchAddress("JpsiPt", &mc_pt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(lts < 2.5) {
	dataHist_vhpt->Fill(data_cos, data_pt);
	dataHist_ab_vhpt->Fill(abs(data_cos), data_pt);
      }
    }

  for(int i = 0; i < m2Evt; i++)
    {
      treeM3->GetEntry(i);
      mcHist_vhpt->Fill(mc_cos, mc_pt);
      mcHist_ab_vhpt->Fill(abs(mc_cos), mc_pt);
    }
  
  // plot all costh histograms
  c->SetRightMargin(0.11);
  c->SetLogz();
  
  dataHist_vhpt->SetStats(0);
  dataHist_vhpt->GetXaxis()->SetTitle("cos#theta_{HX}");
  dataHist_vhpt->GetYaxis()->SetTitle("p_{T} (GeV)");
  dataHist_vhpt->Draw("COLZ");
  c->SaveAs("plots/data_2d_3.pdf");
  c->Clear();

  mcHist_vhpt->SetStats(0);
  mcHist_vhpt->GetXaxis()->SetTitle("cos#theta_{HX}");
  mcHist_vhpt->GetYaxis()->SetTitle("p_{T} (GeV)");
  mcHist_vhpt->Draw("COLZ");
  c->SaveAs("plots/mc_2d_3.pdf");
  c->Clear();

  dataHist_ab_vhpt->SetStats(0);
  dataHist_ab_vhpt->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  dataHist_ab_vhpt->GetYaxis()->SetTitle("p_{T} (GeV)");
  dataHist_ab_vhpt->Draw("COLZ");
  c->SaveAs("plots/data_2d_abs_3.pdf");
  c->Clear();

  mcHist_ab_vhpt->SetStats(0);
  mcHist_ab_vhpt->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  mcHist_ab_vhpt->GetYaxis()->SetTitle("p_{T} (GeV)");
  mcHist_ab_vhpt->Draw("COLZ");
  c->SaveAs("plots/mc_2d_abs_3.pdf");
  c->Clear();
  
  ofstream f_ent_vhpt;
  f_ent_vhpt.open("text_output/nEntries_3.txt");
 
  f_ent_vhpt << "[pTmin, pTmax]; [cosmin, cosmax]: data entries / mc entries" << endl << endl;
  for(int pt = 0; pt < nPtBins_vhpt; pt++) {
    for(int cos = 0; cos < 40; cos++) {
      f_ent_vhpt << "[" << ptBins_vhpt[pt] << ", " << ptBins_vhpt[pt+1] << "]; ";
      f_ent_vhpt << "[" << -1.+cos*0.05 << ", " << -1.+(cos+1.)*0.05 << "]: ";
      f_ent_vhpt << dataHist_vhpt->GetBinContent(cos+1, pt+1) << " / " << mcHist_vhpt->GetBinContent(cos+1, pt+1) << endl;
    }
    f_ent_vhpt << endl;
  }
  f_ent_vhpt.close();
  
  // get the ratio histogram for each bin
  TH2D *ratioHist_vhpt = new TH2D("ratioH_vhpt", "ratioH", 40, -1.0, 1.0, nPtBins_vhpt, ptBins_vhpt);

  c->SetLogz(0);
  ratioHist_vhpt = (TH2D*)dataHist_vhpt->Clone(Form("ratioHist_vhpt"));
  ratioHist_vhpt->Sumw2();
  ratioHist_vhpt->Divide(mcHist_vhpt);
  ratioHist_vhpt->SetTitle("PR/MC");
  ratioHist_vhpt->Draw("COLZ");
  c->SaveAs("plots/ratio_2d_3.pdf");
  c->Clear();

  TH2D *ratioHist_ab_vhpt = new TH2D("ratioH_ab_vhpt", "ratioH abs", 20 , 0, 1.0, nPtBins_vhpt, ptBins_vhpt);

  c->SetLogz(0);
  ratioHist_ab_vhpt = (TH2D*)dataHist_ab_vhpt->Clone(Form("ratioHist_ab_vhpt"));
  ratioHist_ab_vhpt->Sumw2();
  ratioHist_ab_vhpt->Divide(mcHist_ab_vhpt);
  ratioHist_ab_vhpt->SetTitle("PR/MC");
  ratioHist_ab_vhpt->Draw("COLZ");
  c->SaveAs("plots/ratio_2d_abs_3.pdf");
  c->Clear();

  TFile *outfile = new TFile("files/ratioHist.root", "recreate");
  ratioHist->Write();
  ratioHist_ab->Write();
  ratioHist_hpt->Write();
  ratioHist_ab_hpt->Write();
  ratioHist_vhpt->Write();
  ratioHist_ab_vhpt->Write();
  outfile->Close();

}
