void getBinEvts()
{
  TH2F **h_PR = new TH2F*[4];
  
  TFile *fin7 = new TFile("../2017/PR_fit/files/histoStore.root");
  h_PR[0] = (TH2F*)fin7->Get("PRH");
  h_PR[0]->SetDirectory(0);
  fin7->Close();

  TFile *fin8 = new TFile("../2018/PR_fit/files/histoStore.root");
  h_PR[1] = (TH2F*)fin8->Get("PRH");
  h_PR[1]->SetDirectory(0);
  fin8->Close();

  TFile *finS = new TFile("../Simult/PR_fit/files/histoStore.root");
  h_PR[2] = (TH2F*)finS->Get("PRH");
  h_PR[2]->SetDirectory(0);
  finS->Close();

  TFile *fin6 = new TFile("../2016/PR_fit/files/histoStore.root");
  h_PR[3] = (TH2F*)fin6->Get("PRH");
  h_PR[3]->SetDirectory(0);
  fin6->Close();

  int nBinsY = h_PR[0]->GetNbinsY();
  const double *yBins = h_PR[0]->GetYaxis()->GetXbins()->GetArray();

  for(int i = 0; i < nBinsY; i++) {
    cout << "[" << yBins[i] << " , " << yBins[i+1] << "]: "; 
    for(int j = 0; j < 4; j++)
      cout << h_PR[j]->Integral(1,20,i+1,i+1) << "; ";
    cout << endl;
  }
}
