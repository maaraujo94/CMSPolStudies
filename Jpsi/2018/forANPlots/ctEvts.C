void ctEvts()
{
  double pt_min[] = {25, 46, 66};
  double pt_max[] = {46, 66, 120};
  
  double m_min[] = {2.92, 3.0, 3.21};
  double m_max[] = {2.95, 3.2, 3.28};

  double lt_min[] = {-0.01, 0.014};
  double lt_max[] = {0.01,  0.05};
  
  double n_PRSR[4];
  double n_NP[4];
  double n_LSB[4];
  double n_RSB[4];

  // open files and read TTrees
  TFile *finD = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/data18_cos.root");
  TTree *tree = (TTree*)finD->Get("data_cos");

  for(int i_pt = 0; i_pt < 3; i_pt++) { // cycle in pt region
    n_PRSR[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[0], lt_max[0], m_min[1], m_max[1]));
    n_NP[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[1], lt_max[1], m_min[1], m_max[1]));
    n_LSB[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[0], lt_max[0], m_min[0], m_max[0]));
    n_RSB[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[0], lt_max[0], m_min[2], m_max[2]));
  }
  n_PRSR[3] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[2], lt_min[0], lt_max[0], m_min[1], m_max[1]));
  n_NP[3] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[2], lt_min[1], lt_max[1], m_min[1], m_max[1]));
  n_LSB[3] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[2], lt_min[0], lt_max[0], m_min[0], m_max[0]));
  n_RSB[3] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[2], lt_min[0], lt_max[0], m_min[2], m_max[2]));

  finD->Close();

  cout << "PR SR results: " << endl;
  for(int i = 0; i < 3; i++) {
    cout << n_PRSR[i];
    if(i < 2) cout << " + ";
    else cout << " = ";
  }
  cout << n_PRSR[3] << " aka " << n_PRSR[0]+n_PRSR[1]+n_PRSR[2] << endl << endl;

  cout << "NP results:" << endl;
  for(int i = 0; i < 3; i++) {
    cout << n_NP[i];
    if(i < 2) cout << " + ";
    else cout << " = ";
  }
  cout << n_NP[3] << " aka " << n_NP[0]+n_NP[1]+n_NP[2] << endl << endl;

  cout << "LSB results:" << endl;
  for(int i = 0; i < 3; i++) {
    cout << n_LSB[i];
    if(i < 2) cout << " + ";
    else cout << " = ";
  }
  cout << n_LSB[3] << " aka " << n_LSB[0]+n_LSB[1]+n_LSB[2] << endl << endl;

  cout << "RSB results:" << endl;
  for(int i = 0; i < 3; i++) {
    cout << n_RSB[i];
    if(i < 2) cout << " + ";
    else cout << " = ";
  }
  cout << n_RSB[3] << " aka " << n_RSB[0]+n_RSB[1]+n_RSB[2] << endl << endl;
}
