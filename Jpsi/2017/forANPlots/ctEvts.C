void ctEvts()
{
  int pt_n = 4;
  
  double pt_min[] = {25, 45, 50, 70};
  double pt_max[] = {45, 50, 70, 120};
  
  double m_min[] = {2.92, 3.0, 3.21};
  double m_max[] = {2.95, 3.2, 3.28};

  double lt_min[] = {-0.005, 0.01};
  double lt_max[] = {0.005,  0.05};
  
  double n_PRSR[pt_n+1];
  double n_NPSR[pt_n+1];
  double n_PLSB[pt_n+1];
  double n_PRSB[pt_n+1];
  double n_NLSB[pt_n+1];
  double n_NRSB[pt_n+1];
  double n_MC[pt_n+1];

  // open files and read TTrees
  TFile *finD = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/data17_cos.root");
  TTree *tree = (TTree*)finD->Get("data_cos");

  for(int i_pt = 0; i_pt < pt_n; i_pt++) { // cycle in pt region
    n_PRSR[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[0], lt_max[0], m_min[1], m_max[1]));
    n_NPSR[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[1], lt_max[1], m_min[1], m_max[1]));
    n_PLSB[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[0], lt_max[0], m_min[0], m_max[0]));
    n_PRSB[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[0], lt_max[0], m_min[2], m_max[2]));
    n_NLSB[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[1], lt_max[1], m_min[0], m_max[0]));
    n_NRSB[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[1], lt_max[1], m_min[2], m_max[2]));
  }
  n_PRSR[pt_n] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[pt_n-1], lt_min[0], lt_max[0], m_min[1], m_max[1]));
  n_NPSR[pt_n] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[pt_n-1], lt_min[1], lt_max[1], m_min[1], m_max[1]));
  n_PLSB[pt_n] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[pt_n-1], lt_min[0], lt_max[0], m_min[0], m_max[0]));
  n_PRSB[pt_n] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[pt_n-1], lt_min[0], lt_max[0], m_min[2], m_max[2]));
  n_NLSB[pt_n] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[pt_n-1], lt_min[1], lt_max[1], m_min[0], m_max[0]));
  n_NRSB[pt_n] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[pt_n-1], lt_min[1], lt_max[1], m_min[2], m_max[2]));

  finD->Close();

  TFile *fin1 = new TFile("../../Store_data_codes/MC17_cos.root");
  TTree *treeM1 = (TTree*)fin1->Get("MC_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/MCm17_cos.root");
  TTree *treeM2 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../../Store_data_codes/MCh17_cos.root");
  TTree *treeM3 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("../../Store_data_codes/MCvh17_cos.root");
  TTree *treeM4 = (TTree*)fin4->Get("MC_cos");

  n_MC[0] = treeM1->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[0], lt_min[0], lt_max[0], m_min[1], m_max[1]));
  n_MC[1] = treeM2->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[1], pt_max[1], lt_min[0], lt_max[0], m_min[1], m_max[1]));
  n_MC[2] = treeM3->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[2], pt_max[2], lt_min[0], lt_max[0], m_min[1], m_max[1]));
  n_MC[3] = treeM4->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[3], pt_max[3], lt_min[0], lt_max[0], m_min[1], m_max[1]));

  n_MC[4] = n_MC[0] + n_MC[1] + n_MC[2] + n_MC[3];

  fin1->Close();
  fin2->Close();
  fin3->Close();
  fin4->Close();

  // output results into a .tex table
  ofstream ftex;
  ftex.open(Form("text_output/data_mc_evts_2017.tex"));
  ftex << "\\begin{tabular}{cc|cccc|c}\n";
  ftex << "\\hline\n";
  ftex << Form("\\multicolumn{2}{c}{2017} & $[%.0f, %.0f]$ GeV & $[%.0f, %.0f]$ GeV & $[%.0f, %.0f]$ GeV & $[%.0f, %.0f]$ GeV & $[%.0f, %.0f]$ GeV \\\\\n", pt_min[0], pt_max[0], pt_min[1], pt_max[1], pt_min[2], pt_max[2], pt_min[3], pt_max[3], pt_min[0], pt_max[3]);
  ftex << "\\hline\n";

  // data output
  ftex << Form("\\multirow{6}{*}{\\rotatebox[origin=c]{90}{Data}} & Peak & %.3f M & %.3f M & %.3f M & %.3f M & %.3f M \\\\\n", n_PRSR[0]/1e6, n_PRSR[1]/1e6, n_PRSR[2]/1e6, n_PRSR[3]/1e6, n_PRSR[4]/1e6);
  ftex << Form("& NP & %.3f M & %.3f M & %.3f M & %.3f M & %.3f M \\\\\n", n_NPSR[0]/1e6, n_NPSR[1]/1e6, n_NPSR[2]/1e6, n_NPSR[3]/1e6, n_NPSR[4]/1e6);
  ftex << Form("& PR LSB & %.1f k & %.1f k & %.1f k & %.1f k & %.1f k \\\\\n", n_PLSB[0]/1e3, n_PLSB[1]/1e3, n_PLSB[2]/1e3, n_PLSB[3]/1e3, n_PLSB[4]/1e3);
  ftex << Form("& PR RSB & %.1f k & %.1f k & %.1f k & %.1f k & %.1f k \\\\\n", n_PRSB[0]/1e3, n_PRSB[1]/1e3, n_PRSB[2]/1e3, n_PRSB[3]/1e3, n_PRSB[4]/1e3);
  ftex << Form("& NP LSB & %.1f k & %.1f k & %.1f k & %.1f k & %.1f k \\\\\n", n_NLSB[0]/1e3, n_NLSB[1]/1e3, n_NLSB[2]/1e3, n_NLSB[3]/1e3, n_NLSB[4]/1e3);
  ftex << Form("& NP RSB & %.1f k & %.1f k & %.1f k & %.1f k & %.1f k \\\\\n", n_NRSB[0]/1e3, n_NRSB[1]/1e3, n_NRSB[2]/1e3, n_NRSB[3]/1e3, n_NRSB[4]/1e3);
  ftex << "\\hline\n";

  // MC output
  ftex << Form("MC & (only Peak) & %.3f M & %.3f M & %.3f M & %.3f M & %.3f M \\\\\n", n_MC[0]/1e6, n_MC[1]/1e6, n_MC[2]/1e6, n_MC[3]/1e6, n_MC[4]/1e6);
  ftex << "\\hline\n";
  
  ftex << "\\end{tabular}\n";
  ftex.close();

}
