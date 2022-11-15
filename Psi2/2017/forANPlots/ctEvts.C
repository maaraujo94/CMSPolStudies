void ctEvts()
{
  int pt_n = 1;
  
  double pt_min[] = {25};
  double pt_max[] = {100};
  
  double m_min[] = {3.4, 3.57, 3.82};
  double m_max[] = {3.52, 3.81, 4.0};

  double lt_min[] = {-0.005, 0.01};
  double lt_max[] = {0.005,  0.05};
  
  double n_PRSR[pt_n];
  double n_NPSR[pt_n];
  double n_PLSB[pt_n];
  double n_PRSB[pt_n];
  double n_NLSB[pt_n];
  double n_NRSB[pt_n];
  double n_MC[pt_n];

  // open files and read TTrees
  TFile *finD = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/data17_cos.root");
  TTree *tree = (TTree*)finD->Get("data_cos");

  for(int i_pt = 0; i_pt < pt_n; i_pt++) { // cycle in pt region
    n_PRSR[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[0], lt_max[0], m_min[1], m_max[1]));
    n_NPSR[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[1], lt_max[1], m_min[1], m_max[1]));
    n_PLSB[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[0], lt_max[0], m_min[0], m_max[0]));
    n_PRSB[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[0], lt_max[0], m_min[2], m_max[2]));
    n_NLSB[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[1], lt_max[1], m_min[0], m_max[0]));
    n_NRSB[i_pt] = tree->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[i_pt], pt_max[i_pt], lt_min[1], lt_max[1], m_min[2], m_max[2]));
  }
  finD->Close();

  TFile *fin2 = new TFile("../../Store_data_codes/MCm17_cos.root");
  TTree *treeM2 = (TTree*)fin2->Get("MC_cos");

  n_MC[0] = treeM2->GetEntries(Form("dimPt > %f && dimPt < %f && lt > %f && lt < %f && Mass > %f && Mass < %f", pt_min[0], pt_max[0], lt_min[0], lt_max[0], m_min[1], m_max[1]));

  fin2->Close();

  // output results into a .tex table
  ofstream ftex;
  ftex.open(Form("text_output/data_mc_evts_2017.tex"));
  ftex << "\\begin{tabular}{cc|c}\n";
  ftex << "\\hline\n";
  ftex << Form("\\multicolumn{2}{c}{2017} & $[%.0f, %.0f]$ GeV \\\\\n", pt_min[0], pt_max[0]);
  ftex << "\\hline\n";

  // data output
  ftex << Form("\\multirow{6}{*}{\\rotatebox[origin=c]{90}{Data}} & Peak & %.3f M \\\\\n", n_PRSR[0]/1e6);
  ftex << Form("& NP & %.3f M \\\\\n", n_NPSR[0]/1e6);
  ftex << Form("& PR LSB & %.1f k \\\\\n", n_PLSB[0]/1e3);
  ftex << Form("& PR RSB & %.1f k\\\\\n", n_PRSB[0]/1e3);
  ftex << Form("& NP LSB & %.1f k\\\\\n", n_NLSB[0]/1e3);
  ftex << Form("& NP RSB & %.1f k \\\\\n", n_NRSB[0]/1e3);
  ftex << "\\hline\n";

  // MC output
  ftex << Form("MC & (only Peak) & %.3f M  \\\\\n", n_MC[0]/1e6);
  ftex << "\\hline\n";
  
  ftex << "\\end{tabular}\n";
  ftex.close();

}
