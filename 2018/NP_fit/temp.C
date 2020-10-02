void temp()
{
  int nBins = 45;
  int nPlot = 6;
  int nRow = 2;

  for(int i = 0; i < nBins; i++) {
    if(i%6 == 0){
      cout << "\\begin{figure}[h!]" << endl;
      cout << "\\centering" << endl;
    }
    
    cout << Form("\\includegraphics[width = 0.35\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/2018/NP_fit/plots/ratio_final/bin_%d.pdf}", i+1) << endl;
    if(i%2 == 1)
      cout << endl;

    if(i%6 == 5) {
      cout << "\\end{figure}" << endl;
      cout << endl << "\\pagebreak" << endl << endl;
    }
  }

  if(nBins%6 != 0)
    cout << "\\end{figure}" << endl;
}
