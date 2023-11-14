fname = "text_output/massPlot_Jpsi.tex"
fout = open(fname, "w")

nJobs = 19
ptBins = [25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 55, 60, 65, 70, 75, 80, 90, 100, 120];

w_v = 0.4

ct = 0
for i in range(nJobs):
    fout.write("\\begin{figure}[h!]\n")
    fout.write("\\centering\n")
    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/bkgFits/plots/mass/fit_2/fit_pt%d.pdf}\n"%(w_v, i))
    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/bkgFits/plots/massNP/fit/fit_pt%d.pdf}\n\n"%(w_v, i))
    
    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/bkgFits/plots/mass/fit_2/pulls_pt%d.pdf}\n"%(w_v, i))
    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/bkgFits/plots/massNP/fit/pulls_pt%d.pdf}\n\n"%(w_v, i))

    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/bkgFits/plots/mass/fit_2/devs_pt%d.pdf}\n"%(w_v, i))
    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/bkgFits/plots/massNP/fit/devs_pt%d.pdf}\n"%(w_v, i))
    fout.write("\\caption{Fitted $\\jpsi$ mass distributions for $\\pt$ bin %.1f-%.1f GeV (top), corresponding pulls distributions (middle) and relative differences (bottom). In the prompt (left) and non-prompt lifetime intervals (right).}\\label{f:m_fit_%d}\n"%(ptBins[i], ptBins[i+1], i))
    fout.write("\\end{figure}\n")
    fout.write("\n\\pagebreak\n\n")
    
fout.close()
