fname = "text_output/massPlot_Psi2.tex"
fout = open(fname, "w")

nJobs = 8
ptBins = [20, 25, 30, 35, 40, 50, 60, 70, 100];

w_v = 0.4

ct = 0
for i in range(nJobs):
    fout.write("\\begin{figure}[h!]\n")
    fout.write("\\centering\n")
    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Simult/bkgFits/plots/mass/fit_1/fit_pt%d.pdf}\n"%(w_v, i))
    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Simult/bkgFits/plots/massNP/fit/fit_pt%d.pdf}\n\n"%(w_v, i))
    
    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Simult/bkgFits/plots/mass/fit_1/pulls_pt%d.pdf}\n"%(w_v, i))
    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Simult/bkgFits/plots/massNP/fit/pulls_pt%d.pdf}\n\n"%(w_v, i))

    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Simult/bkgFits/plots/mass/fit_1/devs_pt%d.pdf}\n"%(w_v, i))
    fout.write("\\includegraphics[width=%f\\textwidth]{/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Simult/bkgFits/plots/massNP/fit/devs_pt%d.pdf}\n"%(w_v, i))
    fout.write("\\caption{Fitted $\\psip$ mass distributions for $\\pt$ bin %.1f-%.1f GeV (top), corresponding pulls distributions (middle) and relative differences (bottom). In the prompt (left) and non-prompt lifetime intervals (right).}\\label{f:m_fit_%d}\n"%(ptBins[i], ptBins[i+1], i))
    fout.write("\\end{figure}\n")
    fout.write("\n\\pagebreak\n\n")
    
fout.close()
