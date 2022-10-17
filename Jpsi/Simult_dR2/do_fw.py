#!/usr/bin/env python

import os, imp

locs = ["bkgFits/newMCmass_0.C", "bkgFits/NPMCmass.C",
        "cosMax/histoSave.C",
        "PR_fit/histoSave.C", "PR_fit/bkgSave.C", "PR_fit/bkgCosth.C",
        "NP_fit/bkgSave.C", "NP_fit/bkgCosth.C"]

bloc = os.getcwd()

# updating code for cuts on samples
for l in locs:
    fin = open("%s/%s"%(bloc, l))
    base = fin.readlines()
    fin.close()

    ct_t = 0

    print "now running "+l+"\n"
    fout = open("%s/%s"%(bloc, l), "w")
    fout.write('#import "../rcut.C"\n\n')
    for line in base:
        fout.write(line)
        if "Double_t" in line and ct_t is 0:
            fout.write("Double_t dR;\n")
            if l is not "PR_fit/bkgSave.C":
                ct_t = 1
        if '"lt"' in line:
            treeL = line.split("->")
            fout.write('%s->SetBranchAddress("DeltaR", &dR);\n'%treeL[0])
        if "GetEntry" in line:
            fout.write("if(dR > r_cut)\n")
    fout.close()

# replacing code for cosMax procedure 
os.system("cp ../dR_cosMax/getCos.C cosMax")
os.system("cp ../dR_cosMax/getCosMin.C cosMax")
os.system("cp ../dR_cosMax/imp_jumpF.C cosMax")
os.system("cp ../dR_cosMax/plotCos.C cosMax")

# updating code for final fit
locF = ["PR_fit", "NP_fit"]
for f in locF:
    fin = open("%s/%s/indFit.C"%(bloc, f))
    base = fin.readlines()
    fin.close()

    ct_l = 0

    fout = open("%s/%s/indFit.C"%(bloc, f), "w")
    fout.write('#import "../rcut.C"\n\n')
    for line in base:
        ct_w = 0
        if "code to" in line:
            ct_w = 1
            fout.write(line)
            fout.write("\n")
            fout.write("// aux func for costheta_min\n")
            fout.write("double cminf(double pt, double a, double b, double c, double d)\n")
            fout.write("{\n")
            fout.write("  if (pt < d) return 0;\n")
            fout.write("  else\n")
            fout.write("    return a*(1.-exp(b+c*pt));\n")
            fout.write("}\n")
        if "maxPar[0]," in line:
            ct_w = 1
            fout.write(line)
            fout.write("\n")
            fout.write('  in.open("../cosMax/cosMinFitRes.txt");\n')
            fout.write("  getline(in, dataS);\n")
            fout.write("  getline(in, dataS);\n")
            fout.write("  double minPar[4];\n")
            fout.write("  in >> minPar[0] >> aux >> minPar[1] >> aux >> minPar[2] >> aux >> minPar[3];\n")
            fout.write("  in.close();\n")
            fout.write("\n")
            fout.write('  TF1 *cosMin = new TF1("cosMin", "cminf(x, [0], [1], [2], [3])", yBins[0]-10, yBins[nBinsY]+10);\n')
            fout.write("  cosMin->SetParameters(minPar[0], minPar[1], minPar[2], minPar[3]);\n")
        if "double cMaxVal" in line:
            ct_w = 1
            fout.write("double cMaxVal = jumpF(cosMax->Eval(pMin));\n")
            fout.write("double cMinVal = jumpF(cosMin->Eval(pMax));\n")
        if "SetRange" in line:
            ct_w = 1
            fout.write("fit1d[i_t]->SetRange(cMinVal, cMaxVal);\n")
        if "SetMaximum" in line:
            ct_w = 1
            fout.write("    pHist[0][i]->SetMaximum(parA[0][i]*1.5);\n")
            fout.write("    if(i == nBinsY-1) pHist[0][i]->SetMaximum(pHist[0][i]->GetMaximum()*1.5);\n")
        if "c_lim" in line:
            ct_w = 1
            line = line.replace("c_lim", "c_lim_Max")
            fout.write(line)
            ct_l += 1
        if ct_l is 4:
            ct_w = 1
            fout.write("    TLine *c_lim_Min = new TLine(cMinVal, 0, cMinVal, pHist[0][i]->GetMaximum());\n")
            fout.write("    c_lim_Min->SetLineStyle(kDashed);\n")
            fout.write("    c_lim_Min->SetLineColor(kBlack);\n")
            fout.write("    c_lim_Min->Draw();\n")
            ct_l += 1
        if ct_w is 0:
            fout.write(line)
