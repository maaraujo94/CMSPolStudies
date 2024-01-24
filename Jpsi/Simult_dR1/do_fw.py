#!/usr/bin/env python

import os, imp

locs = ["bkgFits/bkgSave.C",
        "cosMax/histoSave.C",
        "PR_fit/histoSave.C", "PR_fit/bkgSave.C"]

bloc = os.getcwd()

for l in locs:
    fin = open("%s/%s"%(bloc, l))
    base = fin.readlines()
    fin.close()

    ct_t = 0

    print "now running "+l+"\n"
    fout = open("%s/%s"%(bloc, l), "w")
    for line in base:
        if '#import "../etacut.C"' in line:
            fout.write('#import "../rcut.C"\n')
        elif "double mPEta" in line:
            fout.write("Double_t dR;\n")
        elif '"muonPEta"' in line:
            continue;
        elif '"muonMEta"' in line:
            treeL = line.split("->")
            fout.write('%s->SetBranchAddress("DeltaR", &dR);\n'%treeL[0])
        elif "abs(mPEta)" in line:
            fout.write("if(dR > r_cut)\n")
        else:
            fout.write(line)
        
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
        if "double cMax[" in line:
            ct_w = 1
            fout.write(line)
            fout.write("\n")
            fout.write("double cMin[nBinsY];\n")
        if "double cMaxVal" in line:
            ct_w = 1
            fout.write("double cMaxVal = jumpF(cosMax->Eval(pMin));\n")
            fout.write("double cMinVal = jumpF(cosMin->Eval(pMax));\n")
        if "cMax[i] =" in line:
            ct_w = 1
            fout.write(line)
            fout.write("\n")
            fout.write("cMin[i] = cMinVal;\n")
        if "SetParameters(p" in line:
            ct_w = 1
            line = line.replace("GetBinContent(1)*1.1", "GetMaximum()")
            fout.write(line)
        if "SetRange(0, cMaxVal)" in line:
            ct_w = 1
            fout.write("fit1d[i_t]->SetRange(cMinVal, cMaxVal);\n")
        if "SetRange(0, cMax[i])" in line:
            ct_w = 1
            fitL = line.split("->")
            fout.write("%s->SetRange(cMin[i], cMax[i]);\n"%fitL[0])
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
