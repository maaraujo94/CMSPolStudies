#!/usr/bin/env python

import os, imp

locs = ["bkgFits/newMCmass_0.C", "bkgFits/NPMCmass.C",
        "cosMax/histoSave.C",
        "PR_fit/histoSave.C", "PR_fit/bkgCosth.C",
        "NP_fit/bkgCosth.C"]

bloc = os.getcwd()

# code to be altered for cuts
for l in locs:
    fin = open("%s/%s"%(bloc, l))
    base = fin.readlines()
    fin.close()

    ct_t = 0
    
    print "now running "+l+"\n"
    fout = open("%s/%s"%(bloc, l), "w")
    fout.write('#import "../effCode.C"\n\n')
    for line in base:
        ct_w = 0
        if "Double_t" in line and ct_t is 0:
            ct_w = 1
            fout.write(line)
            fout.write("double mPPt, mMPt, mPEta, mMEta;\n")
            fout.write("double effP, effM;\n")
            ct_t = 1
        if '"lt", &mc' in line:
            ct_w = 1
            fout.write(line)
            treeL = line.split("->")
            fout.write('%s->SetBranchAddress("muonPEta", &mPEta);\n'%treeL[0])
            fout.write('%s->SetBranchAddress("muonMEta", &mMEta);\n'%treeL[0])
            fout.write('%s->SetBranchAddress("muonPPt", &mPPt);\n'%treeL[0])
            fout.write('%s->SetBranchAddress("muonMPt", &mMPt);\n'%treeL[0])
        if "GetEntry" in line and "treeD" not in line:
            ct_w = 1
            fout.write(line)
            fout.write("	  effP = f_eff(mPPt, mPEta);\n")
            fout.write("	  effM = f_eff(mMPt, mMEta);\n")
        if "Fill(mc" in line or "Fill(abs(cos(mc" in line or "Fill(cos(mc" in line:
            ct_w = 1
            line = line.replace(");", ", effP*effM);")
            fout.write(line)
        if ct_w is 0:
            fout.write(line)
    fout.close()

locR = ["PR_fit/bkgSave.C",
        "PR_fit/ltBkg2d.C", "PR_fit/plotLtPars2d.C",
        "PR_fit/ltBkg.C", "PR_fit/ltPerPt.C", "PR_fit/ltPerPt_muFix.C", "PR_fit/ltPerPt_bFix.C", "PR_fit/plotLtPars.C",
        "PR_fit/fnpProp.C",
        "NP_fit/bkgSave.C"] 

# code to be removed 
for l in locR:
    os.system("rm %s"%l)

locBkg = ["PR_fit/mBkg.C", "PR_fit/fbkgProp.C",
          "NP_fit/mBkg.C", "NP_fit/fbkgProp.C",
          "PR_fit/fNPcorr.C"]

# code to change input source
for l in locBkg:
    fin = open("%s/%s"%(bloc, l))
    base = fin.readlines()
    fin.close()

    ct_t = 0
    ct_w = 0
    
    print "now running "+l+"\n"
    fout = open("%s/%s"%(bloc, l), "w")
    for line in base:
        if "files/mStore" in line:
            subS = l.split("/")[0]
            line = line.replace("files/mStore", "../../Simult/%s/files/mStore"%subS)
            fout.write(line)
        elif "files/NPFrac.root" in line and "update" not in line:
            subS = l.split("/")[0]
            line = line.replace("files/NPFrac.root", "../../Simult/%s/files/NPFrac.root"%subS)
            fout.write(line)
        elif "files/NPFrac.root" in line and "update" in line:
            line = line.replace("update", "recreate")
            fout.write(line)
        else:
            fout.write(line)
    fout.close()
