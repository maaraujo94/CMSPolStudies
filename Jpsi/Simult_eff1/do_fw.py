#!/usr/bin/env python

import os, imp

locs = ["cosMax/histoSave.C",
        "PR_fit/histoSave.C"]

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

locBkg = ["NP_fit",
          "PR_fit"]

# code to change input source
for l in locBkg:
    fin = open("%s/%s/bkgSub.C"%(bloc, l))
    base = fin.readlines()
    fin.close()
    
    print "now running "+l+"\n"
    fout = open("%s/%s/bkgSub.C"%(bloc, l), "w")
    for line in base:
        if "TFile *inBkg" in line:
            fout.write('TFile *inBkg = new TFile("../../Simult/%s/files/bkgCosModel.root");\n'%l)
        elif "TFile *inFracSB" in line:
            if l is "NP_fit":
                fout.write('TFile *inFracSB = new TFile("../../Simult/bkgFits/files/bkgFrac_NP.root");\n')
            else:
                fout.write('TFile *inFracSB = new TFile("../../Simult/bkgFits/files/bkgFrac.root");\n')                
        elif "TFile *inFracNP" in line:
            fout.write('TFile *inFracNP = new TFile("../../Simult/PR_fit/files/NPFrac.root");')
        else:
            fout.write(line)
    fout.close()
