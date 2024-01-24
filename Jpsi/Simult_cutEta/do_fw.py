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
    fout.write('#import "../etacut.C"\n\n')
    for line in base:
        fout.write(line)
        if "Double_t" in line and ct_t is 0:
            fout.write("double mPEta, mMEta;\n")
            ct_t = 1
        if '"lt"' in line:
            treeL = line.split("->")
            fout.write('%s->SetBranchAddress("muonPEta", &mPEta);\n'%treeL[0])
            fout.write('%s->SetBranchAddress("muonMEta", &mMEta);\n'%treeL[0])
        if "GetEntry" in line:
            fout.write("if((abs(mPEta) < eta_lo || abs(mPEta) > eta_hi) && (abs(mMEta) < eta_lo || abs(mMEta) > eta_hi))\n")
    fout.close()
