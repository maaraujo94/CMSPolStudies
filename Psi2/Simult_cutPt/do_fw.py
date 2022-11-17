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
            fout.write('#import "../ptcut.C"\n')
        elif "double mPEta" in line:
            fout.write("double mPPt, mMPt, mPEta, mMEta;\n")
        elif '"muonMEta"' in line:
            fout.write(line)
            treeL = line.split("->")
            fout.write('%s->SetBranchAddress("muonPPt", &mPPt);\n'%treeL[0])
            fout.write('%s->SetBranchAddress("muonMPt", &mMPt);\n'%treeL[0])
        elif "abs(mPEta)" in line:
            fout.write("if((abs(mPEta) > eta_lim || mPPt > pt_cut) && (abs(mMEta) > eta_lim || mMPt > pt_cut))\n")
        else:
            fout.write(line)
        
    fout.close()
    
