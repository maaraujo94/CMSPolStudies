#!/usr/bin/env python

import os, imp

os.system("ls */*.C > locs.txt")
fin = open("locs.txt")
locs = fin.readlines()
fin.close()

for i, line in enumerate(locs):
    locs[i] = line[:-1]

bloc = os.getcwd()

for l in locs:
    fin = open("%s/%s"%(bloc, l))
    base = fin.readlines()
    fin.close()

    fout = open("%s/%s"%(bloc, l), "w")
    for line in base:
        if "S_" in line:
            line = line.replace("S_", "17_")
        if "Run 2" in line:
            line = line.replace("Run 2", "2017")
        fout.write(line)
    fout.close()

os.system("rm locs.txt")
