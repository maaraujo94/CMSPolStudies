#!/bin/bash

cd PR_fit

root -l -b -q histoSave.C

cd ../bkgFits

root -l -b -q fbkgProp.C
root -l -b -q fbkgProp.C

cd ../PR_fit

root -l -b -q fnpProp.C
root -l -b -q fNPcorr.C
root -l -b -q genDist.C

cd ../NP_fit

root -l -b -q genDist.C
root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C

cd ../PR_fit

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C
