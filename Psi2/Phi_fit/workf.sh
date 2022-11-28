#!/bin/bash

echo "saving the PR histos for fitting"

cd ../PR_fit
root -l -b -q histoSave.C

echo "fitting mass bkg costh dists\n"

root -l -b -q bkgCosth.C

echo "getting mass bkg fitted dists in reg binning\n"

root -l -b -q genDist.C

echo "getting bkg frac uncertainties\n"

root -l -b -q fbkgProp.C
root -l -b -q fnpProp.C

echo "running NP framework\n"

cd ../NP_fit

root -l -b -q bkgCosth.C
root -l -b -q genDist.C
root -l -b -q fbkgProp.C
root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C

echo "running bkg sub with corrected NP\n"

cd ../PR_fit

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C
