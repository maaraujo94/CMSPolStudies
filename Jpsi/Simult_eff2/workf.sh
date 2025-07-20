#!/bin/bash

echo "getting costh limits with fine binning\n"

cd cosMax
root -l -b -q histoSave.C
root -l -b -q getCos.C

echo "saving the PR histos for fitting"

cd ../PR_fit
root -l -b -q histoSave.C

echo "running NP framework\n"

cd ../NP_fit

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C

echo "running bkg sub with corrected NP\n"

cd ../PR_fit

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C
