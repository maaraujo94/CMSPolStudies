#!/bin/bash

echo "running MC fits\n"

cd bkgFits

#root -l -b -q newMCmass_0.C
#root -l -b -q newMCmass_1.C
#root -l -b -q newMCmass_2.C
#root -l -b -q newMCmass_3.C
#root -l -b -q newMCmass_4.C
#root -l -b -q newMCmass_5.C
#root -l -b -q newMCmass_G.C
#root -l -b -q plotMMseq.C
#root -l -b -q plotnalpha.C
#root -l -b -q plotMMG.C

echo "running NP MC fits\n"

#root -l -b -q NPMCmass.C

echo "getting costh limits with fine binning\n"

#cd ../cosMax
#root -l -b -q histoSave.C
#root -l -b -q getCos.C

echo "saving the PR histos for fitting"

cd ../PR_fit
root -l -b -q histoSave.C

echo "fitting lifetime, mass distributions\n"

root -l -b -q mBkg.C

echo "getting mass bkg dists in reg binning\n"

root -l -b -q bkgCosth.C
root -l -b -q getfL.C
root -l -b -q genDist.C

echo "getting bkg frac uncertainties\n"

root -l -b -q fbkgProp.C

echo "running NP framework\n"

cd ../NP_fit

root -l -b -q mBkg.C
root -l -b -q bkgCosth.C
root -l -b -q getfL.C
root -l -b -q genDist.C
root -l -b -q fbkgProp.C
root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C

echo "running bkg sub with corrected NP\n"

cd ../PR_fit

root -l -b -q fNPcorr.C
root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C
