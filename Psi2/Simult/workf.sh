#!/bin/bash

echo "running the mass fits"
cd bkgFits

root -l -b -q newMCmass_0.C
root -l -b -q newMCmass_1.C
root -l -b -q newMCmass_2.C
root -l -b -q newMCmass_3.C
root -l -b -q newMCmass_4.C
root -l -b -q plotMMseq.C
root -l -b -q plotnalpha.C
root -l -b -q newMCmass_G.C
root -l -b -q plotMMG.C

root -l -b -q bkgSave.C
root -l -b -q newDatamass_0.C
root -l -b -q newDatamass_1.C
root -l -b -q plotDMPars.C
root -l -b -q fbkgProp.C
root -l -b -q newNPmass.C
root -l -b -q plotDMPars_NP.C
root -l -b -q fbkgProp_NP.C

echo "getting costh limits with fine binning"

cd ../cosMax

root -l -b -q histoSave.C
root -l -b -q getCos.C

echo "run PR fit framework - part 1"
cd ../PR_fit

root -l -b -q histoSave.C

root -l -b -q bkgSave.C 
#root -l -b -q ltBkg2d.C
#root -l -b -q plotLPars2d.C
#root -l -b -q ltBkg.C
root -l -b -q fnpProp.C
root -l -b -q fNPcorr.C

root -l -b -q getfL.C
root -l -b -q genDist.C

echo "run NP fit framework"
cd ../NP_fit

root -l -b -q getfL.C 
root -l -b -q genDist.C

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C

echo "run the PR fit framework - part 2"

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C
