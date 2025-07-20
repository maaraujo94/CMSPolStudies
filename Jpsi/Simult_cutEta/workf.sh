#!/bin/bash

echo "storing the main histos"
cd PR_fit

root -l -b -q histoSave.C

echo "running the mass fits"
cd ../bkgFits

root -l -b -q bkgSave.C
root -l -b -q newDatamass_2.C
root -l -b -q plotDMPars.C
root -l -b -q fbkgProp.C
root -l -b -q newNPmass.C
root -l -b -q plotDMPars_NP.C
root -l -b -q fbkgProp_NP.C

echo "getting costh limits with fine binning"

cd ../cosMax

root -l -b -q histoSave.C
root -l -b -q getCos.C

echo "run SB lifetime fits"
cd ../SBLtFits

root -l -b -q bkgSave.C
root -l -b -q bkgSave_N.C
root -l -b -q ltBkg2d.C
root -l -b -q ltBkg_N.C
root -l -b -q plotLtPars_N.C
root -l -b -q store_SB.C

echo "run PR fit framework - part 1"
cd ../PR_fit

root -l -b -q bkgSave.C
root -l -b -q ltBkg2d.C
root -l -b -q plotLtPars2d.C
root -l -b -q fnpProp.C

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
cd ../PR_fit

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C
