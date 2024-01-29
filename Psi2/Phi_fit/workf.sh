#!/bin/bash

echo "storing the main histos"
cd PR_fit

root -l -b -q histoSave.C

echo "running the mass fits"
cd ../bkgFits

root -l -b -q fbkgProp.C
root -l -b -q fbkgProp_NP.C

echo "run PR fit framework - part 1"
cd ../PR_fit

root -l -b -q fnpProp.C
root -l -b -q genDist.C

echo "run NP fit framework"
cd ../NP_fit

root -l -b -q genDist.C
root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C

echo "run the PR fit framework - part 2"
cd ../PR_fit

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C
