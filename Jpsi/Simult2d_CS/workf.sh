#!/bin/bash

echo "storing the main histos"
cd PR_fit

root -l -b -q histoSave.C
root -l -b -q genDist.C

echo "run NP fit framework"
cd ../NP_fit

root -l -b -q genDist.C

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C

echo "run the PR fit framework - part 2"

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C
