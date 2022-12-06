#!/bin/bash

cd cosMax
root -l -b -q histoSave.C
root -l -b -q getCos.C

cd ../PR_fit
root -l -b -q histoSave.C

cd ../NP_fit

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C

cd ../PR_fit

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C
