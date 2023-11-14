#!/bin/bash

echo "getting costh limits with fine binning\n"

cd cosMax
root -l -b -q histoSave.C
root -l -b -q getCos.C

echo "storing base histos"

cd ../PR_fit
root -l -b -q histoSave.C

echo "running NP framework\n"

cd ../NP_fit

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C

echo "running PR framework\n"

cd ../PR_fit

root -l -b -q bkgSub.C
root -l -b -q indFit.C
root -l -b -q plotRes.C
