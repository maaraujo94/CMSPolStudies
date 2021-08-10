#!/bin/bash

echo "saving the histos for fitting"
echo

# root -l -q histoSave.C

# echo "getting costh limits with fine binning\n"

# cd ../cosMax
# root -l -q histoSave.C
# root -l -q getCos.C

# echo "fitting lifetime, mass distributions\n"

# cd ../PR_fit
# root -l -q bkgSave.C
# root -l -q ltBkg.C
# root -l -q mBkg.C

# echo "fitting mass bkg costh dists\n"

# root -l -q bkgCosth.C
# root -l -q fitBkgCosth.C
# root -l -q fitBkgCosth2d.C
# root -l -q plotCosPars_both.C

# echo "getting mass bkg fitted dists in reg binning\n"

# root -l -q getfL.C
# root -l -q genDist.C

# echo "costh bkg subtraction\n"

# root -l -q fitFrac.C
# root -l -q bkgSub.C

# echo "final fitting and plotting\n"

# root -l -q indFit.C
# root -l -q plotRes.C
