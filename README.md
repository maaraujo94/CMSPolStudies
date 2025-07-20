# CMSPolStudies
Framework for CMS S-wave quarkonium polarization studies: Jpsi and Psi2  

Main framework is in the two folders **Jpsi** and **Psi2**. Folder **Comb** just contains codes for plotting final results of both states simultaneously

## Main directories (for each state):

- **Store_data_codes**: The processed root files with the main variables for data and MC events that pass trigger and single muon cuts must be copied over to this area
  - **codes**: contains the codes to run on LXPLUS to process the full datafiles (ang_\*.C). Also contains codes to combine the 2017 and 2018 samples (comb_\*.C). 

- **2017** / **2018**: runs the main analysis for either of the two years

- **Simult**: runs the main analysis, including background subtraction, to get lambda_theta; also contains systematics calculations

- **Simult_[name]**: several directories that rerun the main analysis under different sampling conditions, for systematics checks

- **Phi_fit**: runs the analysis in the phi dimension, auxiliary to the systematics

## Workflow of the analysis

Each folder running the analysis (main or for systematics) has the same structure (some sections may be unnecessary and thus absent for some systematics):
- *ptbins.C* contains the pT binning used
- *workf.sh* runs all the codes in the proper order. Given how long some can take to run, and issues if something fails, I'd suggest checking the order here but running each by hand
- **bkgFits** contains the codes for the PR and NP mass fits
- **cosMax** contains the codes to determine the maximum |costh| to fit in each pt bin
- **NP_fit** contains the codes for the NP background subtraction and lth determination
- **PR_fit** contains the codes for storing the data and MC pT:|costh| histograms, fitting the lifetime, the PR background subtraction and lth determination
- **SBLtFits** contains the codes for fitting the SB lifetime and interpolating the results to fix the NP continuum dimuon term in the SR lifetime fit

### Systematics

The **Simult** folder additionally has a **Systematics** folder, containing the systematics checks that don't require additional cuts to the data and MC sampling (some only present in **Jpsi** or **Psi2**):
- **deltaR** processes the results from the rho systematics checks and reruns baseline fits with proper |costh| limits for comparison
  - **getDR_FW** has the codes for the determination of the fiducial cuts (DeltaR_pT)
- **fBg_evt_ct** runs the systematics for comparison of continuum dimuon background fraction using event counting in the sidebands
- **Lt_Nbg_[plus/minus]** runs the systematics for variation of the N parameter of the lifetime fits
- **MC_lth** has the closure test of the polarization fit procedure using MC
- **M[PR/NP]_alpha_[plus/minus]** runs the systematics for variation of the alpha parameter of the PR and NP mass fits

Also in **Systematics** is the **mainDiffs** folder to determine the final systematics
- *plotAlts.C* and *plotModel.C* plot the resulting lth variation of the checks
- *calcSys.C* determines the final systematic uncertainties
