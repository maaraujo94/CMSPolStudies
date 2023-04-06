# CMSPolStudies
framework for CMS S-wave quarkonium polarization studies: Jpsi and Psi2

## Main directories (for each state):

- Store_data_codes: codes to run on LXPLUS, stores the main variables for data and MC events that pass trigger and single muon cuts

- Simult: runs the main analysis, including background subtraction, to get lambda_theta; also contains the systematics checks that don't require additional cuts to the data and MC sampling; also contains codes for auxiliary plots

- 2017 / 2018: runs the main analysis for either of the two years

- Simult_[name]: several directories that rerun the main analysis under different sampling conditions, for systematics checks
