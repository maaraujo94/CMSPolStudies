// single muon efficiency
double f_eff(double pt, double eta)
{
  double beta = 1.698;
  double ptz = 3.723;

  double f = 1.0;

  if(abs(eta) < 0.2)
    f /= 1. + exp(-beta*(pt-ptz));

  return 1./f;
}
