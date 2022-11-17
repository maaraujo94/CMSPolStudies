// weight for lambda_phi correction
double f_phi(double phi)
{
  double beta = -0.005;
  
  double f = 1 + beta*cos(2*phi*TMath::Pi()/180.);
    
  return f;
}
