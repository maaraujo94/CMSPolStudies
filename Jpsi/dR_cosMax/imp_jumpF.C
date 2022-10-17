double jumpF(double b_val)
{
  double b_round = floor(b_val*10.)/10.;
  if(b_val-b_round > 0.025) {
    b_round += 0.05;
    if(b_val - b_round > 0.025)
      b_round += 0.05;
  }

  return b_round;
}
