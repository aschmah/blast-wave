#include <TMath.h>
#include <string>

void calFitRange(double beta = 0.68)
{
  const int NumOfParticles = 9;
  double mass[NumOfParticles] = {0.13957,0.493677,0.938272,1.019460,1.67245,1.86962,3.096916,9.46030,1.875612};
  string pid[NumOfParticles]  = {"pi","K","p","phi","Omega","D0","J/psi","Y","d"};

  double gamma = 1.0/TMath::Sqrt(1.0-beta*beta);
  for(int i_pid = 0; i_pid < NumOfParticles; ++i_pid)
  {
    double pT = mass[i_pid]*gamma*beta;
    cout << "radial momentum for " << pid[i_pid].c_str() << " is " << pT+1.0 << endl;
  }
}
