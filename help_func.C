int TowerIPhi(double phi)
{
  double two_pi = 2.0*3.1415926535;
  if(phi<0.0) phi = two_pi+phi;
  return int(phi*NPHI/two_pi);
}

int TowerITheta(double theta)
{
  return int((1+TMath::Cos(theta))*NTHETA/2.);
}

int MarkEarlyPhotons(int N, float* x, float* y, float* z, float* t, int* process, int* pe, bool* ph_vec)
{
//Fill in iphi-itheta towers with photon number (only photons passing QE&P cuts)
////QE&P = quantum efficiency and creation process
  std::vector<int> map[NPHI][NTHETA];
  for(int i=0;i!=N;++i)
  {
//    if(process) continue;
    if(pe[i]==0) continue;

    TVector3 vec(x[i],y[i],z[i]);
    int iphi=TowerIPhi(vec.Phi());
    int itheta=TowerITheta(vec.Theta());
    map[iphi][itheta].push_back(i);
  }
// at this point map is created from photons passing QE&P cuts
// let's mark earliest in each sphere segment
  for(int ip=0;ip!=NPHI;++ip)
    for(int it=0;it!=NTHETA;++it)
    {
      if(map[ip][it].size()==0) continue;
      int n1st=map[ip][it][0];
      for(int i=1;i<map[ip][it].size();++i)
      {
        if(t[map[ip][it][i]]<t[n1st])
          n1st=map[ip][it][i];
      }
      ph_vec[n1st]=1;
    }


  return 0;
}
