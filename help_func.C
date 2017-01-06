// These functions are needed for the 1st round of reconstruction to
// get an estimate on vertex posision. 

//#include <map>
#include <iterator>
#include <algorithm>

//map<int, double> INDEX;

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

// Divides sphere in NPHIxNTHETA segments and takes 1st photon from each segment
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


int MarkEarlyPhotons_(int N, float* x, float* y, float* z, float* t, int* process, int* pe, int* ph_vec, int num=1)
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

//------------ make a loop to select first num photons
      for(int n=1;n<=num;++n)
      {
        if(map[ip][it].size()==0) break; // this may happen if num > initial size of map[ip][it]
        int I=0;
        int n1st=map[ip][it][I];
        for(int i=1;i<map[ip][it].size();++i)
        {
          if(t[map[ip][it][i]]<t[n1st])
          {
            n1st=map[ip][it][i];
	    I=i;
          }
        }
        ph_vec[n1st]=n;
        map[ip][it].erase(map[ip][it].begin()+I);// current earlies is marked. remove it from futher consideartion

      }
//================

    }


  return 0;
}



int FillIndex(char* fName)
{
  int wl;
  double n_ref;
  ifstream infile(fName);
  while(infile>>wl>>n_ref)
  {
    cout<<"nref = "<<n_ref<<endl;
    INDEX[wl] = n_ref;
  }
  
  //now caculate ref. index for group velocity
  double dn=0;
  double dl=0;
  double n=0;
  double l=0;
  map<int, double>::iterator last = INDEX.end();
  last--;
  for(map<int, double>::iterator m=INDEX.begin(); m!=INDEX.end(); ++m)
  {
    cout<<m->first<<"   "<<m->second<<endl;
    n = m->second;
    l = m->first;
    map<int, double>::iterator next = m;
    map<int, double>::iterator prev = m;
    next++;
    prev--;
    if(m==INDEX.begin())
    {
      dn = next->second - m->second;
      dl = next->first - m->first;
    } else 
    if(m==last)
    {
      dn = m->second - prev->second;
      dl = m->first - prev->first;
    } else
    {
      dn = next->second - prev->second;
      dl = next->first - prev->first;
    }
    INDEX_GR[l] = n - l*(dn/dl); // group velocity will be C_VAC/INDEX_GR[]
  }
  //some output for testing:
  cout<<"INDEX_GR[400] = "<<INDEX_GR[400]<<endl;
  for(int i=400;i!=450;++i)
  {
    cout<<i<<"   "<<INDEX_GR[i]<<"   "<<INDEX[i]<<endl;
  }
  return 0;
}

double Velocity(double lambda)
{
  double v=0.;
  return v;
}

//The followinf is copied from M.Wetstein's WChSanBox package
double _c=C_VAC; //this line added for compatibility

Double_t CalcMaxAlpha(Double_t dD, Double_t dT, Double_t ni)
{
   cout<<"dD = "<<dD<<"   dT = "<<dT<<"   ni = "<<ni<<endl;
   cout<<"dD/dT = "<<dD/dT<<"   c/n = "<<_c/ni<<"   dv = "<<dD/dT-_c/ni<<endl;
   double qA=4*(ni*ni*ni*ni)*dD*dD;
   double qB=-8*(ni*ni)*(_c*dT)*dD;
   double qC= (4*dT*dT*_c*_c) - (4*((ni*ni) - 1)*( (ni*ni*dD*dD) - (dT*dT*_c*_c) ));
//   cout<<"qA = "<<qA<<"   qBi"

   double maxa=-555555;
   if( (qB*qB) > (4*qA*qC) ){
     double cosa = (-qB + sqrt(qB*qB - 4*qA*qC))/(2*qA);
     maxa = acos(cosa);

     double cosalphap = cos(maxa + 0.01);
     double qA1=((ni*ni)-1);
     double qB1=2*((_c*dT)-(ni*ni*dD*cosalphap));
     double qC1=(ni*ni*dD*dD)-(_c*_c*dT*dT);

     double cosalpham = cos(maxa - 0.01);
     qA1=((ni*ni)-1);
     qB1=2*((_c*dT)-(ni*ni*dD*cosalpham));
     qC1=(ni*ni*dD*dD)-(_c*_c*dT*dT);

   } else{
   }
   return maxa;
}

Double_t CalcS1(Double_t alpha, Double_t dD, Double_t dT, Double_t ni)
{
   double theS1=-55555.;

   if(alpha>=0){
     double cosalpha = cos(alpha);

     double qA=((ni*ni)-1);
     double qB=2*((_c*dT)-(ni*ni*dD*cosalpha));
     double qC=(ni*ni*dD*dD)-(_c*_c*dT*dT);
     double sqrtterm=qB*qB - 4*qA*qC;
     if( fabs(sqrtterm)<0.0001 ) sqrtterm=0;
     if( sqrtterm>=0 ){

       theS1 = (-qB + sqrt(sqrtterm))/(2*qA);
     }
     else{
     }
   }
   return theS1;
}

Double_t CalcS2(Double_t alpha, Double_t dD, Double_t dT, Double_t ni)
{
  double mins1 = CalcS1(alpha,dD,dT,ni);
  double ds2x = dD - mins1*cos(alpha);
  double ds2y = mins1*sin(alpha);
  double theS2 = sqrt( ds2x*ds2x + ds2y*ds2y );

  return theS2;
}

