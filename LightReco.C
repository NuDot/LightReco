/* Andrey Elagin, September 5, 2014
 * LightReco is a light-weight standalone vertex and (in the near future) directionality
 * reconsruction code for 0vbb-decay events in liquid scintillator. The code is based on
 * quadruplet-based vertex-finding method by Michael Smy. Many lines are directly copied
 * from WCSimAnalysis package.
 *
 * See README.txt for instructions on how to use and notes on significant updates.
 */
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TError.h"
#include "TMinuit.h"
#include "TRandom.h"

#include "globals.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <iterator>

int EVT_NUM=100; //controls maximum number of events to be processed
double R_SPHERE=650; //sphere diameter [cm]
double N_REF=1.53; //average index of refraction
double C_VAC=29.9792458; //speed of light in vacuum [cm/ns]
int NSeedsTarget=400; //number of quadruplets
double TSIGMA=0.5; //total time spread (including detector TTS chromatic dispersions)
static const int NMAX_PHOT=100000; //
static const int NPHI=8;
static const int NTHETA=8;
double RECO_DT=3;
double MOM_DT=1.0;
double MAX_FIT_DIGITS=50;
int RECO_MODE=0;
double VTX_SMEAR=0.0;
double VTX_SHIFT_X=5.0;
double VTX_SHIFT_Y=0;
double VTX_SHIFT_Z=0;

using namespace std;

map<int, double> INDEX;
map<int, double> INDEX_GR;

#include "help_func.C"

double fBaseFOM=100.0; //Figure of merit. Borrowed from WCSim: the higher it is the better
double meanTime=0.;
double seedTime=0.;

TRandom RND;
TRandom rndVtx;

// store seed vertex calculated from quaruplets
vector<double> vSeedVtxX;
vector<double> vSeedVtxY;
vector<double> vSeedVtxZ;
vector<double> vSeedVtxTime;
vector<int> vSeedDigitList;

// store photon hits after filtering cuts (e.g. position dependent cut to increase cherenkov fraction)
vector<double> fDigitX;
vector<double> fDigitY;
vector<double> fDigitZ;
vector<double> fDigitT;
vector<double> fDigitQ;
vector<double> fDigitPE;
vector<double> fDigitW;
vector<double> fDigitV;
vector<double> fDigitVgr;
vector<double> fDigitN;
vector<double> fDigitNgr;
vector<double> fDelta; // time residual

int fNDigits=0;
int fThisDigit=0;
int fLastEntry=0;
int fCounter=0;
int fMinTime=0;



//this is for diagnostics, not finished yet
TFile fFOM("fFOM.root","recreate");
TH1F* hT = new TH1F("hT","hT",100,-10,10);
TH1F* hDT0 = new TH1F("hDT0","hDT0",100,-5,5);
TH1F* hDT = new TH1F("hDT","hDT",100,-5,5);

// This function solves system of 4 equations with 4 unknowns to find the seed vertex for each quadruple
int FindVertex(Double_t x0, Double_t y0, Double_t z0, Double_t t0, Double_t x1, Double_t y1, Double_t z1, Double_t t1, Double_t x2, Double_t y2, Double_t z2, Double_t t2, Double_t x3, Double_t y3, Double_t z3, Double_t t3, Double_t& vxm, Double_t& vym, Double_t& vzm, Double_t& vtm, Double_t& vxp, Double_t& vyp, Double_t& vzp, Double_t& vtp)
{
  vxm = -99999.9;
  vym = -99999.9;
  vzm = -99999.9;
  vtm = -99999.9;

  vxp = -99999.9;
  vyp = -99999.9;
  vzp = -99999.9;
  vtp = -99999.9;

  // speed of light in water
  // ======================= 
  Double_t c = C_VAC/N_REF;

  // causality checks
  // ================
  if( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) >= c*c*(t1-t0)*(t1-t0)
   && (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) >= c*c*(t2-t1)*(t2-t1)
   && (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2) >= c*c*(t3-t2)*(t3-t2)
   && (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0) >= c*c*(t2-t0)*(t2-t0)
   && (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1) >= c*c*(t3-t1)*(t3-t1)
   && (x3-x0)*(x3-x0) + (y3-y0)*(y3-y0) + (z3-z0)*(z3-z0) >= c*c*(t3-t0)*(t3-t0) ){

    // [Note: for causality, require that |x_{i}-x_{j}| >= c*|t_{i}-t_{j}|
    //        for each pair of points]
    Double_t dx1 = x1-x0;  Double_t dy1 = y1-y0;  Double_t dz1 = z1-z0;  Double_t dt1 = c*(t1-t0);
    Double_t dx2 = x2-x0;  Double_t dy2 = y2-y0;  Double_t dz2 = z2-z0;  Double_t dt2 = c*(t2-t0);
    Double_t dx3 = x3-x0;  Double_t dy3 = y3-y0;  Double_t dz3 = z3-z0;  Double_t dt3 = c*(t3-t0);

    Double_t epsilon = 1.0e-7;

    // check that points don't all lie in a plane
    if( !( fabs(dx1)<epsilon && fabs(dx2)<epsilon && fabs(dx3)<epsilon )
     && !( fabs(dy1)<epsilon && fabs(dy2)<epsilon && fabs(dy3)<epsilon )
     && !( fabs(dz1)<epsilon && fabs(dz2)<epsilon && fabs(dz3)<epsilon )
     && !( fabs(dx1)<epsilon && fabs(dy1)<epsilon && fabs(dz1)<epsilon )
     && !( fabs(dx2)<epsilon && fabs(dy2)<epsilon && fabs(dz2)<epsilon )
     && !( fabs(dx3)<epsilon && fabs(dy3)<epsilon && fabs(dz3)<epsilon ) ){

    // [Note: this is a problem for detectors with flat faces!]

      Double_t Mdata[9] = { dx1, dy1, dz1,
                            dx2, dy2, dz2,
                            dx3, dy3, dz3 };

      Double_t Qdata[3] = { 0.5*( dx1*dx1 + dy1*dy1 + dz1*dz1 - dt1*dt1 ),
                            0.5*( dx2*dx2 + dy2*dy2 + dz2*dz2 - dt2*dt2 ),
                            0.5*( dx3*dx3 + dy3*dy3 + dz3*dz3 - dt3*dt3 ) };

      Double_t Tdata[3] = { dt1,
                            dt2,
                            dt3 };

      TMatrixD M(3,3,Mdata);
      TMatrixD Q(3,1,Qdata);
      TMatrixD T(3,1,Tdata);

      if( M.Determinant() != 0.0 ){

        TMatrixD A(3,1);
        TMatrixD B(3,1);

        M.Invert();
        A.Mult(M,T);
        B.Mult(M,Q);

        Double_t ax = A(0,0);
        Double_t ay = A(1,0);
        Double_t az = A(2,0);

        Double_t bx = B(0,0);
        Double_t by = B(1,0);
        Double_t bz = B(2,0);

        Double_t ab = ax*bx + ay*by + az*bz;
        Double_t a2 = ax*ax + ay*ay + az*az;
        Double_t b2 = bx*bx + by*by + bz*bz;

        Double_t qa = a2-1.0;
        Double_t qb = 2.0*ab;
        Double_t qc = b2;

        // check for solutions
        if( qb*qb-4.0*qa*qc>0.0 ){

          // The common vertex is given by a quadratic equation, which has two solutions.
          // Typically, one solution corresponds to photons travelling forwards in time,
          // and the other solution corresponds to photons travelling backwards in time.
          // However, sometimes there appear to be two valid solutions.

          Double_t ctm = ( -qb - sqrt(qb*qb-4.0*qa*qc) ) / ( 2.0*qa );
          Double_t ctp = ( -qb + sqrt(qb*qb-4.0*qa*qc) ) / ( 2.0*qa );

          Double_t tm = t0 + ctm/c;
          Double_t xm = x0 + ctm*ax + bx;
          Double_t ym = y0 + ctm*ay + by;
          Double_t zm = z0 + ctm*az + bz;
          Bool_t foundVertexM = 0;

          if( tm<t0 && tm<t1
           && tm<t2 && tm<t3 ){
            vxm = xm;
            vym = ym;
            vzm = zm;
            vtm = tm;
            foundVertexM = 1;
          }

          Double_t tp = t0 + ctp/c;
          Double_t xp = x0 + ctp*ax + bx;
          Double_t yp = y0 + ctp*ay + by;
          Double_t zp = z0 + ctp*az + bz;
          Bool_t foundVertexP = 0;

          if( tp<t0 && tp<t1
           && tp<t2 && tp<t3 ){
            vxp = xp;
            vyp = yp;
            vzp = zp;
            vtp = tp;
            foundVertexP = 1;
          }

        }
        else
         {
           std::cout << "qb*qb-4.0*qa*qc<0.0." << std::endl;
         }
      }
      else
       {
         std::cout << "M.Determinant() == 0.0. " << std::endl;
       }

    }
    else
     {
       std::cout << "The 4 points lie on a plane. " << std::endl;
     }
  }
  else
  {
    std::cout << "Causality check for this 4-group of digits has not been passed. " << std::endl;
  }


  return 0;
}

// Randomly select photon hit (x,y,z,t) for a quadruple
int ChooseNextDigit(Double_t& xpos, Double_t& ypos, Double_t& zpos, Double_t& time)
{
  xpos=0; ypos=0; zpos=0; time=0;
  // ROOT random number generator
  Double_t r = RND.Rndm();

  cout<<"I'm inside ChooseNextDitig"<<endl;
  // pseudo-random number generator
  Int_t numEntries = vSeedDigitList.size();
/*  cout<<"fCounter = "<<fCounter<<"   fLastEntry = "<<fLastEntry<<endl;
  fCounter++;
  if( fCounter>=fNDigits ) fCounter = 0;
  fThisDigit = vSeedDigitList.at(fLastEntry);

  Double_t t0 = 0.5 + fDigitT[fCounter] - fMinTime;
  Double_t q0 = 0.5 + fDigitQ[fCounter];

  Double_t t1 = 0.5 + fDigitT[fThisDigit] - fMinTime;
  Double_t q1 = 0.5 + fDigitQ[fThisDigit];

  Double_t tq = 100.0*(t0*q0+t1*q1);
  Double_t r = tq - TMath::Floor(tq);
*/
//  r = gRandom->Uniform(); // Christoph Aberle, August 14, 2013: use of a proper RN generator since I saw that quadruplets were duplicated with the pseudo-random number generator used in the lines above

  fLastEntry = (Int_t)(r*numEntries);
  
  std::cout<<"fLastEntry = "<<fLastEntry<<"   r = "<<r<<"   numEntries = "<<numEntries<<std::endl;

  // return the new digit
  fThisDigit = vSeedDigitList.at(fLastEntry);
  cout<<"fThisDigit = "<<fThisDigit<<endl;
  xpos = fDigitX[fThisDigit];
  ypos = fDigitY[fThisDigit];
  zpos = fDigitZ[fThisDigit];
  time = fDigitT[fThisDigit];
  cout<<"xpos = "<<xpos<<"   ypos = "<<ypos<<"   zpos = "<<zpos<<"   time = "<<time<<endl;
  return fThisDigit;
}

// take digit with # digit
int ChooseNextDigit(Double_t& xpos, Double_t& ypos, Double_t& zpos, Double_t& time, int digit)
{
  fThisDigit=digit; // to allow simple copy-paste from random version of ChooseNextDigit

  cout<<"fThisDigit = "<<fThisDigit<<endl;
  xpos = fDigitX[fThisDigit];
  ypos = fDigitY[fThisDigit];
  zpos = fDigitZ[fThisDigit];
  time = fDigitT[fThisDigit];
  cout<<"xpos = "<<xpos<<"   ypos = "<<ypos<<"   zpos = "<<zpos<<"   time = "<<time<<endl;
  return fThisDigit;
}

// Make a quadruple
int ChooseNextQuadruple(Double_t& x0, Double_t& y0, Double_t& z0, Double_t& t0, Double_t& x1, Double_t& y1, Double_t& z1, Double_t& t1, Double_t& x2, Double_t& y2, Double_t& z2, Double_t& t2, Double_t& x3, Double_t& y3, Double_t& z3, Double_t& t3)
{
  int code=0; // 0 -if OK, 1 -if failed to chose 4 different digits 
  int digit0=0;
  int digit1=-1;
  int digit2=-2;
  int digit3=-3;
  int counter1=0;
  int counter2=0;
  int counter3=0;
  digit0 = ChooseNextDigit(x0,y0,z0,t0);
//  digit1 = ChooseNextDigit(x1,y1,z1,t1);
//  digit2 = ChooseNextDigit(x2,y2,z2,t2);
//  digit3 = ChooseNextDigit(x3,y3,z3,t3);

//check that selected digits for the quadruple are not identical
//if they are after 100 attempts then be it, the quadruple will not
//be used later
  while( (digit1<0 || digit1==digit0) && counter1<100)
  {
    digit1 = ChooseNextDigit(x1,y1,z1,t1);
    counter1++;
  }

  while( (digit2<0 || digit2==digit0 || digit2==digit1) && counter2<100)
  {
    digit2 = ChooseNextDigit(x2,y2,z2,t2);
    counter2++;
  }

  while( (digit3<0 || digit3==digit0 || digit3==digit1 || digit3==digit2) && counter3<100)
  {
    digit3 = ChooseNextDigit(x3,y3,z3,t3);
    counter3++;
  }

  if(counter1>=100 || counter2>=100 || counter3>=100) code=1;

  return code;
}

int ChooseNextQuadruple(Double_t& x0, Double_t& y0, Double_t& z0, Double_t& t0, Double_t& x1, Double_t& y1, Double_t& z1, Double_t& t1, Double_t& x2, Double_t& y2, Double_t& z2, Double_t& t2, Double_t& x3, Double_t& y3, Double_t& z3, Double_t& t3, int digit)
{
  ChooseNextDigit(x0,y0,z0,t0,digit);
  ChooseNextDigit(x1,y1,z1,t1,digit+1);
  ChooseNextDigit(x2,y2,z2,t2,digit+2);
  ChooseNextDigit(x3,y3,z3,t3,digit+3);

  return 0;
}

// Calculate NSeedsTarget verticies
int CalcVertexSeeds()
{
  // reset list of seeds
  // ===================
  vSeedVtxX.clear();
  vSeedVtxY.clear();
  vSeedVtxZ.clear();
  vSeedVtxTime.clear();

  Double_t x0 = 0.0;
  Double_t y0 = 0.0;
  Double_t z0 = 0.0;
  Double_t t0 = 0.0;

  Double_t x1 = 0.0;
  Double_t y1 = 0.0;
  Double_t z1 = 0.0;
  Double_t t1 = 0.0;

  Double_t x2 = 0.0;
  Double_t y2 = 0.0;
  Double_t z2 = 0.0;
  Double_t t2 = 0.0;

  Double_t x3 = 0.0;
  Double_t y3 = 0.0;
  Double_t z3 = 0.0;
  Double_t t3 = 0.0;

  double fVtxX1, fVtxY1, fVtxZ1, fVtxTime1;
  double fVtxX2, fVtxY2, fVtxZ2, fVtxTime2;

  int counter=0;  
  cout<<"I'm inside CalcVertexSeeds"<<endl;
//  while( vSeedVtxX.size()<NSeedsTarget && counter<100*NSeedsTarget ) // uncomment for random quadruplets
  for(int i=0;i!=fDigitX.size()-3;i++) // comment for random quadruplets
  {
    cout<<"counter = "<<endl;
    ChooseNextQuadruple(x0,y0,z0,t0,
	                            x1,y1,z1,t1,
        	                    x2,y2,z2,t2,
	                            x3,y3,z3,t3, //); uncomment for random quadruplets
				    i); // comment for random quadruplets
    cout<<"counter = "<<counter<<endl; 
    std::cout << "   digit0: (x,y,z,t)=(" << x0 << "," << y0 << "," << z0 << "," << t0 << ") " << std::endl;
    std::cout << "   digit1: (x,y,z,t)=(" << x1 << "," << y1 << "," << z1 << "," << t1 << ") " << std::endl;
    std::cout << "   digit2: (x,y,z,t)=(" << x2 << "," << y2 << "," << z2 << "," << t2 << ") " << std::endl;
    std::cout << "   digit3: (x,y,z,t)=(" << x3 << "," << y3 << "," << z3 << "," << t3 << ") " << std::endl;

    FindVertex(x0,y0,z0,t0,
                              x1,y1,z1,t1,
                              x2,y2,z2,t2,
                              x3,y3,z3,t3,
                              fVtxX1,fVtxY1,fVtxZ1,fVtxTime1,
                              fVtxX2,fVtxY2,fVtxZ2,fVtxTime2);

     std::cout << "   result1: (x,y,z,t)=(" << fVtxX1 << "," << fVtxY1 << "," << fVtxZ1 << "," << fVtxTime1 << ") " << std::endl
               << "   result2: (x,y,z,t)=(" << fVtxX2 << "," << fVtxY2 << "," << fVtxZ2 << "," << fVtxTime2 << ") " << std::endl;
     cout<<"-------------"<<endl;

    if(fVtxX1==-99999.9 && fVtxX2==-99999.9) //no solutions try anouther quadruple
    { 
      counter++;
      continue;
    }


    bool inside_det;
    if(sqrt(fVtxX1*fVtxX1+fVtxY1*fVtxY1+fVtxZ1*fVtxZ1)<R_SPHERE) inside_det=true; else inside_det=false;
    std::cout<<"Solution1: inside_det = "<<inside_det<<"  X = "<<fVtxX1<<"  Y = "<<fVtxY1<<"  Z = "<<fVtxZ1<<"  T = "<<fVtxTime1<<std::endl;
//    if(!inside_det)
//    {
//      std::cout<<"Solution outside detector"<<std::endl;
//      continue;
//    }

    // add first digit
    if( inside_det ){
      vSeedVtxX.push_back(fVtxX1);
      vSeedVtxY.push_back(fVtxY1);
      vSeedVtxZ.push_back(fVtxZ1);
      vSeedVtxTime.push_back(fVtxTime1);
      std::cout << "New vertex seed 1: x= " << fVtxX1 << ", " << fVtxY1 << ", " << fVtxZ1 << ", " <<fVtxTime1 << std::endl;
    }


    if(sqrt(fVtxX2*fVtxX2+fVtxY2*fVtxY2+fVtxZ2*fVtxZ2)<R_SPHERE) inside_det=true; else inside_det=false;
    std::cout<<"Solution2: inside_det = "<<inside_det<<"  X = "<<fVtxX2<<"  Y = "<<fVtxY2<<"  Z = "<<fVtxZ2<<"  T = "<<fVtxTime2<<std::endl;
//    if(!inside_det)
//    {
//      std::cout<<"Solution outside detector"<<std::endl;
//      continue;
//    }

    // add second digit
    if( inside_det ){
      vSeedVtxX.push_back(fVtxX2);
      vSeedVtxY.push_back(fVtxY2);
      vSeedVtxZ.push_back(fVtxZ2);
      vSeedVtxTime.push_back(fVtxTime2);
      std::cout << "New vertex seed 2: x= " << fVtxX2 << ", " << fVtxY2 << ", " << fVtxZ2 << ", " <<fVtxTime2 << std::endl;
    }

    counter++;

  }

  std::cout << "The number of calculated seeds in vSeedVtxX, vSeedVtxY, vSeedVtxZ, vSeedVtxTime for this event is = " << vSeedVtxX.size() << std::endl;


  return 0;
}

// calculates time residuals for a given vertex (vertex time excluded by default)
int FillResiduals(Double_t vtxX, Double_t vtxY, Double_t vtxZ, Double_t vtxT=0)
{
  fDelta.clear();
  for( Int_t idigit=0; idigit<fDigitX.size(); idigit++ )
  {
    Double_t dx = fDigitX[idigit]-vtxX;
    Double_t dy = fDigitY[idigit]-vtxY;
    Double_t dz = fDigitZ[idigit]-vtxZ;
    Double_t ds = sqrt(dx*dx+dy*dy+dz*dz);
    double fPointResidual = fDigitT[idigit] - ds/(C_VAC/N_REF) - vtxT;
    fDelta.push_back(fPointResidual);
  }
  return 0;
}

// mode==0 -> sorts by time
// mode==1 -> sorts by residuals
int SortDigits(int mode,double X=0, double Y=0, double Z=0, double T=0)
{
  if(mode==1) FillResiduals(X,Y,Z,T);

  double tmpX;
  double tmpY;
  double tmpZ;
  double tmpT;
  double tmpQ;
  double tmpPE;
  double tmpW;
  double tmpV;
  double tmpVgr;
  double tmpN;
  double tmpNgr;
  double tmpDelta;
 
  int flag=1;

  for(int i=0;(i!=fDigitX.size()) && flag;i++)
  {
    flag = 0;
    for(int j=0;j!=fDigitX.size()-1;j++)
    if( (fDigitT[j]>fDigitT[j+1]) && mode==0 ||
	(fDelta[j]>fDelta[j+1]) && mode==1     )
    {
//      cout<<"i = "<<i<<"   j = "<<j<<endl;
      flag=1;

      tmpX = fDigitX[j];
      tmpY = fDigitY[j];
      tmpZ = fDigitZ[j];
      tmpT = fDigitT[j];
      tmpQ = fDigitQ[j];
//      tmpPE = fDigitPE[j];
      tmpW = fDigitW[j];
      tmpV = fDigitV[j];
      tmpVgr = fDigitVgr[j];
      tmpN = fDigitN[j];
      tmpNgr = fDigitNgr[j];
      if(mode==1) tmpDelta=fDelta[j];

      fDigitX[j]=fDigitX[j+1];
      fDigitY[j]=fDigitY[j+1];
      fDigitZ[j]=fDigitZ[j+1];
      fDigitT[j]=fDigitT[j+1];
      fDigitQ[j]=fDigitQ[j+1];
//      fDigitPE[j]=fDigitPE[j+1];
      fDigitW[j]=fDigitW[j+1];
      fDigitV[j]=fDigitV[j+1];
      fDigitVgr[j]=fDigitVgr[j+1];
      fDigitN[j]=fDigitN[j+1];
      fDigitNgr[j]=fDigitNgr[j+1];
      if(mode==1) fDelta[j]=fDelta[j+1];

      fDigitX[j+1]=tmpX;
      fDigitY[j+1]=tmpY;
      fDigitZ[j+1]=tmpZ;
      fDigitT[j+1]=tmpT;
      fDigitQ[j+1]=tmpQ;
//      fDigitPE[j+1]=tmpPE;
      fDigitW[j+1]=tmpW;
      fDigitV[j+1]=tmpV;
      fDigitVgr[j+1]=tmpVgr;
      fDigitN[j+1]=tmpN;
      fDigitNgr[j+1]=tmpNgr;
      if(mode==1) fDelta[j+1]=tmpDelta;
    }
  }
  return 0;
}


// Calculates FOM by looking at time residuals
// make sure fDelta has been filled before calling this function
int TimePropertiesLnL(double & vtx_time, double & fom)
{
  double A = 1.0 / ( 2.0*TSIGMA*sqrt(0.5*TMath::Pi()) );
  double Preal=0.;
  double P=0.;
  fom=0.;
  double chi2=0.;
  double ndof=0.0;

  cout<<"fDigitX.size() = "<<fDigitX.size()<<endl;
  fNDigits=fDigitX.size();
  if(RECO_MODE==0 || RECO_MODE==2) fNDigits=MAX_FIT_DIGITS;
  for( Int_t idigit=0; idigit<fNDigits; idigit++ )
  {
    double delta = fDelta[idigit] - vtx_time;
    Preal = A*exp(-(delta*delta)/(2.0*TSIGMA*TSIGMA));  
//    P = (1.0-Pnoise)*Preal + Pnoise;
//    chi2 += -2.0*log(P);
//    cout<<"delta = "<<delta<<"   Preal = "<<Preal<<"   fDelta[idigit] = "<<fDelta[idigit]<<"   vtx_time = "<<vtx_time<<endl;
    chi2 += -2.0*log(Preal);
    ndof += 1.0;
  }
  cout<<"chi2 = "<<chi2<<"   ndof = "<<ndof<<endl;
  if( ndof>0.0 ){
    fom = fBaseFOM - 1.0*chi2/ndof;
//    fom = chi2/ndof;
  }
 
  return 0;
}

// To be used by Minuit
static void vertex_time_lnl(Int_t&, Double_t*, Double_t& f, Double_t* par, Int_t)
{
  Double_t vtx_time = par[0];
  Double_t fom = 0.0;
//  time_fit_itr();
  TimePropertiesLnL(vtx_time,fom);
  f = -fom; // note: need to maximize this fom
  return;
}

//Minuit manipulation are decribed here
void FitPointTimePropertiesLnL(Double_t& fit_time, Double_t& fom)
{
  Int_t err = 0;
  Int_t flag = 0;

//  Double_t seedTime = meanTime;

//  Double_t fitTime = 0.0;
  Double_t fitTimeErr = 0.0;

  TMinuit* fMinuitTimeFit = new TMinuit();
  fMinuitTimeFit->SetMaxIterations(5000);

  Double_t* arglist = new Double_t[10];
  arglist[0]=1;  // 1: standard minimization
                 // 2: try to improve minimum

  fMinuitTimeFit->mncler();
  fMinuitTimeFit->SetFCN(vertex_time_lnl);
  fMinuitTimeFit->mnexcm("SET STR",arglist,1,err);
  fMinuitTimeFit->mnparm(0,"vtx_time",seedTime,0.1,-100.0,100.0,err);  //Negative times have to be possible since the true time is at zero                 
  flag = fMinuitTimeFit->Migrad();
  fMinuitTimeFit->GetParameter(0,fit_time,fitTimeErr);

  std::cout <<"fitTime = " << fit_time << ", fitTimeErr = " << fitTimeErr << std::endl;

  delete [] arglist;
  delete fMinuitTimeFit;

  return;
}

// Out of all seed verticies selects the one with the best FOM
// Important: the fit time is discarded and solution from FindVertex(...) is used
int SelectBestSeed(int evt_num)
{
  std::ostringstream oss;
  oss<<evt_num;
  string iter_evt = oss.str();
  
  TH1F* ht = (TH1F*) hT->Clone(("hT_"+iter_evt).c_str());
  TH1F* hdt0 = (TH1F*) hDT0->Clone(("hDT0_"+iter_evt).c_str());
  TH1F* hdt = (TH1F*) hDT->Clone(("hDT_"+iter_evt).c_str());

  Int_t bestSeed = -1;
  Double_t bestFOM = -1.0;

  for(int i=0;i!=vSeedVtxX.size();++i)
  {
    // loop over digits
    // ================
    double Swx=0.;
    double Sw=0.;
//    double meanTime=0.;
    fDelta.clear();
    for( Int_t idigit=0; idigit<fDigitX.size(); idigit++ )
    {
      Double_t dx = fDigitX[idigit]-vSeedVtxX[i];
      Double_t dy = fDigitY[idigit]-vSeedVtxY[i];
      Double_t dz = fDigitZ[idigit]-vSeedVtxZ[i];
      Double_t ds = sqrt(dx*dx+dy*dy+dz*dz);
      //need to check what is the proper variable out of the two listed below
      Double_t time0 = fDigitT[idigit] - 0; //this is what was done in WCSim for the JINST paper
      Double_t time = fDigitT[idigit] - vSeedVtxTime[i];
    
      double fPointResidual0 = time0 - ds/(C_VAC/N_REF);
      double fPointResidual = time - ds/(C_VAC/N_REF);
// TEMP test: use true velocity:
//      double fPointResidual0 = time0 - ds/fDigitVgr[idigit];
//      double fPointResidual = time - ds/fDigitVgr[idigit];
//      cout<<"Lambda = "<<fDigitW[idigit]<<"   index_ref = "<<C_VAC/fDigitV[idigit]<<endl;
      hdt0->Fill(fPointResidual0);
      hdt->Fill(fPointResidual);

      fDelta.push_back(fPointResidual0); //this is what was done in WCSim for the JINST paper
//      fDelta.push_back(fPointResidual);
      double weight = 1.0/(TSIGMA*TSIGMA);
      Swx += time*weight; // here is some room for upgrade id TSIGMA is not always the same 
      Sw += weight;
    }
    ht->Fill(vSeedVtxTime[i]);
    meanTime=Swx/Sw;
    seedTime=vSeedVtxTime[i];

    double vtx_time=0.;
    double fom=0.;
    FitPointTimePropertiesLnL(vtx_time,fom); //do time fit
    TimePropertiesLnL(vtx_time,fom); //calculate fom for fitted time value (possibly this can be done during the fit)
    
    if( fom>bestFOM ){
        bestSeed = i;
        bestFOM = fom;
    }
    cout<<"fom = "<<fom<<"   vtx_time = "<<vtx_time<<endl;  
  }
  fFOM.cd();
  ht->Write();
  hdt0->Write();
  hdt->Write();

  return bestSeed;
}

void PointVertexChi2(Double_t vtxX, Double_t vtxY, Double_t vtxZ, Double_t& vtxTime, Double_t& fom)
{
  fDelta.clear();
  fNDigits=fDigitX.size();
  if(RECO_MODE==0 || RECO_MODE==2) fNDigits=MAX_FIT_DIGITS;
  for( Int_t idigit=0; idigit<fNDigits; idigit++ )
  {
    Double_t dx = fDigitX[idigit]-vtxX;
    Double_t dy = fDigitY[idigit]-vtxY;
    Double_t dz = fDigitZ[idigit]-vtxZ;
    Double_t ds = sqrt(dx*dx+dy*dy+dz*dz);
    double fPointResidual0 = fDigitT[idigit] - ds/(C_VAC/N_REF);
    fDelta.push_back(fPointResidual0);
  }
  TimePropertiesLnL(vtxTime,fom);
  return;
}

static void point_vertex_chi2(Int_t&, Double_t*, Double_t& f, Double_t* par, Int_t)
{
  Double_t vtxX     = par[0]; // centimetres
  Double_t vtxY     = par[1];
  Double_t vtxZ     = par[2];
  Double_t vtime    = par[3];

  Double_t fom = 0.0;

  PointVertexChi2(vtxX,vtxY,vtxZ,vtime,fom);

  f = -fom; // note: need to maximize this fom
  return;
}

int FitPointVertexWithMinuit(double & fX, double & fY, double & fZ, double & fT, int evt_num)
{
  Int_t err = 0;
  Int_t flag = 0;

//  Double_t fitXpos = 0.0;
//  Double_t fitYpos = 0.0;
//  Double_t fitZpos = 0.0;

  Double_t fitXposErr = 0.0;
  Double_t fitYposErr = 0.0;
  Double_t fitZposErr = 0.0;
  Double_t fitTErr = 0.0;

  TMinuit* fMinuitPointVertex = new TMinuit();
  fMinuitPointVertex->SetMaxIterations(5000);

  Double_t* arglist = new Double_t[10];
  arglist[0]=1;  // 1: standard minimization
                 // 2: try to improve minimum
                 //
                 //   // re-initialize everything...
  fMinuitPointVertex->mncler();
  fMinuitPointVertex->SetFCN(point_vertex_chi2);
  fMinuitPointVertex->mnexcm("SET STR",arglist,1,err);
  fMinuitPointVertex->mnparm(0,"x",fX,0.1,fX-50.0,fX+50.0,err);
  fMinuitPointVertex->mnparm(1,"y",fY,0.1,fY-50.0,fY+50.0,err);
  fMinuitPointVertex->mnparm(2,"z",fZ,0.1,fZ-50.0,fZ+50.0,err);
  fMinuitPointVertex->mnparm(3,"z",fT,0.05,fT-5.0,fT+5.0,err);

  flag = fMinuitPointVertex->Migrad();
  fMinuitPointVertex->GetParameter(0,fX,fitXposErr);
  fMinuitPointVertex->GetParameter(1,fY,fitYposErr);
  fMinuitPointVertex->GetParameter(2,fZ,fitZposErr);
  fMinuitPointVertex->GetParameter(3,fT,fitTErr);

  delete [] arglist;
  delete fMinuitPointVertex;

  return 0;
}


// Main function. Loops over input TTree, filters photon hits and calls everything above
int LightReco(char* fInputName, char* fOutputName, int fRecoIt=0, char* fFirstRecoName="f.root")
{
  RECO_MODE=fRecoIt;
  //define reconstrution output here
  TFile f_out(fOutputName,"recreate");
  int evt_num=0;
  double recoVtxX;
  double recoVtxY;
  double recoVtxZ;
  double recoVtxTime;
  double bsVtxX;
  double bsVtxY;
  double bsVtxZ;
  double bsVtxTime;
  double trueVtxX;
  double trueVtxY;
  double trueVtxZ;
  int Nphot;
  float maxalpha[100000];
  float true_alpha[100000];
  float minS1[100000];
  float distS2[100000]; //this S2 is for Matt's isochron
  double S0=0.;
  double S1=0.;
  double S2=0.; //this S2 is 2nd moment
  double S3=0.;
  int recoN_che=0;
  int recoN_sci=0;
  int momNpe=0;
  int momNpe_che=0;
  int momNpe_sci=0;
  int momN_che=0;
  int momN_sci=0;
  float recoDT[NMAX_PHOT];
  float momDT[NMAX_PHOT];
  float momDT_che[NMAX_PHOT];
  float momDT_sci[NMAX_PHOT];

  Float_t edep;

  TTree* reco_out_ntuple = new TTree("ntuple","ntuple");
  reco_out_ntuple->Branch("evt_num",&evt_num,"evt_num/I");
  reco_out_ntuple->Branch("recoVtxX",&recoVtxX,"recoVtxX/D");
  reco_out_ntuple->Branch("recoVtxY",&recoVtxY,"recoVtxY/D");
  reco_out_ntuple->Branch("recoVtxZ",&recoVtxZ,"recoVtxZ/D");
  reco_out_ntuple->Branch("recoVtxTime",&recoVtxTime,"recoVtxTime/D");
  reco_out_ntuple->Branch("bsVtxX",&bsVtxX,"bsVtxX/D");
  reco_out_ntuple->Branch("bsVtxY",&bsVtxY,"bsVtxY/D");
  reco_out_ntuple->Branch("bsVtxZ",&bsVtxZ,"bsVtxZ/D");
  reco_out_ntuple->Branch("bsVtxTime",&bsVtxTime,"bsVtxTime/D");
  reco_out_ntuple->Branch("trueVtxX",&trueVtxX,"trueVtxX/D");
  reco_out_ntuple->Branch("trueVtxY",&trueVtxY,"trueVtxY/D");
  reco_out_ntuple->Branch("trueVtxZ",&trueVtxZ,"trueVtxZ/D");
  reco_out_ntuple->Branch("Nphot",&Nphot,"Nphot/I");
  reco_out_ntuple->Branch("maxalpha",maxalpha,"maxalpha[Nphot]/F");
  reco_out_ntuple->Branch("true_alpha",true_alpha,"true_alpha[Nphot]/F");
  reco_out_ntuple->Branch("minS1",minS1,"minS1[Nphot]/F");
  reco_out_ntuple->Branch("distS2",distS2,"distS2[Nphot]/F"); //this S2 is for Matt's isochron
  reco_out_ntuple->Branch("edep",&edep,"edep/F");
  reco_out_ntuple->Branch("S0",&S0,"S0/D");
  reco_out_ntuple->Branch("S1",&S1,"S1/D");
  reco_out_ntuple->Branch("S2",&S2,"S2/D"); //this S2 is 2nd moment
  reco_out_ntuple->Branch("S3",&S3,"S3/D");
  reco_out_ntuple->Branch("recoN_che",&recoN_che,"recoN_che/I");
  reco_out_ntuple->Branch("recoN_sci",&recoN_sci,"recoN_sci/I");
  reco_out_ntuple->Branch("momNpe",&momNpe,"momNpe/I");
  reco_out_ntuple->Branch("momNpe_che",&momNpe_che,"momNpe_che/I");
  reco_out_ntuple->Branch("momNpe_sci",&momNpe_sci,"momNpe_sci/I");
  reco_out_ntuple->Branch("momN_che",&momN_che,"momN_che/I");
  reco_out_ntuple->Branch("momN_sci",&momN_sci,"momN_sci/I");
  reco_out_ntuple->Branch("recoDT",recoDT,"recoDT[Nphot]/F");
  reco_out_ntuple->Branch("momDT",momDT,"momDT[momNpe]/F");
  reco_out_ntuple->Branch("momDT_che",momDT_che,"momDT_che[momNpe_che]/F");
  reco_out_ntuple->Branch("momDT_sci",momDT_sci,"momDT_sci[momNpe_sci]/F");
  //====================

 
  std::vector<double> x_vec; //these 4 vectors are for the Moment Analysis
  std::vector<double> y_vec;
  std::vector<double> z_vec;
  std::vector<double> sl_vec;
  
  //results of the first itteration of reco is defined here
  double X0, Y0, Z0;
//  int NNN=1;
  TTree* recoTree;
  if(fRecoIt==2)
  {
    TFile* recoFile = new TFile(fFirstRecoName);
    recoTree = (TTree*)recoFile->Get("ntuple");
    recoTree->SetBranchAddress("recoVtxX",&X0);
    recoTree->SetBranchAddress("recoVtxY",&Y0);
    recoTree->SetBranchAddress("recoVtxZ",&Z0);
  }
  //=====================

  

  //input from Sphere1 is defined here
  TFile* Sphere1File = new TFile(fInputName);
  TTree* Hits_Tree = (TTree*)Sphere1File->Get("epgTree");

  const int MAX_phot=100000;
  int N_phot_v=0;
  Float_t x_hit_v[MAX_phot];
  Float_t y_hit_v[MAX_phot];
  Float_t z_hit_v[MAX_phot];
  Float_t cos_theta_v[MAX_phot];
  Float_t photon_wavelength_v[MAX_phot];
  Float_t true_time_v[MAX_phot];
  int PE_creation_v[MAX_phot];
  Float_t PE_time_v[MAX_phot];
  Float_t detection_coverage_included_v[MAX_phot];
  Float_t true_time_corrected_v[MAX_phot];
  Float_t PE_time_corrected_v[MAX_phot];
  Float_t cos_theta_reco_v[MAX_phot];
  Float_t theta_reco_v[MAX_phot];
  Float_t phi_reco_v[MAX_phot];
  int process_v[MAX_phot];
  double trueVtxX_v;
  double trueVtxY_v;
  double trueVtxZ_v;

  TBranch  *N_phot_b = 0;
  TBranch  *x_hit_b = 0;
  TBranch  *y_hit_b = 0;
  TBranch  *z_hit_b = 0;
  TBranch  *cos_theta_b = 0;
  TBranch  *photon_wavelength_b = 0;
  TBranch  *true_time_b = 0;
  TBranch  *PE_creation_b = 0;
  TBranch  *PE_time_b = 0;
  TBranch  *true_time_corrected_b = 0;
  TBranch  *PE_time_corrected_b = 0;
  TBranch  *process_b = 0;
  TBranch  *trueVtxX_b = 0;
  TBranch  *trueVtxY_b = 0;
  TBranch  *trueVtxZ_b = 0;


  TBranch  *edep_b=0;

  Hits_Tree->SetBranchAddress("N_phot", &N_phot_v, &N_phot_b);
  Hits_Tree->SetBranchAddress("x_hit", x_hit_v, &x_hit_b);
  Hits_Tree->SetBranchAddress("y_hit", y_hit_v, &y_hit_b);
  Hits_Tree->SetBranchAddress("z_hit", z_hit_v, &z_hit_b);
  Hits_Tree->SetBranchAddress("cos_theta", cos_theta_v, &cos_theta_b);
  Hits_Tree->SetBranchAddress("photon_wavelength", photon_wavelength_v, &photon_wavelength_b);
  Hits_Tree->SetBranchAddress("true_time", true_time_v, &true_time_b);
  Hits_Tree->SetBranchAddress("PE_creation", PE_creation_v, &PE_creation_b);
  Hits_Tree->SetBranchAddress("PE_time", PE_time_v, &PE_time_b);
  Hits_Tree->SetBranchAddress("true_time_corrected", true_time_corrected_v, &true_time_corrected_b);
  Hits_Tree->SetBranchAddress("PE_time_corrected", PE_time_corrected_v, &PE_time_corrected_b);
  Hits_Tree->SetBranchAddress("process", process_v, &process_b);
  Hits_Tree->SetBranchAddress("edep", &edep, &edep_b);

  Hits_Tree->SetBranchAddress("trueVtxX", &trueVtxX_v, &trueVtxX_b);
  Hits_Tree->SetBranchAddress("trueVtxY", &trueVtxY_v, &trueVtxY_b);
  Hits_Tree->SetBranchAddress("trueVtxZ", &trueVtxZ_v, &trueVtxZ_b);
  //===============

  FillIndex("../data/IndexOfRefraction_KamLAND.txt");

  //ready to loop over input and create digits
  int N_Entries_Hits_Tree = Hits_Tree->GetEntries();
  int nev = EVT_NUM<N_Entries_Hits_Tree ? EVT_NUM : N_Entries_Hits_Tree; 
  for(int i=0;i<nev;i++)
  {
    S0=-999;
    S1=-999;
    S2=-999;
    S3=-999;
    recoN_che=0;
    recoN_sci=0;
    momNpe=0;
    momNpe_che=0;
    momNpe_sci=0;
    momN_che=0;
    momN_sci=0;

    fThisDigit=0;
    fDigitX.clear();
    fDigitY.clear();
    fDigitZ.clear();
    fDigitT.clear();
    fDigitW.clear();
    fDigitV.clear();
    fDigitVgr.clear();
    fDigitQ.clear();
    fDigitPE.clear();
    fDigitN.clear();
    fDigitNgr.clear();
    vSeedDigitList.clear();

    x_vec.clear();
    y_vec.clear();
    z_vec.clear();
    sl_vec.clear();

    int currentEvent=i;
    cout<<"currentEvent = "<<currentEvent<<endl;

    Hits_Tree->GetEntry(currentEvent);

    
    bool early_ph_vec[NMAX_PHOT];
    if(fRecoIt==1)
    {
      for(int ii=0;ii!=N_phot_v;++ii)
        early_ph_vec[ii]=0;
      MarkEarlyPhotons(N_phot_v,x_hit_v,y_hit_v,z_hit_v,PE_time_v,process_v,PE_creation_v,  early_ph_vec);
    }
  
//    if(fRecoIt==2) recoTree->GetEntry(currentEvent-100*(NNN-1));
    if(fRecoIt==2) recoTree->GetEntry(currentEvent);

    for(int iphot=0;iphot!=N_phot_v;iphot++)
    {
//         if(process_v==0) continue; //!Sphere1
//        if(PE_creation_v[iphot]==0) continue; //!Sphere1

      if(fRecoIt==0)
      {
//	if(process_v[iphot]==1) continue;
//        if(photon_wavelength_v[iphot]<413.||photon_wavelength_v[iphot]>448.) continue;	
	if(PE_creation_v[iphot]==0) continue;
	if(PE_time_v[iphot]>34.) continue;
      }
       
      if(fRecoIt==1)
        if(!early_ph_vec[iphot]) continue;
      
      if(fRecoIt==2)
      {
        if(PE_creation_v[iphot]==0) continue;
        double distL = TMath::Sqrt(TMath::Power((x_hit_v[iphot] - X0*10), 2) + (y_hit_v[iphot]-Y0*10)*(y_hit_v[iphot]-Y0*10) + (z_hit_v[iphot]-Z0*10)*(z_hit_v[iphot]-Z0*10))/10;
        double light_vel = C_VAC/N_REF;
        double TPredicted = distL/light_vel;
        recoDT[iphot] = (PE_time_v[iphot] - TPredicted);
        if((PE_time_v[iphot] - TPredicted) > RECO_DT) continue;
      } 
       
//first thing: count how many cherenkov and scintillation photons are used by the reco
//
	if(process_v[iphot]==1) recoN_che++;
        if(process_v[iphot]==0) recoN_sci++;

        fDigitX.push_back(x_hit_v[iphot]/10.); //!Sphere1
        fDigitY.push_back(y_hit_v[iphot]/10.);//!Sphere1
        fDigitZ.push_back(z_hit_v[iphot]/10.);//!Sphere1
	fDigitT.push_back(PE_time_v[iphot]);// - min_PE_time; //!Sphere1
	fDigitQ.push_back(1);
        fDigitW.push_back(photon_wavelength_v[iphot]);
        fDigitN.push_back(INDEX[(int)photon_wavelength_v[iphot]]);
        fDigitNgr.push_back(INDEX_GR[(int)photon_wavelength_v[iphot]]);
	fDigitV.push_back(C_VAC/INDEX[(int)photon_wavelength_v[iphot]]);
        fDigitVgr.push_back(C_VAC/INDEX_GR[(int)photon_wavelength_v[iphot]]);
	vSeedDigitList.push_back(fThisDigit);
        fThisDigit++;
    }

    cout<<"Photon filtering for event #"<<i<<" has just finished. fDigits are ready."<<endl;

    if(fRecoIt==0)  
      SortDigits(0); //sort by time of arrival
    if(fRecoIt==2)
      SortDigits(1,X0,Y0,Z0); //sort by time residuals for a given vertex


    cout<<"NDigits = "<<fDigitX.size()<<endl;

// the double loop below (commented) is just a check whether time ordering worked properly
/*    for(int i=0;i!=fDigitX.size()-1;i++)
    {
      cout<<fDigitT[i]<<"	";
      for(int j=i+1;j!=fDigitX.size();j++)
      {
        if(fDigitT[i]>fDigitT[j])
        {
	  cout<<endl;
          cout<<"WARNING: wrong time ordering"<<endl;
        }
      }
    }
*/

    if(fRecoIt>=10) //option 10 doesn't do any real reconstructionm just smear true vertex with resolution
    {
      recoVtxX = trueVtxX_v + rndVtx.Gaus(VTX_SHIFT_X,VTX_SMEAR);
      recoVtxY = trueVtxY_v + rndVtx.Gaus(VTX_SHIFT_Y,VTX_SMEAR);
      recoVtxZ = trueVtxZ_v + rndVtx.Gaus(VTX_SHIFT_Z,VTX_SMEAR);


      trueVtxX = trueVtxX_v;
      trueVtxY = trueVtxY_v;
      trueVtxZ = trueVtxZ_v;
      Nphot = fDigitX.size();
    } else
    {

  
    //calculate seed vertices
    CalcVertexSeeds();

    //now select the best vertex
    int best_seed = SelectBestSeed(currentEvent);
    cout<<"best_seed = "<<best_seed<<"  vSeedVtxX.size() = "<<vSeedVtxX.size()<<endl;
    //and save reco data
    recoVtxX = vSeedVtxX[best_seed];
    recoVtxY = vSeedVtxY[best_seed];
    recoVtxZ = vSeedVtxZ[best_seed];
    recoVtxTime = vSeedVtxTime[best_seed];

    bsVtxX = vSeedVtxX[best_seed];
    bsVtxY = vSeedVtxY[best_seed];
    bsVtxZ = vSeedVtxZ[best_seed];
    bsVtxTime = vSeedVtxTime[best_seed];

    FitPointVertexWithMinuit(recoVtxX, recoVtxY, recoVtxZ, recoVtxTime,currentEvent);

    trueVtxX = trueVtxX_v;
    trueVtxY = trueVtxY_v;
    trueVtxZ = trueVtxZ_v;
    Nphot = fDigitX.size();
    }

    if(fRecoIt!=2 && fRecoIt<10)
    {
      reco_out_ntuple->Fill();
      evt_num++;
      continue;
    }
//    reco_out_ntuple->Fill();



// lines below are Moment Analysis
// they are using reconstructed vertex but do not affect
// the vertex reconstruction itself
// comment if you only do vertex recon

    for(int iphot=0;iphot!=N_phot_v;iphot++)
    {
//      if(process_v[iphot]==0) continue; //che Light only
//      if(process_v[iphot]==1) continue; //sci Light only

      if(PE_creation_v[iphot]==0) continue;

      double distL = TMath::Sqrt(TMath::Power((x_hit_v[iphot] - recoVtxX*10.), 2) + (y_hit_v[iphot]-recoVtxY*10.)*(y_hit_v[iphot]-recoVtxY*10.) + (z_hit_v[iphot]-recoVtxZ*10.)*(z_hit_v[iphot]-recoVtxZ*10.))/10.; //distance in cm
//      double distL = TMath::Sqrt(TMath::Power((x_hit_v[iphot] - trueVtxX*10.), 2) + (y_hit_v[iphot]-trueVtxY*10.)*(y_hit_v[iphot]-trueVtxY*10.) + (z_hit_v[iphot]-trueVtxZ*10.)*(z_hit_v[iphot]-trueVtxZ*10.))/10.; //distance in cm
      double light_vel = C_VAC/N_REF;
      double TPredicted = distL/light_vel;
//      cout<<"distL = "<<distL<<"   light_vel = "<<light_vel<<endl;

      momDT[momNpe] = (PE_time_v[iphot] - TPredicted);

      if(process_v[iphot]==1) 
      {
	momDT_che[momNpe_che] = momDT[momNpe];
	momNpe_che++;
      }
      if(process_v[iphot]==0) 
      {
	momDT_sci[momNpe_sci] = momDT[momNpe];
	momNpe_sci++;
      }
      momNpe++;

      if((PE_time_v[iphot] - TPredicted) > MOM_DT) continue; 

//count how many cherenkov and scintillation photons are used by the moment analysis
      if(process_v[iphot]==1) momN_che++;
      if(process_v[iphot]==0) momN_sci++;


      //coordinat transfer: draw a sphere around vertex and find intersaction with
      //a vector A, where in old coordinates A=HIT-VTX
      // ALL THE FOLLOWING IN METERS
      distL=distL/100; //distL is already in cm
      x_vec.push_back(6.5*(x_hit_v[iphot]/1000.-recoVtxX/100.)/distL);
      y_vec.push_back(6.5*(y_hit_v[iphot]/1000.-recoVtxY/100.)/distL);
      z_vec.push_back(6.5*(z_hit_v[iphot]/1000.-recoVtxZ/100.)/distL);
      
    }
    cout<<"recoVtxX = "<<recoVtxX<<"   trueVtxX = "<<trueVtxX<<endl;
    cout<<"recoVtxY = "<<recoVtxY<<"   trueVtxY = "<<trueVtxY<<endl;
    cout<<"recoVtxZ = "<<recoVtxZ<<"   trueVtxZ = "<<trueVtxZ<<endl;
   
    cout<<"x_vec.size() = "<<x_vec.size()<<endl;
    cout<<"x[0] = "<<x_vec[0]<<"   y[0] = "<<y_vec[0]<<"   z[0] = "<<z_vec[0]<<endl;
    cout<<"x[10] = "<<x_vec[10]<<"   y[10] = "<<y_vec[10]<<"   z[10] = "<<z_vec[10]<<endl;
    analysis(x_vec, y_vec, z_vec, evt_num, sl_vec);
    S0 = sl_vec[0];
    S1 = sl_vec[1];
    S2 = sl_vec[2];
    S3 = sl_vec[3];

/*
// lines below are for testing Matt's Isochron
// ignore them - no effect on vertex reconstruction
    for(int it=0; it!=fDigitX.size();++it)
    { 
      TVector3 fVec(fDigitX[it]-recoVtxX,fDigitY[it]-recoVtxY,fDigitZ[it]-recoVtxZ);
//      TVector3 fVec(fDigitX[it]-0,fDigitY[it]-0,fDigitZ[it]-0);
      double dD = fVec.Mag();
      double dT = fDigitT[it] - recoVtxTime;
      double n_eff = N_REF;//fDigitNgr[it];//C_VAC/(dD/dT);
      cout<<"dD = "<<dD<<"   dT = "<<dT<<"   n_eff = "<<n_eff<<"   C_VAC/(dD/dT)"<<C_VAC/(dD/dT)<<endl;
      maxalpha[it] = (float) CalcMaxAlpha(dD,dT,n_eff);//INDEX[(int)fDigitW[it]]);//N_REF);
      TVector3 fXaxis(1,0,0);
      true_alpha[it] = fXaxis.Angle(fVec);
      double alpha = maxalpha[it];//true_alpha[it];
      minS1[it] = (float) CalcS1(alpha,dD,dT,n_eff);//INDEX[(int)fDigitW[it]]);//N_REF);
      distS2[it] = (float) CalcS2(alpha,dD,dT,n_eff);//INDEX[(int)fDigitW[it]]);//N_REF);
    }
//========= end of Isochron test =============
*/
    reco_out_ntuple->Fill();
    evt_num++;

  } //end i-loop over Hits_Tree entries

  f_out.cd();
  reco_out_ntuple->Write();
  return 0;
}

//end
