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

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>

int EVT_NUM=100; //controls maximum number of events to be processed
double R_SPHERE=650; //sphere diameter [cm]
double N_REF=1.53; //average index of refraction
double C_VAC=29.9792458; //speed of light in vacuum [cm/ns]
int NSeedsTarget=400; //number of quadruplets
double TSIGMA=0.5; //total time spread (including detector TTS chromatic dispersions)
static const int NMAX_PHOT=100000; //
static const int NPHI=12;
static const int NTHETA=12;


map<int, double> INDEX;

#include "help_func.C"

double fBaseFOM=100.0; //Figure of merit. Borrowed from WCSim: the higher it is the better
double meanTime=0.;
double seedTime=0.;

TRandom RND;

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
  while( vSeedVtxX.size()<NSeedsTarget && counter<100*NSeedsTarget )
  {
    cout<<"counter = "<<endl;
    ChooseNextQuadruple(x0,y0,z0,t0,
	                            x1,y1,z1,t1,
        	                    x2,y2,z2,t2,
	                            x3,y3,z3,t3);
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
  for( Int_t idigit=0; idigit<fDigitX.size(); idigit++ )
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
    
//      double fPointResidual0 = time0 - ds/(C_VAC/N_REF);
//      double fPointResidual = time - ds/(C_VAC/N_REF);
// TEMP test: use true velocity:
      double fPointResidual0 = time0 - ds/fDigitV[idigit];
      double fPointResidual = time - ds/fDigitV[idigit];
//      cout<<"Lambda = "<<fDigitW[idigit]<<"   index_ref = "<<C_VAC/fDigitV[idigit]<<endl;
      hdt0->Fill(fPointResidual0);
      hdt->Fill(fPointResidual);

      fDelta.push_back(fPointResidual0);
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

// Main function. Loops over input TTree, filters photon hits and calls everything above
int LightReco(char* fInputName, char* fOutputName, int fRecoIt=0, char* fFirstRecoName="f.root")
{
  //define reconstrution output here
  TFile f_out(fOutputName,"recreate");
  double recoVtxX;
  double recoVtxY;
  double recoVtxZ;
  double recoVtxTime;
  TTree* reco_out_ntuple = new TTree("ntuple","ntuple");
  reco_out_ntuple->Branch("recoVtxX",&recoVtxX,"recoVtxX/D");
  reco_out_ntuple->Branch("recoVtxY",&recoVtxY,"recoVtxY/D");
  reco_out_ntuple->Branch("recoVtxZ",&recoVtxZ,"recoVtxZ/D");
  reco_out_ntuple->Branch("recoVtxTime",&recoVtxTime,"recoVtxTime/D");
  //====================
  
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
  //===============

  FillIndex("../data/IndexOfRefraction_KamLAND.txt");

  //ready to loop over input and create digits
  int N_Entries_Hits_Tree = Hits_Tree->GetEntries();
  int nev = EVT_NUM<N_Entries_Hits_Tree ? EVT_NUM : N_Entries_Hits_Tree; 
  for(int i=0;i<nev;i++)
  {

    fThisDigit=0;
    fDigitX.clear();
    fDigitY.clear();
    fDigitZ.clear();
    fDigitT.clear();
    fDigitW.clear();
    fDigitV.clear();
    fDigitQ.clear();
    fDigitPE.clear();
    vSeedDigitList.clear();

    int currentEvent=i;
    cout<<"currentEvent = "<<currentEvent<<endl;

    Hits_Tree->GetEntry(currentEvent);

    
    bool early_ph_vec[NMAX_PHOT];
    if(fRecoIt==1)
    {
      for(int i=0;i!=N_phot_v;++i)
        early_ph_vec[i]=0;
      MarkEarlyPhotons(N_phot_v,x_hit_v,y_hit_v,z_hit_v,PE_time_v,process_v,PE_creation_v,  early_ph_vec);
    }
  
//    if(fRecoIt==2) recoTree->GetEntry(currentEvent-100*(NNN-1));
    if(fRecoIt==2) recoTree->GetEntry(currentEvent);

    for(int iphot=0;iphot!=N_phot_v;iphot++)
    {
//         if(process_v==0) continue; //!Sphere1
//         //         if(PE_creation_v[iphot]==0) continue; //!Sphere1

      if(fRecoIt==0)
      {
	//if(process_v[iphot]==0) continue;
	if(PE_creation_v[iphot]==0) continue;
	if(PE_time_v[iphot]>34.) continue;
      }
       
      if(fRecoIt==1)
        if(!early_ph_vec[iphot]) continue;
      
      if(fRecoIt==2)
      {
        double distL = TMath::Sqrt(TMath::Power((x_hit_v[iphot] - X0*10), 2) + (y_hit_v[iphot]-Y0*10)*(y_hit_v[iphot]-Y0*10) + (z_hit_v[iphot]-Z0*10)*(z_hit_v[iphot]-Z0*10));
        double light_vel = 300./1.53;
        double TPredicted = distL/light_vel;
        if((PE_time_v[iphot] - TPredicted) > 3.0) continue;
      } 
        fDigitX.push_back(x_hit_v[iphot]/10.); //!Sphere1
        fDigitY.push_back(y_hit_v[iphot]/10.);//!Sphere1
        fDigitZ.push_back(z_hit_v[iphot]/10.);//!Sphere1
	fDigitT.push_back(PE_time_v[iphot]);// - min_PE_time; //!Sphere1
	fDigitQ.push_back(1);
        fDigitW.push_back(photon_wavelength_v[iphot]);
//	int lambda = (int)
	fDigitV.push_back(C_VAC/INDEX[(int)photon_wavelength_v[iphot]]);
	vSeedDigitList.push_back(fThisDigit);
        fThisDigit++;
    }

    cout<<"Photon filtering for event #"<<i<<" has just finished. fDigits are ready."<<endl;
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
    reco_out_ntuple->Fill();

  } //end i-loop over Hits_Tree entries

  f_out.cd();
  reco_out_ntuple->Write();
  return 0;
}

