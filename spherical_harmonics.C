//Validation 
// Ylm's are verified by http://keisan.casio.com/exec/system/1180573409
// for l=0,1,2,3,4 and m=0 at theta=30 and phi=45
// Plm's are verified for l=1..4, m=1..l
//

#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_legendre.h"

#include "TMath.h"
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TVector3.h"

#include "globals.h"

#include <iostream>
#include <iomanip> 
#include <string>


using namespace std;


/*
struct th2_info
{
  int Nentries;
  int bin_num_x;
  int bin_num_y;
  double bin_size_x;
  double bin_size_y;
  double min_x;
  double min_y;
  double max_x;
  double max_y;

};

TH2F* hThetaPhi;
TH2F* htp;
int NL=11;

th2_info h2_par;
*/
double Nlm(int l, int m)
{
  double f1 = (2*l+1)/(4*TMath::Pi());
  double f2 = TMath::Factorial(l-m)/TMath::Factorial(l+m);
  return sqrt(f1*f2);
}

double Plm(int l, int m, double x)
{

  if(l<0 || l>4 || m<0)
  {
    cout<<"WARNING: Plm is currently not defined for l="<<l<<"  and m="<<m<<endl;
    return 0;
  } 
  
  if(l==0 && m==0)
    return 1;

  if(l==1 && m==0)
    return x;
  if(l==1 && m==1)
    return -sqrt(1-x*x);

  if(l==2 && m==0)
    return 0.5*(3*x*x-1);
  if(l==2 && m==1)
    return -3*x*sqrt(1-x*x) ;
  if(l==2 && m==2)
    return 3*(1-x*x);

  if(l==3 && m==0)
    return 0.5*(5*x*x*x-3*x);
  if(l==3 && m==1)
    return -1.5*(5*x*x-1)*sqrt(1-x*x);
  if(l==3 && m==2)
    return 15*x*(1-x*x);
  if(l==3 && m==3)
    return -15*(1-x*x)*sqrt(1-x*x);

  if(l==4 && m==0)
    return (1./8.)*(35*x*x*x*x-30*x*x+3);
  if(l==4 && m==1)
    return -(5./2.)*(7*x*x*x-3*x)*sqrt(1-x*x);
  if(l==4 && m==2)
    return (15./2.)*(7*x*x-1)*(1-x*x);
  if(l==4 && m==3)
    return -105*x*(1-x*x)*sqrt(1-x*x);
  if(l==4 && m==4)
    return 105*(1-x*x)*(1-x*x);

/*  if(l==5 && m==0) return P50(x);
  if(l==5 && m==0) return P51(x);
  if(l==5 && m==0) return P52(x);
  if(l==5 && m==0) return P53(x);
  if(l==5 && m==0) return P54(x);
  if(l==5 && m==0) return P55(x);
*/
}

double PlmCos(int l, int m, double theta)
{
  return Plm(l,m,TMath::Cos(theta));
}

double PlmCos_gsl(int l, int m, double theta)
{
  return gsl_sf_legendre_Plm(l,m,TMath::Cos(theta));
}

double Ylm(int l, int m, double theta, double phi)
{
  if(m==0)
    return Nlm(l,m)*PlmCos_gsl(l,m,theta);
  if(m>0)
    return sqrt(2.)*Nlm(l,m)*PlmCos_gsl(l,m,theta)*TMath::Cos(m*phi);
  if(m<0)
    return sqrt(2.)*Nlm(l,abs(m))*PlmCos_gsl(l,abs(m),theta)*TMath::Sin(abs(m)*phi);
}


th2_info create_h2_par(TH2F* h2)
{
  th2_info h2_par;
  h2_par.Nentries = (int)h2->Integral();
  h2_par.bin_num_x = h2->GetNbinsX();
  h2_par.bin_num_y = h2->GetNbinsY();
  h2_par.bin_size_x = (double)h2->GetXaxis()->GetBinWidth(1);
  h2_par.bin_size_y = (double)h2->GetYaxis()->GetBinWidth(1);
  h2_par.min_x = (double)h2->GetXaxis()->GetBinLowEdge(1);
  h2_par.min_y = (double)h2->GetYaxis()->GetBinLowEdge(1);
  h2_par.max_x = (double)h2->GetXaxis()->GetBinLowEdge(h2_par.bin_num_x) + h2_par.bin_size_x;
  h2_par.max_y = (double)h2->GetYaxis()->GetBinLowEdge(h2_par.bin_num_y) + h2_par.bin_size_y;
  return h2_par;
}


double Func(double theta, double phi, int l=0, int m=0)
{
//  double pi=TMath::Pi();
//  return (theta*theta+pi/3.)*(phi+pi/2.);
/*  if(theta<0.8*TMath::Pi())
  {
  //  return 1;
    if(phi<=TMath::Pi())
      return 1;
    if(phi>TMath::Pi())
      return 0;
  } else
  {
    return 0;
  }
*/
  int bin_x = int((theta-h2_par.min_x)/h2_par.bin_size_x)+1; //bins count from 1
  int bin_y = int((phi-h2_par.min_y)/h2_par.bin_size_y)+1; //bins count from 1
//  double pdf = hThetaPhi->GetBinContent(bin_x,bin_y)/(h2_par.Nentries*h2_par.bin_size_x*h2_par.bin_size_y);  
  double pdf = htp->GetBinContent(bin_x,bin_y)/(h2_par.Nentries*h2_par.bin_size_x*h2_par.bin_size_y);
  return pdf;

}

double I2D(double (*func)(double,double,int,int), double x1, double x2, double y1, double y2, int l=0, int m=0)
{
  TRandom rndGen;
  double sum=0.;
  double sum2=0.;
  double I=0.;
  double dI=0.;
  double err=0.;//relative error:  err=dI/I
  double NUM=0.;
  double dx=x2-x1;
  double dy=y2-y1;

  for(int i=1;i!=10000000;i++)
  {
    double x=x1+rndGen.Rndm()*dx;
    double y=y1+rndGen.Rndm()*dy;
    double val=func(x,y,l,m)*dx*dy;
    sum+=val;
    sum2+=val*val;
    NUM=double(i);
    I=sum/NUM;
    dI=sqrt((sum2/NUM-sum*sum/(NUM*NUM))/NUM);
    err=dI/I;
    if(fabs(I)<0.001 && i>1000)
    {
      I=0.000001;
      dI=0.000001;
      err=0.000001;
    }
    if(fabs(err)<0.04 && i>1000) break;
  }
  cout<<"I = "<<I<<"   dI = "<<dI<<"   err = "<<err<<"   N = "<<NUM<<endl;
  if(fabs(err)>0.05) 
    cout<<"WARNING: integral does not converge!"<<endl;
  return I;
}

double SinFunc(double theta, double phi, int l=0, int m=0 )
{
  return Func(theta,phi)*TMath::Sin(theta);
}

double SinFunc2(double theta, double phi, int l=0, int m=0)
{
  return SinFunc(theta,phi)*Func(theta,phi);
}

double SinFuncYlm(double theta, double phi, int l=0, int m=0)
{
  return SinFunc(theta,phi)*Ylm(l,m,theta,phi);
}

double SurfI()
{
  return I2D(&SinFunc, 0.,TMath::Pi(),0.,2.*TMath::Pi());
}

double alm(int l, int m)
{
  double dI=0;
  return I2D(&SinFuncYlm, 0.,TMath::Pi(),0.,2.*TMath::Pi(),l,m);
}

double Fl(int l, double theta, double phi)
{
  cout<<"enter Fl:"<<endl;
  double sum=0.;
  for(int m=-l;m<=l;m++)
  {
    double a=alm(l,m);
    cout<<"l = "<<l<<"    m = "<<m<<"   alm = "<<a<<endl;
    sum+=a*Ylm(l,m,theta,phi);
  }
  return sum;
}

double Fl2(double theta, double phi, int l)
{
  cout<<"enter Fl2:"<<endl;
  double fl_val=Fl(l,theta,phi);
  return fl_val*fl_val;
}

double SinFl2(double theta, double phi, int l, int m=0)
{
  cout<<"enter SinFl2:"<<endl;
  return TMath::Sin(theta)*Fl2(theta,phi,l);
}

double FlNorm(int l)
{
  double dI=0;
  return I2D(&SinFl2, 0.,TMath::Pi(),0.,2.*TMath::Pi(),l);
}

double Sl(int l)
{
  double sum=0.;
  for(int m=-l;m<=l;m++)
  {
    double a=alm(l,m);
    cout<<"l = "<<l<<"    m = "<<m<<"   alm = "<<a<<endl;
    sum+=a*a;
  }
  return sum;// /(2.*l+1);//*4.*TMath::Pi()/(2.*l+1);  
}

double SumSl(int b)
{
  double sum=0.;
  for(int l=0;l<=b;l++)
  {
    double sl=Sl(l);
    cout<<"	S("<<l<<") = "<<sl<<endl;
    sum+=sl;
  }
  return sum;
}

double FuncNorm()
{
  double dI=0;
  return I2D(&SinFunc2, 0.,TMath::Pi(),0.,2.*TMath::Pi());
}

int analysis()
{
  string case_name="all_pXmX_100_10_20";
  TFile f_in(Form("fThetaPhi_%s.root",case_name.c_str()));
//  TFile f_in("fBkg_all.root");
  hThetaPhi = (TH2F*)f_in.Get("hThetaPhi");
  h2_par = create_h2_par(hThetaPhi);


  TH1F* hL = new TH1F("hL","hL",11,0,11);
  hL->SetName(Form("hL_%s",case_name.c_str()));

  double sum=0.;
  for(int l=0;l!=11;l++)
  {
    double sl=Sl(l);
    cout<<"     S("<<l<<") = "<<sl<<endl;
    sum+=sl;
    hL->SetBinContent(l+1,sl);
  }

  TCanvas c;
  c.cd();
  hL->Draw();
  //c.Print("hL.pdf");
  c.Print(Form("hL_%s.pdf",case_name.c_str()));
  TFile f_out("fL_simple_.root","UPDATE");
  f_out.cd();
  hL->Write();
  f_out.Close();

  double norm=FuncNorm();
  cout<<"Norm = "<<norm<<"   SumSl = "<<sum<<"   epsilon = "<<(norm-sum)/norm<<endl;
  return 0;  

} 


int analysis(vector<double> & x, vector<double> & y, vector<double> & z, int evt, vector<double> & fSL)
{
  //string case_name="bkgTl208_center_allLight_35ns_100";
  //string case_name="Se82_center_allLight_33ns_1k"
  string case_name="tmp";//"1el_2p995MeV_center_allLight_33ns_400";

  TH3F* hxyz = new TH3F(Form("hxyz_%d",evt),Form("hxyz_%d",evt),70,-7,7, 70,-7,7, 70,-7,7);
  htp = new TH2F(Form("htp_%d",evt),Form("htp_%d",evt),Ntheta,0,3.1415926535, Nphi,0,6.283185307);
  for(int i=0;i!=x.size();i++)
  {
    TVector3 hit_vec(x[i],y[i],z[i]);
    double theta = hit_vec.Theta();
    double phi = hit_vec.Phi();
    if(phi<0.) phi=2*pi+phi;
    htp->Fill(theta,phi);
    hxyz->Fill(x[i],y[i],z[i]);
  }
  h2_par = create_h2_par(htp);
  cout<<"htp_Nev = "<<h2_par.Nentries<<endl;
  cout<<"NL = "<<NL<<endl;

  TH1F* hL = new TH1F(Form("hL_%d",evt),Form("hL_%d",evt),NL,0,NL);

  double sum=0.;
  for(int l=0;l!=NL;l++)
  {
    double sl=Sl(l);
    cout<<"     S("<<l<<") = "<<sl<<endl;
    sum+=sl;
    hL->SetBinContent(l+1,sl);
    fSL.push_back(sl);
  }

  string fOption;
  if(evt==0) fOption="RECREATE";
  else fOption="UPDATE";

  TFile f_out(Form("fL_%s.root",case_name.c_str()),fOption.c_str());
  f_out.cd();
  hL->Write();
  htp->Write();
  hxyz->Write();
  f_out.Close();
  delete hL;
  delete htp;

  return 0;
}






