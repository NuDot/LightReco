#include "TH2F.h"

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

const double pi=3.1415926535;

int NL=5;
//int NL=6;
const int Ntheta=10;//7;//14;//10;
const int Nphi=20;//14;//28;//20;


th2_info h2_par;
 
int analysis(vector<double> &, vector<double> &, vector<double> &, int, vector<double> &, double &, char*);
int fourie_analysis(vector<double> &, vector<double> &, vector<double> &, int, vector<double> &);
double Psi1(vector<double> &, vector<double> &, vector<double> &);
double Psi2(vector<double> &, vector<double> &, vector<double> &);
double Psi3(vector<double> &, vector<double> &, vector<double> &);

