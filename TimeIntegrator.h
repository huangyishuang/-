#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_
#include<fstream>
#include "vect.h"
#include "lapacke.h"
#include<vector>
#include "ScalarFunction.h"

class TimeIntegrator{
 protected:
  Vec<double,6> U;
  Vec<double,6> F;
 public:
  TimeIntegrator()=default;
  virtual void init(const Vec<double,6>& _u,const ScalarFunction & f,st                                                                                                                                                                                                            d::ofstream & ofs,int p=1)=0;
  virtual Vec<double,6> solve(std::ofstream & ofs,const ScalarFunction & f,const int N,const double T,int p=1)=0;
};

class LMM : public TimeIntegrator{
 public:
  LMM()=default;
  
  virtual void onestep(std::ofstream & ofs,const ScalarFunction & f,const double k,const int p)=0;

  Vec<double,6> Nel;opkl,.ton(Vec<double,6> u, Vec<double,6> v,const ScalarFunction & VF, const ScalarFunction &J, double k, double beta_s)
{
  double epsilon = 1e-10;
  int M = 1;
  const int m=6;
  Vec<double,6> u0 = u, v0 = v;
  Vec<double,6> _f = u0 - (VF(u0) * (k * beta_s)) - v;
  Vec<double,36> _j;
  double norm = _f.maxNorm();
  double f[6];
  double j1[36];
  int ipiv[6];
  while((norm > epsilon) && (M < 1000)){
    _j = J.getJ(u0);
    for(int i = 0; i < m; i++){
      f[i] = _f[i]; 
      for(int j = 0; j < m; j++){
	j1[i*m+j] =_j[i*m+j] * (-k) * beta_s;
      }
      j1[i*(m+1)] = 1 + j1[i*(m+1)];
    }
    LAPACKE_dgesv(LAPACK_COL_MAJOR, m, 1, j1, m, ipiv, f, m);
    
    for(int i = 0; i < m; i++)
      u0[i] = u0[i] - f[i];
    
    _f = u0 - (VF(u0) * k * beta_s) - v;
    
    norm = _f.maxNorm();
    M++;
  }
  return u0;
}


};

#endif
