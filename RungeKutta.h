#ifndef _RUNGE_KUTTA_H_
#define _RUNGE_KUTTA_H_
#include "TimeIntegrator.h"

class RungeKutta:public TimeIntegrator{
 public:
  RungeKutta()=default;
 RungeKutta(const Vec<double,6>& _u,std::ofstream & ofs)
    {
      U=_u;
      ofs<<U[0]<<'\t'<<U[1]<<std::endl;
    }

 void init(const Vec<double,6>& _u,const ScalarFunction & f,std::ofstream & ofs,int p=1)
    {
      U=_u;
      ofs<<U[0]<<'\t'<<U[1]<<std::endl;
    }
  
  Vec<double,6> solve(std::ofstream & ofs,const ScalarFunction & f,const int _N,const double T,int p=1){
    double k=T/_N;
    Vec<double,6> temp1,temp2,temp3,temp4;
    for(int i=0;i<_N;i++){
      temp1=f(U);
      temp2=f(U+temp1*(k/2));
      temp3=f(U+temp2*(k/2));
      temp4=f(U+temp3*k);
      U=U+(temp1+temp2*2+temp3*2+temp4)*(k/6);
      U[2]=0.0;U[5]=0.0;
      ofs<<U[0]<<'\t'<<U[1]<<std::endl;
    }
    return U;
  };
  static TimeIntegrator* Creator(){
    return new RungeKutta();
  }
};


#endif
