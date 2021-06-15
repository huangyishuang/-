#ifndef _ADAMSMOULTON_H_
#define _ADAMSMOULTON_H_
#include<fstream>
#include "vect.h"
#include "ScalarFunction.h"
#include "TimeIntegrator.h"
#include "lapacke.h"
//#include "AMcoef.h"
#include<cmath>


class AdamsMoulton : public LMM{
 protected:
  std::vector<Vec<double,6>> Fcontainer;
  std::vector<std::vector<double> > coe={{0},{1.0/2,1.0/2},{-1.0/12,8.0/12,5.0/12},{1.0/24,-5.0/24,19.0/24,9.0/24},{-19.0/720,106.0/720,-264.0/720,646.0/720,251.0/720}};
 public:
  AdamsMoulton()=default;

  /*AdamsMoulton(const Vec<double,6>& _u,const ScalarFunction & f,std::ofstream & ofs,const int p)
    {
      U=_u;
      F=f(U);
      Fcontainer.resize(p);
      Fcontainer[0]=F;
      ofs<<U[0]<<'\t'<<U[1]<<std::endl;
      }*/

    void init(const Vec<double,6>& _u,const ScalarFunction & f,std::ofstream & ofs,int p=1)
    {
      U=_u;
      F=f(U);
      Fcontainer.erase(Fcontainer.begin(),Fcontainer.end());
      //Fcontainer.resize(p);
      //Fcontainer[0]=F;
      Fcontainer.push_back(F);
      ofs<<U[0]<<'\t'<<U[1]<<std::endl;
    }
    
  void onestep(std::ofstream & ofs,const ScalarFunction & f,const double k,const int p){
    for(int i=1;i<p;i++){
      U=U+F*k;
      F=f(U);
      //Fcontainer[i]=F;
      Fcontainer.push_back(F);
      ofs<<U[0]<<'\t'<<U[1]<<std::endl;
      U[2]=0.0;U[5]=0.0;
    }
  };


  Vec<double,6> solve(std::ofstream & ofs,const ScalarFunction & f,const int _N,const double T,int p=1){
    double k=T/_N;
    onestep(ofs,f,k,p);
    int N=_N-p+1;
    Jacobi J;
    for(int i=0;i<N;i++){
      Vec<double,6> sum;
      for(int j = 0; j < p; j++)
	sum = (sum + (Fcontainer[j] * coe[p][j]));
      sum = sum * k + U;
      //std::cout << "hello"  << std::endl;
      U = Newton(U, sum, f, J, k, coe[p][p]);
      
      F = f(U);
      ofs << U[0] << " " << U[1] << std::endl;
      Fcontainer.erase(Fcontainer.begin());
      Fcontainer.push_back(F);
      
      //U=calculateU(f,k,p);
      //F=f(U);
      /*for(int i=0;i<p-1;i++)
	Fcontainer[i]=Fcontainer[i+1];
	Fcontainer[p-1]=F;
      Fcontainer.erase(Fcontainer.begin());
      Fcontainer.push_back(f(U));
      U[2]=0.0;U[5]=0.0;
      ofs<<U[0]<<'\t'<<U[1]<<std::endl;*/
    }
    return U;
  };
  static TimeIntegrator* Creator(){
    return new AdamsMoulton();
  }
};


#endif
