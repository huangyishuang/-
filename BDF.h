#ifndef _BDF_H_
#define _BDF_H_
#include<fstream>
#include "vect.h"
#include <vector>
#include "ScalarFunction.h"
#include "TimeIntegrator.h"
#include "lapacke.h"
#include<cmath>

class BDF:public LMM{
 protected:
  std::vector<Vec<double,6>> Ucontainer;
  std::vector<std::vector<double> > coe={{0},{-1.0},{1.0/3,-4.0/3},{-2.0/11,9.0/11,-18.0/11},{3.0/25,-16.0/25,36.0/25,-48.0/25}};
  std::vector<double> beta ={0,1.0,2.0/3,6.0/11,12.0/25};
 public:
  BDF()=default;
  /*BDF(const Vec<double,6>& _u,const ScalarFunction & f,std::ofstream & ofs,const int p)
   {
     U=_u;
      F=f(U);
      Ucontainer.resize(p);
      Ucontainer[0]=U;
      ofs<<U[0]<<'\t'<<U[1]<<std::endl;
      };*/

  void init(const Vec<double,6>& _u,const ScalarFunction & f,std::ofstream & ofs,int p=1)
   {
     U=_u;
      F=f(U);
      Ucontainer.erase(Ucontainer.begin(),Ucontainer.end());
      //Ucontainer.resize(p);
      //Ucontainer[0]=U;
      Ucontainer.push_back(U);
      ofs<<U[0]<<'\t'<<U[1]<<std::endl;
   };
  
  void onestep(std::ofstream & ofs,const ScalarFunction & f,const double k,const int p){
    for(int i=1;i<p;i++){
      U=U+F*k;
      F=f(U);
      Ucontainer.push_back(U);
      //Ucontainer[i]=U;
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
      //U=calculateU(f,k,p);
      Vec<double,6> sum;
      for(int j = 0; j < p; j++)
	sum = sum - (Ucontainer[j] * coe[p][j]);
      U = Newton(U, sum, f, J, k, beta[p]);

      
      Ucontainer.erase(Ucontainer.begin());
      Ucontainer.push_back(U);
      //for(int i=0;i<p-1;i++)
      //	Ucontainer[i]=Ucontainer[i+1];
      //Ucontainer[p-1]=U;
      U[2]=0.0;U[5]=0.0;
      ofs<<U[0]<<'\t'<<U[1]<<std::endl;
    }
    return U;
  };
  static TimeIntegrator* Creator(){
    return new BDF();
  }
  
};


#endif
