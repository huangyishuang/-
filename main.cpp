#include "vect.h"
#include "ScalarFunction.h"
#include "AdamsMoulton.h"
#include "AdamsBashforth.h"
#include "RungeKutta.h"
#include "BDF.h"
#include "TimeIntegratorFactory.h"
#include<iostream>
#include<fstream>
#include<string>
#include<time.h>
#include<cmath>
#include <iomanip>
using namespace std;

int main(int argc,char *argv[]){
  
  TimeIntegratorFactory Factory=TimeIntegratorFactory::Instance();
  Factory.Register("AdamsBashforth",&AdamsBashforth::Creator);
  Factory.Register("AdamsMoulton",&AdamsMoulton::Creator);
  Factory.Register("BDF",&BDF::Creator);
  Factory.Register("RungeKutta",&RungeKutta::Creator);

  ifstream is("input.data");
  string str;int pp;double t;

  double v[6]={0.994,0,0,0,-2.0015851063790825224,0};
  double dv;
  for(int i=0;i<6;i++){
    is>>dv;
    v[i]=dv;
  }
    double T;
    is>>T;
    /*double v[6]={0.87978,0,0,0,-0.3797,0};
      double T=19.14045706162071;*/
    
  double mu=1.0/81.45;
  Vec<double,6> vec(v);
  Func f(mu);
  
  ofstream ofs("result.data");
  ofstream ofs1("report.data");
  Vec<double,6> approU;
  double error;
  while(is>>str>>pp){
    ofs1<<str<<" "<<pp<<endl;
    
    //ofstream ofs("./data/result_"+str+to_string(pp)+".data");
    ofs1<<std::setw(8)<<"step length n "<<std::setw(13)<<"Solution error "<<std::setw(13)<<" convergence rate "<<std::setw(10)<<" time "<<std::endl;
    double e1=1.0,e2;double rate;
    //for(int N=5000;N<1280000+1;N=N*2){
    for(int N=50000;N<50000+1;N=N+1){
    clock_t start,end;
    start=clock();
    TimeIntegrator* obj=Factory.CreateObject(str);
    obj->init(vec,f,ofs,pp);
    approU=obj->solve(ofs,f,N,T,pp);
    e2=(approU-vec).errorNorm();
    rate=log(e1/e2)/log(2);
    end=clock();
    t=double(end-start)/CLOCKS_PER_SEC;
    ofs1<<std::setw(8)<< N <<std::setw(15)<<e2<<std::setw(15)<<rate<<std::setw(17)<<t<<std::endl;
    e1=e2;
    }
  }
  return 0;
}
