#ifndef _SCALARFUNCTION_H_
#define _SCALARFUNCTION_H_
#include<cmath>
#include "vect.h"

class ScalarFunction{
 
 public:
  virtual ~ScalarFunction(){}
  virtual Vec<double,6> operator()(const Vec<double,6>& u) const =0;
  virtual Vec<double,36> getJ(const Vec<double,6> &u) const =0;
};

class Func:public ScalarFunction
{
 private:
  double mu;
 public:
 Func(double _mu):mu(_mu){};
  Vec<double,36> getJ(const Vec<double,6> &u) const override{
    Vec<double,36> temp;
    return temp;
    }
  Vec<double,6> operator()(const Vec<double,6>& u) const override{
    Vec<double,6> temp;
    temp[0]=u[3];temp[1]=u[4];temp[2]=u[5];
    temp[3]=2*u[4]+u[0]-(mu*(u[0]+mu-1))/pow((u[1]*u[1]+u[2]*u[2]+pow((u[0]+mu-1),2)),3.0/2)-(1-mu)*(u[0]+mu)/pow((u[1]*u[1]+u[2]*u[2]+pow(u[0]+mu,2)),3.0/2);
    temp[4]=-2*u[3]+u[1]-mu*u[1]/pow((u[1]*u[1]+u[2]*u[2]+pow((u[0]+mu-1),2)),3.0/2)-(1-mu)*u[1]/pow((u[1]*u[1]+u[2]*u[2]+pow(u[0]+mu,2)),3.0/2);
    temp[5]= -mu*u[2]/pow((u[1]*u[1]+u[2]*u[2]+pow((u[0]+mu-1),2)),3.0/2)-(1-mu)*u[2]/pow((u[1]*u[1]+u[2]*u[2]+pow(u[0]+mu,2)),3.0/2);
    return temp;
  };
};

class Jacobi : public ScalarFunction {
 public:
  Vec<double,6> operator()(const Vec<double,6>& u) const override{
    Vec<double,6> temp;
    return temp;
  };
  Vec<double,36> getJ(const Vec<double,6> &U) const override
  {
    double mu=1.0/81.45;
    Vec<double,36> u;
    double u1=U[0];double u2=U[1];double u3=U[2];
    double u4=U[3];double u5=U[4];double u6=U[5];
    double u22=pow(u2,2);
    double u33=pow(u3,2);
    double u1mu_1=pow(u1+mu-1,2);
    double u1mu=pow(u1+mu,2);

      u[18] = 1;
    u[25] = 1;
    u[32] = 1;
    
    u[3] = 1-(mu*(pow(u22+u33+u1mu_1,3.0/2)-u1mu_1*3*pow(u22+u33+u1mu_1,1.0/2)))/pow(u22+u33+u1mu_1,3)-(1-mu)*(pow(u22+u33+u1mu,3.0/2)-u1mu*3*pow(u22+u33+u1mu,1.0/2))/pow(u22+u33+u1mu,3);
    u[9] =3*mu*(u1+mu-1)*u2/pow(u22+u33+u1mu_1,5.0/2)+3*(1-mu)*(u1+mu)*u2/pow(u22+u33+u1mu,5.0/2);
    u[15] =3*mu*(u1+mu-1)*u3/pow(u22+u33+u1mu_1,5.0/2)+3*(1-mu)*(u1+mu)*u3/pow(u22+u33+u1mu,5.0/2);
    u[27] = 2;
    u[4] = u[9];
    u[10] =1-(mu*(pow(u22+u33+u1mu_1,3.0/2)-3*u22*pow(u22+u33+u1mu_1,1.0/2)))/pow(u22+u33+u1mu_1,3)-(1-mu)*(pow(u22+u33+u1mu,3.0/2)-3*u22*pow(u22+u33+u1mu,1.0/2))/pow(u22+u33+u1mu,3);
    u[16] =3*mu*u2*u3/pow(u22+u33+u1mu_1,5.0/2)+3*(1-mu)*u2*u3/pow(u22+u33+u1mu,5.0/2);
    
    u[22] = -2;
    u[5] = u[15];
    u[11] = u[16];
    u[17] = -mu*(pow(u22+u33+u1mu_1,3.0/2)-3*u33*pow(u22+u33+u1mu_1,1.0/2))/pow(u22+u33+u1mu_1,3)-(1-mu)*(pow(u22+u33+u1mu,3.0/2)-3*u33*pow(u22+u33+u1mu,1.0/2))/pow(u22+u33+u1mu,3);
    
    return u;
  }
  };


#endif
