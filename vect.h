#ifndef _VEC_H_
#define _VEC_H_

#include<iostream>
#include<vector>
#include<cmath>

template<class T,int D>
class Vec{
  T u[D];
 public:
  Vec(){
    for(int i=0;i<D;i++)
      u[i]=0;
  };
  Vec(double *v){
    for(int i=0;i<D;i++)
      u[i]=v[i];
  };
  Vec(const Vec<T,D> & v){
    for(int i=0;i<D;i++)
      u[i]=v[i];
  };
  int size(){
    return D;
  }
  T& operator[](const int i){
    if(i<0 || i>D-1) 
    std::cout<<" operator[] wrong"<<std::endl;
    return u[i];
  };
  const T& operator[](const int i) const{
    if(i<0 || i>D-1) std::cout<<" operator[] wrong"<<std::endl;
    return u[i];
  };
  Vec<T,D>& operator=(const Vec<T,D> & v){
    for(int i=0;i<D;i++)
      u[i]=v[i];
    return *this;
  };
  Vec<T,D> operator+(const Vec<T,D> & v) const{
    Vec<T,D> temp;
    for(int i=0;i<D;i++)
      temp[i]=u[i]+v[i];
    return temp;
  };

  Vec<T,D> operator-(const Vec<T,D> & v) const{
    Vec<T,D> temp;
    for(int i=0;i<D;i++)
      temp[i]=u[i]-v[i];
    return temp;
  };
  
  Vec<T,D> operator*(double k) const{
    Vec<T,D> temp;
    for(int i=0;i<D;i++)
      temp[i]=u[i]*k;
    return temp;
  };

  T maxNorm() const{
    double M=fabs(u[0]);
    for(int i=1;i<D;i++){
      if(fabs(u[i])>M) M=fabs(u[i]);
    }
    return M;
  };

  T errorNorm() const{
    double M=fabs(u[0]);
    if(fabs(u[1])>M) M=fabs(u[1]);
    return M;
  };

  friend Vec<T,D> operator-(const Vec<T,D> &a){
    Vec<T,D> temp;
    for(int i=0;i<D;i++)
      temp[i]=-a[i];
    return temp;
  };
  
};

#endif
