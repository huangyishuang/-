#ifndef _TIMEINTEGRATORFACTORY_H_
#define _TIMEINTEGRATORFACTORY_H_

#include <map>
#include <string>
#include "TimeIntegrator.h"

class TimeIntegratorFactory{
  typedef TimeIntegrator* (*ProductCreator)();
  typedef std::map<std::string, ProductCreator> Map;
  Map m;
public:
  ~TimeIntegratorFactory(){
    m.clear();
  }
  static TimeIntegratorFactory& Instance(){
    static TimeIntegratorFactory instance;
    return instance;
  }
  bool Register(const std::string& id, ProductCreator creator){
    return m.insert(Map::value_type(id, creator)).second;
  }
  bool Unregister(const std::string& id){
    return m.erase(id) == 1;
  }
  TimeIntegrator* CreateObject(const std::string& id){
    typename Map::const_iterator i = m.find(id);
    if(i != m.end()){
      return (i->second)();
    }
    std::cout << "No product find for " << id << std::endl;
    return NULL;
  }
};

#endif
