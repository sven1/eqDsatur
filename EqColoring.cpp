#include "EqColoring.hpp"

EqColoring::EqColoring() : Coloring(){
  std::cout << "Test B" << std::endl;
}

EqColoring::EqColoring(const Parameters &parm) : Coloring(parm){

}

EqColoring::EqColoring(const Graph &g) : Coloring(g){

}

EqColoring::EqColoring(const Graph &g, const Parameters &parm) : Coloring(g, parm){

}

EqColoring::~EqColoring(){

}

bool EqColoring::node(Graph &g){
  return true;
}

bool EqColoring::pruningRulePaper(){
  return true;
}

bool EqColoring::pruningRuleFF(){
  return true;
}

bool EqColoring::checkEquitability() const{
  for(unsigned int i = 0; i < cls.n.size(); i++){
    for(unsigned int j = i+1; j < cls.n.size(); j++){
      if(std::abs(cls.n[i] - cls.n[j]) > 1){
        return false; 
      }  
    }
  }

  return true;
}

bool EqColoring::checkEqColoring() const{
  if(!Coloring::checkColoring()){
    return false;
  }

  if(!checkEquitability()){
    return false;
  }

  return true;
}
