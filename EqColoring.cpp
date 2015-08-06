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

bool EqColoring::pruningRulePaper(){
  return true;
}

bool EqColoring::pruningRuleFF(){
  return true;
}

bool EqColoring::checkEqColoring(const Graph &g){
  return true;
}
