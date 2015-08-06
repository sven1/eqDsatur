#include "Coloring.hpp"

Coloring::Coloring(){
  std::cout << "Test A" << std::endl;
}

Coloring::Coloring(const Parameters &parm){
  this->parm = parm;
}

Coloring::Coloring(const Graph &g){
  this->g = g;
}

Coloring::Coloring(const Graph &g, const Parameters &parm){
  this->g = g;
  this->parm = parm;
}

Coloring::~Coloring(){

}

long Coloring::naiveUB(const Graph &g){
  return 0;
}

long Coloring::eqlLB(const Graph &g){
  return 0;
}

bool Coloring::node(Graph &g){
  return true;
}

void Coloring::printGraph(const Graph &g){

}

bool Coloring::checkColoring(const Graph &g){
  return true;
}
