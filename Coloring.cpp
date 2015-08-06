#include "Coloring.hpp"

Coloring::Coloring(){
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

void Coloring::printVertexInfo() const{
  vertexIter vIt1, vIt2;

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    printVertexInfo(*vIt1);

    std::cout << std::endl;
  }
}

void Coloring::printVertexInfo(const Vertex &v) const{
    std::cout << pm.i[v] << ": ";
    std::cout << pm.r[v] << " ";
    std::cout << pm.c[v] << " ";
    std::cout << pm.cl[v] << " ";
    std::cout << pm.g[v];
}

void Coloring::printAdjMatrix() const{
  vertexIter vIt1, vIt2;

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    printNeighbours(*vIt1);

    std::cout << std::endl;
  }
}

void Coloring::printNeighbours(const Vertex &v) const{
  std::cout << pm.i[v] << ":";

  for(unsigned int i = 0; i < pm.n[v].size(); i++){
    std::cout << " " << pm.n[v][i];
  }
}

void Coloring::printFBC(const Vertex &v) const{
  std::cout << pm.i[v] << ":";

  for(unsigned int i = 0; i < pm.fbc[v].size(); i++){
    std::cout << " " << pm.fbc[v][i];
  }
}

void Coloring::printFBC() const{
  vertexIter vIt1, vIt2;

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    printFBC(*vIt1);

    std::cout << std::endl;
  }
}

bool Coloring::checkColoring() const{
  vertexIter vIt1, vIt2;

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    if(!checkColoring(*vIt1)){
      return false;
    }
  }

  return true;
}

bool Coloring::checkColoring(const Vertex &v) const{
  outEdgeIter oIt1, oIt2;

  for(tie(oIt1,oIt2) = out_edges(v,g); oIt1 != oIt2; oIt1++){
    if(pm.c[source(*oIt1,g)] == pm.c[target(*oIt1,g)]){
      return false;   
    }
  }

  return true;
}
