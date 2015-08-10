#include "Coloring.hpp"

Coloring::Coloring(){
  Heuristic::constructRandomGraph(g, 10, 0.5);

  if(!initVar()){
    exit(0);
  }

  printAdjMatrix();

  std::vector<Vertex> clq;
  findMaxClique(clq, true, true);
  printClique(clq);
}

Coloring::Coloring(const Parameters &parm){
  if(!initVar()){
    exit(0);
  }
  
  this->parm = parm;
}

Coloring::Coloring(const Graph &g){
  if(!initVar()){
    exit(0);
  }
  
  this->g = g;
}

Coloring::Coloring(const Graph &g, const Parameters &parm){
  if(!initVar()){
    exit(0);
  }

  this->g = g;
  this->parm = parm;
}

Coloring::~Coloring(){

}

bool Coloring::compareDegree(Vertex v, Vertex w){
  return (in_degree(v, g) < in_degree(w, g));
}

bool Coloring::putNodeWithParm(Vertex v, std::vector<Vertex> &tmp, bool uncolored, bool inNoOtherClique){
  if(uncolored == true && inNoOtherClique == true){
    if(pm.c[v] == 0 && pm.cl[v] == 0){
      tmp.push_back(v);
    }else{
      return false;
    }
  }else if(uncolored == true && inNoOtherClique == false){
    if(pm.c[v] == 0){
      tmp.push_back(v);
    }else{
      return false;
    }
  }else if(uncolored == false && inNoOtherClique == true){
    if(pm.cl[v] == 0){
      tmp.push_back(v);
    }else{
      return false;
    }
  }else{
    tmp.push_back(v);
  }

  return true;
}

bool Coloring::findMaxClique(std::vector<Vertex> &clq, Vertex v, bool uncolored, bool inNoOtherClique){
  vertexIter vIt1, vIt2;
  adjaIter aIt1, aIt2;
  std::vector<Vertex> neighbours, tmp;
  Vertex w;
  
  if(!putNodeWithParm(v, clq, uncolored, inNoOtherClique)){
    return true;
  }

  for(tie(aIt1,aIt2) = adjacent_vertices(v,g); aIt1 != aIt2; aIt1++){
    putNodeWithParm(*aIt1, neighbours, uncolored, inNoOtherClique);
  }

  while(neighbours.size() != 0){
    std::sort(neighbours.begin(), neighbours.end(), boost::bind(&Coloring::compareDegree, this, _1, _2));

    w = neighbours.back();
    neighbours.pop_back();
    clq.push_back(w);

    for(tie(aIt1,aIt2) = adjacent_vertices(w,g); aIt1 != aIt2; aIt1++){
      if(std::find(neighbours.begin(), neighbours.end(), *aIt1) != neighbours.end()){
        tmp.push_back(*aIt1);  
      } 
    }

    neighbours = tmp;
    tmp.clear();
  }

  return true;
}

void Coloring::printClique(const std::vector<Vertex> &clq) const{
  std::cout << "Clique: ";

  for(unsigned int i = 0; i < clq.size(); i++){
    std::cout << pm.i[clq[i]] << " ";
  }

  std::cout << std::endl;
}

bool Coloring::findMaxClique(std::vector<Vertex> &clq, bool uncolored, bool inNoOtherClique){
  vertexIter vIt1, vIt2;
  adjaIter aIt1, aIt2;

  std::vector<Vertex> tmpClique, maxClique;

  maxClique.clear();

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    tmpClique.clear();

    findMaxClique(tmpClique, *vIt1, uncolored, inNoOtherClique);

    if(tmpClique.size() > maxClique.size()){
      maxClique = tmpClique;
    }
  }

  if(!checkClique(maxClique)){
    std::cout << "not a clique!" << std::endl;
  }

  clq = maxClique;

  return true;
}

bool Coloring::checkClique(const std::vector<Vertex> &clq) const{
  for(unsigned int i = 0; i < clq.size(); i++){
    for(unsigned int j = i+1; j < clq.size(); j++){
      if(pm.n[clq[j]][pm.i[clq[i]]] == 0){
        return false;
      }
    }
  }

  return true;
}

bool Coloring::initPropMap(){
  pm.r = get(rank_t(),g);
  pm.c = get(vertex_color, g);
  pm.cl = get(clique_t(), g);
  pm.n = get(nachbarn_t(), g);
  pm.fbc = get(fb_colors_t(), g);
  pm.i = get(vertex_index, g);
  pm.g = get(grad_t(), g);

  return true;
}

bool Coloring::initNeighbours(){
  vertexIter vIt1, vIt2;
  adjaIter aIt1, aIt2;

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    for(tie(aIt1,aIt2) = adjacent_vertices(*vIt1,g); aIt1 != aIt2; aIt1++){
      pm.n[*vIt1][pm.i[*aIt1]] = 1;
    }
  }

  return true;
}

bool Coloring::setBounds(int LB, int UB){
  b.LB = LB;
  b.UB = UB;

  return true;
}

bool Coloring::setCurr(int c, int r, Vertex node, int uncoloredVertices){
  curr.color = c;
  curr.rank = r;
  curr.node = node;
  curr.uncoloredVertices = uncoloredVertices;

  return true;
}

bool Coloring::setBacktracking(bool status, Vertex node){
  bt.status = status;
  bt.toNode = node;

  return true;
}

bool Coloring::initVar(){
  if(!initPropMap()){
    std::cout << "error while loading property map" << std::endl;
    
    return false;
  }

  if(!setBounds(0, 10)){
    std::cout << "error while setting init bounds" << std::endl;
    
    return false;
  }

  if(!initResize()){
    std::cout << "problem by resize init vector" << std::endl;
    
    return false;
  }

  if(!setBacktracking(false, 0)){
    std::cout << "error while setting init backtracking" << std::endl;
    
    return false;
  }

  if(!setCurr(0, 0, 0, 0)){
    std::cout << "error while setting init current information" << std::endl;

    return false;
  }

  if(!initNeighbours()){
    std::cout << "init neighbours vector failed!" << std::endl;
    
    return false;
  }

  return true;
}

bool Coloring::initResize(){
  vertexIter vIt1, vIt2;

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    pm.n[*vIt1].resize(b.UB, 0);
    pm.fbc[*vIt1].resize(b.UB, 0);
  }

  cc.n.resize(b.UB, 0);

  return true;
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


