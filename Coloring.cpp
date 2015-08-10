#include "Coloring.hpp"

Coloring::Coloring(){
  parm.variant = 'R';

  if(parm.variant == 'R'){
    Heuristic::constructRandomGraph(g, 10, 0.2);
  }else if(parm.variant == 'G'){
    Input::readGraph("res/queen7_7.col", g);
  }

  setParm(num_vertices(g));

  if(!initVar()){
    exit(0);
  }

  printAll();
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

bool Coloring::initCliques(){
  findMaxClique(startClique, true, true);
  putInClique(startClique, 1);
  cl.nodesInClique += startClique.size();
  cl.nCliques = 1;
  colorClique(startClique, 1);
  findIndepCliques(indClq, true, true, 1);

  return true;
}

void Coloring::printAll() const{
  printAdjMatrix();
  printClique(startClique);
  printIndepCliques(indClq);
  printCliqueInfo();
  printVertexInfo();
}

void Coloring::printCliqueInfo() const{
  std::cout << "nodes in clique: " << cl.nodesInClique << std::endl;
  std::cout << "number of cliques: " << cl.nCliques << std::endl;
  std::cout << "status new clique: " << cl.newClique << std::endl;
}

bool Coloring::checkIndepClique(std::vector<std::vector<Vertex> > &indClq){
  return true;
}

bool Coloring::colorClique(std::vector<Vertex> &clq, int startColor){
  int color = startColor;

  for(unsigned int i = 0; i < clq.size(); i++){
    pm.c[clq[i]] = color;
    pm.r[clq[i]] = curr.rank + 1;
    color++;
    curr.uncoloredVertices--;
    curr.rank++;
  }

  return true;
}

bool Coloring::setParm(long n, double p, long npr, long tl, long th, std::string res, long nrg, char variant){
  parm.n = n;
  parm.ressource = res;
  parm.nPruningRule = npr;
  parm.threshold = th;
  parm.timeLimit = tl;
  parm.p = p;
  parm.nRandomGraphs = nrg;
  parm.variant = variant;

  return true;
}

bool Coloring::compareDegree(Vertex v, Vertex w){
  return (in_degree(v, g) < in_degree(w, g));
}

bool Coloring::putInClique(std::vector<Vertex> &clq, int toClique){
  for(unsigned int i = 0; i < clq.size(); i++){
    pm.cl[clq[i]] = toClique;
  }

  return true;
}

void Coloring::printIndepCliques(const std::vector<std::vector<Vertex> > &indClq) const{
  for(unsigned int i = 0; i < indClq.size(); i++){
    printClique(indClq[i]);
  }
}

bool Coloring::findIndepCliques(std::vector<std::vector<Vertex> > &indClq, bool uncolored, bool inNoOtherClique, int iFirstClique){
  adjaIter aIt1, aIt2;
  vertexIter vIt1, vIt2;
  int toVisit = curr.uncoloredVertices, toClique = iFirstClique;
  std::vector<Vertex> tmp;

  while(toVisit > 0){
    tmp.clear();
    
    findMaxClique(tmp, uncolored, inNoOtherClique);
    indClq.push_back(tmp);

    toClique++;
    putInClique(tmp, toClique);

    toVisit -= tmp.size();
    cl.nodesInClique += tmp.size();

    //mark neighbours as visited -> Clique: -1
    for(unsigned int i = 0; i < tmp.size(); i++){
      for(tie(aIt1,aIt2) = adjacent_vertices(tmp[i],g); aIt1 != aIt2; aIt1++){
        if(pm.cl[*aIt1] == 0 && pm.c[*aIt1] == 0){
          pm.cl[*aIt1] = -1;
          toVisit--;
        }
      }
    }
  }

  cl.nCliques = toClique;

  //remove visited Nodes
  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    if(pm.cl[*vIt1] == -1){
      pm.cl[*vIt1] = 0;
    }
  }

  return true;
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

bool Coloring::setClique(long nodesInClique, long nCliques, bool newClique){
  cl.nodesInClique = nodesInClique;
  cl.nCliques = nCliques;
  cl.newClique = newClique;

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

  if(!setBounds(0, 49)){
    std::cout << "error while setting init bounds" << std::endl;
    
    return false;
  }

  if(!initResize()){
    std::cout << "problem by resize init vector" << std::endl;
    
    return false;
  }

  if(!setClique(0, 0, false)){
    std::cout << "setting init clique" << std::endl;
    
    return false;
  }

  if(!setBacktracking(false, 0)){
    std::cout << "error while setting init backtracking" << std::endl;
    
    return false;
  }

  if(!setCurr(0, 0, 0, parm.n)){
    std::cout << "error while setting init current information" << std::endl;

    return false;
  }

  if(!initNeighbours()){
    std::cout << "init neighbours vector failed!" << std::endl;
    
    return false;
  }

  if(!initCliques()){
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


