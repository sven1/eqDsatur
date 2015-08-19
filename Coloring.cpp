#include "Coloring.hpp"

Coloring::Coloring(){
  if(!initVar()){
    exit(0);
  }

  printAll();
}

Coloring::Coloring(const Parameters &parm){
  this->parm = parm;

  if(this->parm.variant == 'R'){
    Heuristic::constructRandomGraph(g, this->parm.n, this->parm.p);
  }else if(this->parm.variant == 'N'){
    Input::readGraph(this->parm.ressource, g);

    this->parm.n = num_vertices(g);
  }

  if(!initVar()){
    exit(0);
  }
  
  printAll();
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

int Coloring::calcLB(){
  return startClique.size();
}

int Coloring::calcUB(){
  Graph tmpG;
  Current tmp;

  copy_graph(g, tmpG);
  tmp = curr;
  greedyColoring(g);
  printAll();
  g = tmpG;
  curr = tmp;

  return parm.n;
}

bool Coloring::initCliques(){
  findMaxClique(startClique, true, true);
  putInClique(startClique);
  colorClique(startClique, 1);
  findIndepCliques(indClq, true, true);

  return true;
}

bool Coloring::setTandM(int T, int M){
  curr.M = M;
  curr.T = T;

  return true;
}

void Coloring::printGraphHeaders() const{
  std::cout << "n: " << parm.n << std::endl;
  std::cout << "how often do FF pruning rule (mod): " << parm.nPruningRule << std::endl;
  std::cout << "time limit: " << parm.timeLimit << std::endl;
  std::cout << "threshold: " << parm.threshold << std::endl;

  if(parm.variant == 'N'){
    std::cout << "ressource file: " << parm.ressource << std::endl;
  }else if(parm.variant == 'R'){
    std::cout << "p: " << parm.p << std::endl;
    std::cout << "number of random graphs: " << parm.nRandomGraphs << std::endl;
  }
}

void Coloring::printCurrent() const{
  std::cout << "amount of used colors: " << curr.nColors << std::endl;
  std::cout << "current node: " << curr.node << std::endl;
  std::cout << "current color: " << curr.color << std::endl;
  std::cout << "current rank: " << curr.rank << std::endl;
  std::cout << "M (cardinality of greatest color class): " << curr.M << std::endl;
  std::cout << "T (amount of greatest color classes): " << curr.T << std::endl;
  std::cout << "current uncoloredVertices: " << curr.uncoloredVertices << std::endl;
}

void Coloring::printColorClass(int i) const{
  std::cout << "color class " << i << ": " << cc.n[i-1] << std::endl;
}

void Coloring::printColorClasses() const{
  for(int i = 1; i <= curr.nColors; i++){
    printColorClass(i);
  }
}

void Coloring::printAll() const{
  printGraphHeaders();
  printCurrent();
  printBounds();
  printAdjMatrix();
  printClique(startClique);
  printIndepCliques(indClq);
  printCliqueInfo();
  printVertexInfo();
  printFBC();
  printColorClasses();
}

void Coloring::printCliqueInfo() const{
  std::cout << "nodes in clique: " << cl.nodesInClique << std::endl;
  std::cout << "number of cliques: " << cl.nCliques << std::endl;
  std::cout << "status new clique: " << cl.newClique << std::endl;
}

bool Coloring::checkIndepClique(std::vector<std::vector<Vertex> > &indClq){
  return true;
}

// passt rank und uncoloredVertices direkt an
bool Coloring::colorClique(std::vector<Vertex> &clq, int startColor){
  int color = startColor;

  for(unsigned int i = 0; i < clq.size(); i++){
    colorVertex(clq[i], color);

    color++;
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

// passt Knoten in Clique an und anzahl cliquen
bool Coloring::putInClique(std::vector<Vertex> &clq){
  cl.nCliques++;

  for(unsigned int i = 0; i < clq.size(); i++){
    pm.cl[clq[i]] = cl.nCliques;
  }

  cl.nodesInClique += clq.size();

  return true;
}

void Coloring::printIndepCliques(const std::vector<std::vector<Vertex> > &indClq) const{
  for(unsigned int i = 0; i < indClq.size(); i++){
    printClique(indClq[i]);
  }
}

bool Coloring::findIndepCliques(std::vector<std::vector<Vertex> > &indClq, bool uncolored, bool inNoOtherClique){
  adjaIter aIt1, aIt2;
  vertexIter vIt1, vIt2;
  int toVisit = curr.uncoloredVertices;
  std::vector<Vertex> tmp;

  indClq.clear();

  while(toVisit > 0){
    tmp.clear();
    
    findMaxClique(tmp, uncolored, inNoOtherClique);
    indClq.push_back(tmp);

    putInClique(tmp);

    toVisit -= tmp.size();

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

    return false;
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

void Coloring::printBounds() const{
  std::cout << "UB: " << b.UB << std::endl;
  std::cout << "LB: " << b.LB << std::endl;
}

bool Coloring::setClique(long nodesInClique, long nCliques, bool newClique){
  cl.nodesInClique = nodesInClique;
  cl.nCliques = nCliques;
  cl.newClique = newClique;

  return true;
}

bool Coloring::setCurr(int c, int r, Vertex node, int uncoloredVertices, int T, int M, long nColors){
  curr.color = c;
  curr.rank = r;
  curr.node = node;
  curr.uncoloredVertices = uncoloredVertices;
  curr.T = T;
  curr.M = M;
  curr.nColors = nColors;

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

  if(!setBounds(0, parm.n)){
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

  if(!setCurr(0, 0, 0, parm.n, 0, 0, 0)){
    std::cout << "error while setting init current information" << std::endl;

    return false;
  }

  if(!initNeighbours()){
    std::cout << "init neighbours vector failed!" << std::endl;
    
    return false;
  }

  if(!initCliques()){
    std::cout << "init cliques failed!" << std::endl;
    
    return false;
  }

  if(!setBounds(calcLB(), b.UB)){
    std::cout << "error while setting calc bounds" << std::endl;
    
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

bool Coloring::passVSS(){
  return true;
}

bool Coloring::node(){
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

bool Coloring::greedyColoring(Graph &g){
  vertexIter vIt1, vIt2;
  int sColor;
  long color = curr.nColors;

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    if(pm.c[*vIt1] == 0){
      sColor = smallestPosColor(*vIt1);

      if(sColor == -1){
        color++;
        colorVertex(*vIt1, color);
      }else{
        colorVertex(*vIt1, sColor);
      }
    }
  }

  if(!checkColoring()){
    std::cout << "No coloring!" << std::endl;
  }

  return true;
}

int Coloring::smallestPosColor(Vertex v) const{
  for(unsigned int i = 0; i < pm.fbc[v].size() && i < curr.nColors; i++){
    if(pm.fbc[v][i] == 0){
      return i + 1;
    }
  }

  return -1;
}

bool Coloring::colorVertex(Vertex v, int color){
  curr.rank++;
  curr.uncoloredVertices--;
  
  pm.c[v] = color;
  pm.r[v] = curr.rank;

  if(color > curr.nColors){
    if(std::abs(color - curr.nColors) > 1){
      std::cout << "differenz der farben zu gross" << std::endl;
    }

    curr.nColors++;
  }

  addFBC(v, color);
  incColorClass(color);

  return true;
}

bool Coloring::updateTandM(int lastColor){
  if(cc.n[lastColor - 1] > curr.M){
    curr.M = cc.n[lastColor - 1];
    curr.T = 1;
  }else if(cc.n[lastColor - 1] == curr.M){
    curr.T++;
  }

  return true;
}

bool Coloring::incColorClass(int color){
  if((int) cc.n.size() > color){
    cc.n[color - 1]++;

    updateTandM(color);

    return true;
  }else{
    std::cout << "not enough color classes" << std::endl;

    return false;
  }
}

bool Coloring::addFBC(Vertex v, int color){
  adjaIter aIt1, aIt2;

  for(tie(aIt1,aIt2) = adjacent_vertices(v,g); aIt1 != aIt2; aIt1++){
    pm.fbc[*aIt1][color - 1]++;
  }

  return true;
}
