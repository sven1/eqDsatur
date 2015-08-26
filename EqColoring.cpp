#include "EqColoring.hpp"

EqColoring::EqColoring() : Coloring(){
}

EqColoring::EqColoring(const Parameters &parm) : Coloring(parm){
  setBounds(b.LB, calcUB());

  initPrevGraphsFF();

  if(parm.eqDsatur == 'N'){
    dsatur();
  }else if(parm.eqDsatur == 'C'){
    dsaturClique();
  }
  
  printAll();

  printCounts();
}

bool EqColoring::initPrevGraphsFF(){
  prevGraphsFF.resize(b.UB);
  pmPrevGraphsFF.resize(b.UB);

  for(unsigned int i = 0; i < pmPrevGraphsFF.size(); i++){
    pmPrevGraphsFF[i].c = get(edge_capacity, prevGraphsFF[i].g);
    pmPrevGraphsFF[i].re = get(edge_reverse, prevGraphsFF[i].g);
    pmPrevGraphsFF[i].rc = get(edge_residual_capacity, prevGraphsFF[i].g);
    pmPrevGraphsFF[i].rf = get(ref_vertex_t(), prevGraphsFF[i].g);
  }

  return true;
}

EqColoring::EqColoring(const Graph &g) : Coloring(g){

}

EqColoring::EqColoring(const Graph &g, const Parameters &parm) : Coloring(g, parm){

}

EqColoring::~EqColoring(){

}

bool EqColoring::node(){
  if(curr.uncoloredVertices == 0){
    b.UB = curr.nColors;

    return true;
  }

  Vertex v = passVSS();


  for(int i = 1; i <= std::min(b.UB - 1, curr.nColors + 1); i++){
    if(pm.fbc[v][i - 1] == 0){
      if(parm.n >= (curr.M - 1) * std::max(curr.nColors, b.LB) + curr.T){
        c.visitedNodes++;

        colorVertex(v, i);

        node(); 
          
        if(bt.status){
          if(bt.toRank == curr.rank){
            bt.status = false;
          }else if(bt.toRank > curr.rank){
            uncolorVertex(v);
            return true;
          }
        }

        uncolorVertex(v);
      }
    }
  }

  checkForBacktracking(v);

  return false;
}

bool EqColoring::pruneFF(){
  for(unsigned int i = curr.nColors; i < b.UB; i++){
    c.nFF++;

    if(!pruneFF(i)){
      return false;
    }
  }

  return true;
}

bool EqColoring::initA1(std::vector<VertexFord> &vert){
  EdgeFord eF1, eF2;

  for(int i = 1; i <= curr.uncoloredVertices; i++){
    eF1 = add_edge(vert[0], vert[i], gf).first;
    eF2 = add_edge(vert[i], vert[0], gf).first;

    pmf.re[eF1] = eF2;
    pmf.re[eF2] = eF1;

    pmf.c[eF1] = 1;
    pmf.c[eF2] = 0;
  }

  return true;
}

bool EqColoring::initA2andA3(std::vector<VertexFord> &vert, int color){
  EdgeFord eF1, eF2;
  VertexFordIter vIt1, vIt2;

  int aUVPos = 1;
  int tmpIndex, colorPos;

  for(unsigned int i = 0; i < indClq.size(); i++){
    for(unsigned int j = 0; j < indClq[i].size(); j++){
      pmf.rf[vert[aUVPos]] = indClq[i][j];

      for(int k = 0; k < color; k++){
        if(pm.fbc[indClq[i][j]][k] == 0){
          tmpIndex = curr.uncoloredVertices + (k + 1) + i * color; 

          eF1 = add_edge(vert[aUVPos], vert[tmpIndex], gf).first;
          eF2 = add_edge(vert[tmpIndex], vert[aUVPos], gf).first;

          pmf.re[eF1] = eF2;
          pmf.re[eF2] = eF1;

          pmf.c[eF1] = 1;
          pmf.c[eF2] = 0;

          //zur Farbe
          colorPos = curr.uncoloredVertices + (k + 1) + cl.nCliques * color; 

          eF1 = add_edge(vert[tmpIndex], vert[colorPos], gf).first;
          eF2 = add_edge(vert[colorPos], vert[tmpIndex], gf).first;

          pmf.re[eF1] = eF2;
          pmf.re[eF2] = eF1;

          pmf.c[eF1] = 1;
          pmf.c[eF2] = 0;
        }
      }

      aUVPos++;
    }
  }

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    if(pm.c[*vIt1] == 0 && pm.cl[*vIt1] == 0){
      pmf.rf[vert[aUVPos]] = *vIt1;

      for(int k = 0; k < color; k++){
        if(pm.fbc[*vIt1][k] == 0){
          tmpIndex = curr.uncoloredVertices + (k + 1) + (cl.nCliques - 1) * color; 

          eF1 = add_edge(vert[aUVPos], vert[tmpIndex], gf).first;
          eF2 = add_edge(vert[tmpIndex], vert[aUVPos], gf).first;

          pmf.re[eF1] = eF2;
          pmf.re[eF2] = eF1;

          pmf.c[eF1] = 1;
          pmf.c[eF2] = 0;

          colorPos = curr.uncoloredVertices + (k + 1) + cl.nCliques * color; 

          eF1 = add_edge(vert[tmpIndex], vert[colorPos], gf).first;
          eF2 = add_edge(vert[colorPos], vert[tmpIndex], gf).first;

          pmf.re[eF1] = eF2;
          pmf.re[eF2] = eF1;

          pmf.c[eF1] = 1;
          pmf.c[eF2] = 0;
        }
      }

      aUVPos++;
    }
  }

  if(aUVPos != curr.uncoloredVertices + 1){
    std::cout << "aUVPos = " << aUVPos << " and uncoloredVertices = " << curr.uncoloredVertices << std::endl;
    std::cout << "error while constructing network" << std::endl;
    std::cout << "not using all uncolored Vertices" << std::endl;

    return false;
  }

  return true;
}

int EqColoring::initA4(std::vector<VertexFord> &vert, int color, std::pair<int, int> counts){
  int sumLB = 0, rU, rL, colorPos, n = counts.first, nNew = counts.second;
  EdgeFord eF1, eF2;
  

  for(int i = 1; i <= color; i++){
    rU = parm.n / color + 1;
    rL = parm.n / color;

    rU -= cc.n[i-1];
    rL -= cc.n[i-1];

    if(rL < 0){
      rL = 0;
    }

    sumLB += rL;

    colorPos = curr.uncoloredVertices + i + cl.nCliques * color; 

    eF1 = add_edge(vert[colorPos], vert[n - 1], gf).first;
    eF2 = add_edge(vert[n - 1], vert[colorPos], gf).first;

    pmf.re[eF1] = eF2;
    pmf.re[eF2] = eF1;

    pmf.c[eF1] = rU - rL;
    pmf.c[eF2] = 0;

    eF1 = add_edge(vert[colorPos], vert[nNew - 1], gf).first;
    eF2 = add_edge(vert[nNew - 1], vert[colorPos], gf).first;

    pmf.re[eF1] = eF2;
    pmf.re[eF2] = eF1;

    pmf.c[eF1] = rL;
    pmf.c[eF2] = 0;
  }

  return sumLB;
}

bool EqColoring::initRespectLB(std::vector<VertexFord> &vert, std::pair<int, int> counts, int sumLB){
  EdgeFord eF1, eF2;
  int n = counts.first, nNew = counts.second;

  eF1 = add_edge(vert[nNew - 2], vert[n - 1], gf).first;
  eF2 = add_edge(vert[n - 1], vert[nNew - 2], gf).first;

  pmf.re[eF1] = eF2;
  pmf.re[eF2] = eF1;

  pmf.c[eF1] = sumLB;
  pmf.c[eF2] = 0;

  eF1 = add_edge(vert[n - 1], vert[0], gf).first;
  eF2 = add_edge(vert[0], vert[n - 1], gf).first;

  pmf.re[eF1] = eF2;
  pmf.re[eF2] = eF1;

  pmf.c[eF1] = INT_MAX;
  pmf.c[eF2] = 0;

  return true;
}

bool EqColoring::removeRespectLB(std::vector<VertexFord> &vert, int color, std::pair<int, int> counts, int sumLB){
  int rL, rU, colorPos, n = counts.first, nNew = counts.second;
  EdgeFord eF1, eF2;

  for(int i = 1; i <= color; i++){
    rU = parm.n / color + 1;
    rL = parm.n / color;

    rU -= cc.n[i-1];
    rL -= cc.n[i-1];

    if(rL < 0){
      rL = 0;
    }

    colorPos = curr.uncoloredVertices + i + cl.nCliques * color; 

    eF1 = edge(vert[colorPos], vert[n - 1], gf).first;
    eF2 = edge(vert[n - 1], vert[colorPos], gf).first;

    pmf.c[eF1] = rU;
    pmf.c[eF2] = 0;

    pmf.rc[eF1] = rU - rL;
    pmf.rc[eF1] = 0;

    eF1 = edge(vert[colorPos], vert[nNew - 1], gf).first;
    eF2 = edge(vert[nNew - 1], vert[colorPos], gf).first;

    pmf.c[eF1] = 0;
    pmf.c[eF2] = 0;

    pmf.rc[eF1] = 0;
    pmf.rc[eF1] = 0;
  }
    
  eF1 = edge(vert[nNew - 2], vert[n - 1], gf).first;
  eF2 = edge(vert[n - 1], vert[nNew - 2], gf).first;

  pmf.c[eF1] = sumLB;
  pmf.c[eF2] = 0;

  pmf.rc[eF1] = 0;
  pmf.rc[eF1] = 0;
    
  eF1 = edge(vert[n - 1], vert[0], gf).first;
  eF2 = edge(vert[0], vert[n - 1], gf).first;

  pmf.c[eF1] = 0;
  pmf.c[eF2] = 0;

  pmf.rc[eF1] = 0;
  pmf.rc[eF1] = 0;

  return true;
}

bool EqColoring::pruneFF(int color){
  gf = prevGraphsFF[color - 1].g;
  pmf = pmPrevGraphsFF[color - 1];
  pmf.rf = get(ref_vertex_t(), gf);

  EdgeFord eF1, eF2;
  VertexFordIter vIt1, vIt2;
  std::vector<VertexFord> vert;

  gf.clear();

  int n = 1 + curr.uncoloredVertices + cl.nCliques * color + color + 1, nNew = n + 2;
  long flow;

  for(int i = 0; i < nNew; i++){
    vert.push_back(add_vertex(gf));
  }

  initA1(vert);
  
  initA2andA3(vert, color);
  
  int sumLB = initA4(vert, color, std::make_pair(n, nNew));

  initRespectLB(vert, std::make_pair(n, nNew), sumLB);

  flow = performEKMF(gf, vert[nNew - 2], vert[nNew - 1]);

  if(!(flow == sumLB)){
    return true;
  }else{
    removeRespectLB(vert, color, std::make_pair(n, nNew), sumLB);
    
    flow = performEKMF(gf, vert[0], vert[n - 1]);

    if(flow == curr.uncoloredVertices){
      return false;
    }else{
      return true;
    }
  }
}

long EqColoring::performEKMF(GraphFord &fg, VertexFord &vs, VertexFord &vt){
  std::vector<default_color_type> col(num_vertices(gf));
  std::vector<Traits::edge_descriptor> pred(num_vertices(gf));

  Traits::vertex_descriptor s = vs, t = vt;
  
  long flow = edmonds_karp_max_flow(gf, s, t, pmf.c, pmf.rc, pmf.re, &col[0], &pred[0]);

  return flow;
}

bool EqColoring::useNewIndepCliques(bool sBetterClique){
  vertexIter vIt1, vIt2;
  Graph tmpG;
  Cliques tmpCl = cl;
  std::vector< std::vector<Vertex> > tmpIndClq = indClq;

  if(sBetterClique == true){
    copy_graph(g, tmpG);
  }

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    if(pm.cl[*vIt1] > 1){
      pm.cl[*vIt1] = 0;
    }
  }

  setClique(startClique.size(), 1, false);
            
  findIndepCliques(indClq, true, true);
  
  if(sBetterClique == true && cl.nodesInClique <= tmpCl.nodesInClique){
    g = tmpG;
    cl = tmpCl;
    indClq = tmpIndClq;

    return false;
  }

  return true;
}

bool EqColoring::updateIndepCliques(Vertex &v){
  if(pm.cl[v] > 1){
    useNewIndepCliques(false);
  }else{
    if(!useNewIndepCliques(true)){
      return false;
    }
  }

  c.newCliques++;
  curr.createNewGraphs = true;

  return true;
}

void EqColoring::updateBackupGraphs(Vertex &v, bool removeVertex){
  for(int i = b.LB; i < b.UB; i++){
    updateBackupGraphsHelp(v, i, removeVertex);
  }
}

void EqColoring::updateBackupGraphsHelp(Vertex &v, int i, bool removeVertex){
  gf = prevGraphsFF[i - 1].g;
  pmf = pmPrevGraphsFF[i - 1];
  pmf.rf = get(ref_vertex_t(), gf);

  EdgeFord eF1;

  for(int j = 1; j <= prevGraphsFF[i - 1].uncoloredVertices; j++){
    if(prevGraphsFF[i - 1].vert[j] == v){
      eF1 = edge(prevGraphsFF[i - 1].vert[0], prevGraphsFF[i - 1].vert[j], gf).first;

      if(removeVertex){
        pmf.c[eF1] = 0;
      }else{
        pmf.c[eF1] = 1;
      }

      break;
    } 
  }
}

bool EqColoring::nodeClique(){
  if(curr.uncoloredVertices == 0){
    std::cout << "neue bound = " << curr.nColors << std::endl;
    b.UB = curr.nColors;

    return true;
  }

  Vertex v = passVSS();

  for(int i = 1; i <= std::min(b.UB - 1, curr.nColors + 1); i++){
    if(pm.fbc[v][i - 1] == 0){
      if(parm.n >= (curr.M - 1) * std::max(curr.nColors, b.LB) + curr.T){
        if(!pruneFF()){
          c.visitedNodes++; 

          colorVertex(v, i);

          updateIndepCliques(v);

          //if(!curr.createNewGraphs){
            //updateBackupGraphs(v, true); 
          //}

          nodeClique(); 
          
          if(bt.status){
            if(bt.toRank == curr.rank){
              bt.status = false;
            }else if(bt.toRank > curr.rank){
              uncolorVertex(v);
              return true;
            }
          }

          uncolorVertex(v);
        }
      }
    }
  }

  checkForBacktracking(v);

  return false;
}

bool EqColoring::dsatur(){
  if(node()){
    return true;
  }else{
    return false;
  }
}

bool EqColoring::dsaturClique(){
  if(nodeClique()){
    return true;
  }else{
    return false;
  }
}

bool EqColoring::pruningRulePaper(){
  return true;
}

bool EqColoring::pruningRuleFF(){
  return true;
}

bool EqColoring::checkEquitability() const{
  for(unsigned int i = 0; i < curr.nColors; i++){
    for(unsigned int j = i+1; j < curr.nColors; j++){
      if(std::abs(cc.n[i] - cc.n[j]) > 1){
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

int EqColoring::naiveUB(){
  vertexIter vIt1, vIt2;
  Vertex v;
  bool haveSwapped = false;
  std::pair<int, int> cMinMax;

  while(!checkEquitability()){
    cMinMax = findMinMaxColorClass(INT_MAX, INT_MIN);

    for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
      v = *vIt1;

      if(cc.n[pm.c[*vIt1]-1] == cMinMax.second){
        for(int j = 0; j < curr.nColors; j++){
          if(cc.n[j] == cMinMax.first){
            haveSwapped = swapNodeToColor(v, j + 1);
            
            if(haveSwapped){
              break;
            }
          }
        }

        if(haveSwapped){
          break;
        }
      }
    }

    if(!haveSwapped){
      curr.nColors++;
      cc.n[pm.c[v]-1]--;
      pm.c[v] = curr.nColors;
      cc.n[curr.nColors-1]++;
    }
  }

  return curr.nColors;
}

bool EqColoring::swapNodeToColor(Vertex v, int color){
  int tmpColor;
  
  tmpColor = pm.c[v];
  pm.c[v] = color;

  if(!checkColoring()){
    pm.c[v] = tmpColor;

    return false;
  }else{
    cc.n[tmpColor - 1]--;
    cc.n[color - 1]++;

    return true;
  }
}

std::pair<int, int> EqColoring::findMinMaxColorClass(int cMin, int cMax){
  std::pair<int, int> cMinMax = std::make_pair(cMin, cMax);

  for(int i = 0; i < curr.nColors; i++){
    if(cc.n[i] > cMinMax.second){
      cMinMax.second = cc.n[i];
    }
    
    if(cc.n[i] < cMinMax.first){
      cMinMax.first = cc.n[i];
    }
  }

  return cMinMax;
}

int EqColoring::calcUB(){
  Graph tmpG;
  Current tmp;
  Colors cc_tmp;
  int UB;

  copy_graph(g, tmpG);
  tmp = curr;
  cc_tmp = cc;

  greedyColoring(g);
  UB = naiveUB();

  g = tmpG;
  curr = tmp;
  cc = cc_tmp;

  return UB;
}
