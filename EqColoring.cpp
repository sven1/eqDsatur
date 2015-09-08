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
  gf.resize(b.UB);
  pmf.resize(b.UB);

  for(unsigned int i = 0; i < pmf.size(); i++){
    pmf[i].c = get(edge_capacity, gf[i].g);
    pmf[i].re = get(edge_reverse, gf[i].g);
    pmf[i].rc = get(edge_residual_capacity, gf[i].g);
    pmf[i].rf = get(ref_vertex_t(), gf[i].g);
  }

  initBackupGraphs();

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

bool EqColoring::initA1(int j){
  EdgeFord eF1, eF2;

  for(int i = 1; i <= curr.uncoloredVertices; i++){
    eF1 = add_edge(gf[j-1].vert[0], gf[j-1].vert[i], gf[j-1].g).first;
    eF2 = add_edge(gf[j-1].vert[i], gf[j-1].vert[0], gf[j-1].g).first;

    pmf[j-1].re[eF1] = eF2;
    pmf[j-1].re[eF2] = eF1;

    pmf[j-1].c[eF1] = 1;
    pmf[j-1].c[eF2] = 0;
  }

  return true;
}

bool EqColoring::initA2andA3(int l){
  EdgeFord eF1, eF2;
  VertexFordIter vIt1, vIt2;

  int aUVPos = 1;
  int tmpIndex, colorPos;
  int color = l;

  for(unsigned int i = 0; i < indClq.size(); i++){
    for(unsigned int j = 0; j < indClq[i].size(); j++){
      pmf[l-1].rf[gf[l-1].vert[aUVPos]] = indClq[i][j];

      for(int k = 0; k < color; k++){
        if(pm.fbc[indClq[i][j]][k] == 0){
          tmpIndex = curr.uncoloredVertices + (k + 1) + i * color; 

          eF1 = add_edge(gf[l-1].vert[aUVPos], gf[l-1].vert[tmpIndex], gf[l-1].g).first;
          eF2 = add_edge(gf[l-1].vert[tmpIndex], gf[l-1].vert[aUVPos], gf[l-1].g).first;

          pmf[l-1].re[eF1] = eF2;
          pmf[l-1].re[eF2] = eF1;

          pmf[l-1].c[eF1] = 1;
          pmf[l-1].c[eF2] = 0;

          //zur Farbe
          colorPos = curr.uncoloredVertices + (k + 1) + cl.nCliques * color; 

          eF1 = add_edge(gf[l-1].vert[tmpIndex], gf[l-1].vert[colorPos], gf[l-1].g).first;
          eF2 = add_edge(gf[l-1].vert[colorPos], gf[l-1].vert[tmpIndex], gf[l-1].g).first;

          pmf[l-1].re[eF1] = eF2;
          pmf[l-1].re[eF2] = eF1;

          pmf[l-1].c[eF1] = 1;
          pmf[l-1].c[eF2] = 0;
        }
      }

      aUVPos++;
    }
  }

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    if(pm.c[*vIt1] == 0 && pm.cl[*vIt1] == 0){
      pmf[l-1].rf[gf[l-1].vert[aUVPos]] = *vIt1;

      for(int k = 0; k < color; k++){
        if(pm.fbc[*vIt1][k] == 0){
          tmpIndex = curr.uncoloredVertices + (k + 1) + (cl.nCliques - 1) * color; 

          eF1 = add_edge(gf[l-1].vert[aUVPos], gf[l-1].vert[tmpIndex], gf[l-1].g).first;
          eF2 = add_edge(gf[l-1].vert[tmpIndex], gf[l-1].vert[aUVPos], gf[l-1].g).first;

          pmf[l-1].re[eF1] = eF2;
          pmf[l-1].re[eF2] = eF1;

          pmf[l-1].c[eF1] = 1;
          pmf[l-1].c[eF2] = 0;

          colorPos = curr.uncoloredVertices + (k + 1) + cl.nCliques * color; 

          eF1 = add_edge(gf[l-1].vert[tmpIndex], gf[l-1].vert[colorPos], gf[l-1].g).first;
          eF2 = add_edge(gf[l-1].vert[colorPos], gf[l-1].vert[tmpIndex], gf[l-1].g).first;

          pmf[l-1].re[eF1] = eF2;
          pmf[l-1].re[eF2] = eF1;

          pmf[l-1].c[eF1] = curr.uncoloredVertices - (cl.nodesInClique - startClique.size());
          pmf[l-1].c[eF2] = 0;
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

int EqColoring::initA4(int l){
  int sumLB = 0, rU, rL, colorPos, n = gf[l-1].n, nNew = gf[l-1].nNew, color = l;
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

    colorPos = gf[l-1].uncoloredVertices + i + cl.nCliques * color; 

    eF1 = add_edge(gf[l-1].vert[colorPos], gf[l-1].vert[n - 1], gf[l-1].g).first;
    eF2 = add_edge(gf[l-1].vert[n - 1], gf[l-1].vert[colorPos], gf[l-1].g).first;

    pmf[l-1].re[eF1] = eF2;
    pmf[l-1].re[eF2] = eF1;

    pmf[l-1].c[eF1] = rU - rL;
    pmf[l-1].c[eF2] = 0;

    eF1 = add_edge(gf[l-1].vert[colorPos], gf[l-1].vert[nNew - 1], gf[l-1].g).first;
    eF2 = add_edge(gf[l-1].vert[nNew - 1], gf[l-1].vert[colorPos], gf[l-1].g).first;

    pmf[l-1].re[eF1] = eF2;
    pmf[l-1].re[eF2] = eF1;

    pmf[l-1].c[eF1] = rL;
    pmf[l-1].c[eF2] = 0;
  }

  return sumLB;
}

bool EqColoring::initRespectLB(int l, int sumLB){
  EdgeFord eF1, eF2;
  int n = gf[l-1].n, nNew = gf[l-1].nNew;

  eF1 = add_edge(gf[l-1].vert[nNew - 2], gf[l-1].vert[n - 1], gf[l-1].g).first;
  eF2 = add_edge(gf[l-1].vert[n - 1], gf[l-1].vert[nNew - 2], gf[l-1].g).first;

  pmf[l-1].re[eF1] = eF2;
  pmf[l-1].re[eF2] = eF1;

  pmf[l-1].c[eF1] = sumLB;
  pmf[l-1].c[eF2] = 0;

  eF1 = add_edge(gf[l-1].vert[n - 1], gf[l-1].vert[0], gf[l-1].g).first;
  eF2 = add_edge(gf[l-1].vert[0], gf[l-1].vert[n - 1], gf[l-1].g).first;

  pmf[l-1].re[eF1] = eF2;
  pmf[l-1].re[eF2] = eF1;

  pmf[l-1].c[eF1] = INT_MAX;
  pmf[l-1].c[eF2] = 0;

  return true;
}

bool EqColoring::removeRespectLB(int l, int sumLB){
  int rL, rU, colorPos, n = gf[l-1].n, nNew = gf[l-1].nNew, color = l;

  EdgeFord eF1, eF2;

  for(int i = 1; i <= color; i++){
    rU = parm.n / color + 1;
    rL = parm.n / color;

    rU -= cc.n[i-1];
    rL -= cc.n[i-1];

    if(rL < 0){
      rL = 0;
    }

    colorPos = gf[l-1].uncoloredVertices + i + cl.nCliques * color; 

    eF1 = edge(gf[l-1].vert[colorPos], gf[l-1].vert[n - 1], gf[l-1].g).first;
    eF2 = edge(gf[l-1].vert[n - 1], gf[l-1].vert[colorPos], gf[l-1].g).first;

    pmf[l-1].c[eF1] = rU;
    pmf[l-1].c[eF2] = 0;

    pmf[l-1].rc[eF1] = rU - rL;
    pmf[l-1].rc[eF1] = 0;

    eF1 = edge(gf[l-1].vert[colorPos], gf[l-1].vert[nNew - 1], gf[l-1].g).first;
    eF2 = edge(gf[l-1].vert[nNew - 1], gf[l-1].vert[colorPos], gf[l-1].g).first;

    pmf[l-1].c[eF1] = 0;
    pmf[l-1].c[eF2] = 0;

    pmf[l-1].rc[eF1] = 0;
    pmf[l-1].rc[eF1] = 0;
  }

  eF1 = edge(gf[l-1].vert[n - 1], gf[l-1].vert[0], gf[l-1].g).first;
  eF2 = edge(gf[l-1].vert[0], gf[l-1].vert[n - 1], gf[l-1].g).first;

  pmf[l-1].c[eF1] = 0;
  pmf[l-1].c[eF2] = 0;

  pmf[l-1].rc[eF1] = 0;
  pmf[l-1].rc[eF1] = 0;

  return true;
}

void EqColoring::initBackupGraphs(){
  for(int i = b.LB; i < b.UB; i++){
    gf[i-1].g.clear();
    gf[i-1].vert.clear();

    gf[i - 1].uncoloredVertices = curr.uncoloredVertices;
    gf[i - 1].n = 1 + curr.uncoloredVertices + cl.nCliques * i + i + 1;
    gf[i - 1].nNew = gf[i - 1].n + 2;

    for(int j = 0; j < gf[i-1].nNew; j++){
      gf[i-1].vert.push_back(add_vertex(gf[i-1].g));
    }
    
    initA1(i);
  
    initA2andA3(i);
    
    int sumLB = initA4(i);

    initRespectLB(i, sumLB);

    gf[i-1].sumLB = sumLB;
  }
}

void EqColoring::resetCap(int i){
  EdgesOutFordIter ei1, ei2;
  VertexFordIter it1, it2;

  for(tie(it1, it2) = vertices(gf[i-1].g); it1 != it2; it1++){
    for(tie(ei1, ei2) = out_edges(*it1, gf[i-1].g); ei1 != ei2; ei1++){
      pmf[i-1].rc[*ei1] = 0; 
    }
  }
}

bool EqColoring::pruneFF(int color){
  EdgeFord eF1, eF2;
  VertexFordIter vIt1, vIt2;

  long flow;

  if(curr.createNewGraphs){
    initBackupGraphs();

    curr.createNewGraphs = false;
  }else{
    resetCap(color);
  }

  flow = performEKMF(color, gf[color-1].vert[gf[color-1].nNew - 2], gf[color-1].vert[gf[color-1].nNew - 1]);

  if(!(flow == gf[color-1].sumLB)){
    return true;
  }else{
    removeRespectLB(color, gf[color-1].sumLB);
    
    flow = performEKMF(color, gf[color-1].vert[0], gf[color-1].vert[gf[color-1].n - 1]);

    if(flow == curr.uncoloredVertices){
      return false;
    }else{
      return true;
    }
  }
}

long EqColoring::performEKMF(int l, VertexFord &vs, VertexFord &vt){
  std::vector<default_color_type> col(num_vertices(gf[l-1].g));
  std::vector<Traits::edge_descriptor> pred(num_vertices(gf[l-1].g));

  Traits::vertex_descriptor s = vs, t = vt;
  
  long flow = edmonds_karp_max_flow(gf[l-1].g, s, t, pmf[l-1].c, pmf[l-1].rc, pmf[l-1].re, &col[0], &pred[0]);

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
  curr.rankNC = curr.rank;
  curr.createNewGraphs = true;

  return true;
}

void EqColoring::updateBackupGraphs(Vertex &v, int k, bool removeVertex){
  for(int i = b.LB; i < b.UB; i++){
    updateBackupGraphsHelp(v, i, k, removeVertex);
  }
}

void EqColoring::updateBackupGraphsHelp(Vertex &v, int i, int k, bool removeVertex){
  EdgeFord eF1, eF2;

  for(int j = 1; j <= gf[i - 1].uncoloredVertices; j++){
    if(pmf[i-1].rf[gf[i - 1].vert[j]] == (int) v){
      eF1 = edge(gf[i - 1].vert[0], gf[i - 1].vert[j], gf[i-1].g).first;

      if(removeVertex){
        pmf[i-1].c[eF1] = 0;

        if(pm.cl[v] != 0){
          int colorPos1 = gf[i-1].uncoloredVertices + k + cl.nCliques * i; 
          int tmpIndex = gf[i-1].uncoloredVertices + k + (pm.cl[v] - 2) * i; 

          eF1 = edge(gf[i-1].vert[tmpIndex], gf[i-1].vert[colorPos1], gf[i-1].g).first;
          eF2 = edge(gf[i-1].vert[colorPos1], gf[i-1].vert[tmpIndex], gf[i-1].g).first;

          pmf[i-1].c[eF1] = 0;
          pmf[i-1].c[eF2] = 0;
        }
      }else{
        pmf[i-1].c[eF1] = 1;

        if(pm.cl[v] != 0){
          int colorPos1 = gf[i-1].uncoloredVertices + k + cl.nCliques * i; 
          int tmpIndex = gf[i-1].uncoloredVertices + k + (pm.cl[v] - 2) * i; 

          eF1 = edge(gf[i-1].vert[tmpIndex], gf[i-1].vert[colorPos1], gf[i-1].g).first;
          eF2 = edge(gf[i-1].vert[colorPos1], gf[i-1].vert[tmpIndex], gf[i-1].g).first;

          pmf[i-1].c[eF1] = 1;
          pmf[i-1].c[eF2] = 0;
        }
      }

      break;
    } 
  }

  if(removeVertex){
    eF1 = edge(gf[i-1].vert[gf[i-1].n - 1], gf[i-1].vert[0], gf[i-1].g).first;
    eF2 = edge(gf[i-1].vert[0], gf[i-1].vert[gf[i-1].n - 1], gf[i-1].g).first;

    pmf[i-1].c[eF1] = INT_MAX;
    pmf[i-1].c[eF2] = 0;

    eF1 = edge(gf[i-1].vert[gf[i-1].nNew - 2], gf[i-1].vert[gf[i-1].n - 1], gf[i-1].g).first;
    eF2 = edge(gf[i-1].vert[gf[i-1].n - 1], gf[i-1].vert[gf[i-1].nNew - 2], gf[i-1].g).first;

    gf[i-1].sumLB--;

    pmf[i-1].c[eF1]--;
    pmf[i-1].c[eF2] = 0;

    int color = i, rU, rL, colorPos;

    for(int z = 1; z <= color; z++){
      rU = parm.n / color + 1;
      rL = parm.n / color;

      rU -= cc.n[z-1];
      rL -= cc.n[z-1];

      if(rL < 0){
        rL = 0;
      }

      colorPos = gf[i-1].uncoloredVertices + z + cl.nCliques * color; 

      eF1 = edge(gf[i-1].vert[colorPos], gf[i-1].vert[gf[i-1].n - 1], gf[i-1].g).first;
      eF2 = edge(gf[i-1].vert[gf[i-1].n - 1], gf[i-1].vert[colorPos], gf[i-1].g).first;

      pmf[i-1].c[eF1] = rU - rL;
      pmf[i-1].c[eF2] = 0;

      eF1 = edge(gf[i-1].vert[colorPos], gf[i-1].vert[gf[i-1].nNew - 1], gf[i-1].g).first;
      eF2 = edge(gf[i-1].vert[gf[i-1].nNew - 1], gf[i-1].vert[colorPos], gf[i-1].g).first;

      pmf[i-1].c[eF1] = rL;
      pmf[i-1].c[eF2] = 0;
    }
  }else{
    eF1 = edge(gf[i-1].vert[gf[i-1].n - 1], gf[i-1].vert[0], gf[i-1].g).first;
    eF2 = edge(gf[i-1].vert[0], gf[i-1].vert[gf[i-1].n - 1], gf[i-1].g).first;

    pmf[i-1].c[eF1] = INT_MAX;
    pmf[i-1].c[eF2] = 0;

    eF1 = edge(gf[i-1].vert[gf[i-1].nNew - 2], gf[i-1].vert[gf[i-1].n - 1], gf[i-1].g).first;
    eF2 = edge(gf[i-1].vert[gf[i-1].n - 1], gf[i-1].vert[gf[i-1].nNew - 2], gf[i-1].g).first;

    gf[i-1].sumLB++;

    pmf[i-1].c[eF1]++;
    pmf[i-1].c[eF2] = 0;

    int color = i, rU, rL, colorPos;

    for(int z = 1; z <= color; z++){
      rU = parm.n / color + 1;
      rL = parm.n / color;

      rU -= cc.n[z-1];
      rL -= cc.n[z-1];

      if(rL < 0){
        rL = 0;
      }

      colorPos = gf[i-1].uncoloredVertices + z + cl.nCliques * color; 

      eF1 = edge(gf[i-1].vert[colorPos], gf[i-1].vert[gf[i-1].n - 1], gf[i-1].g).first;
      eF2 = edge(gf[i-1].vert[gf[i-1].n - 1], gf[i-1].vert[colorPos], gf[i-1].g).first;

      pmf[i-1].c[eF1] = rU - rL;
      pmf[i-1].c[eF2] = 0;

      eF1 = edge(gf[i-1].vert[colorPos], gf[i-1].vert[gf[i-1].nNew - 1], gf[i-1].g).first;
      eF2 = edge(gf[i-1].vert[gf[i-1].nNew - 1], gf[i-1].vert[colorPos], gf[i-1].g).first;

      pmf[i-1].c[eF1] = rL;
      pmf[i-1].c[eF2] = 0;
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

          if(!curr.createNewGraphs){
            updateBackupGraphs(v, i, true); 
          }

          nodeClique(); 

          if(bt.status){
            if(bt.toRank == curr.rank){
              bt.status = false;
            }else if(bt.toRank > curr.rank){
              uncolorVertex(v);

              checkUpdateBackupGraphs(v, i);

              return true;
            }
          }

          uncolorVertex(v);

          checkUpdateBackupGraphs(v, i);
        }
      }
    }
  }

  checkForBacktracking(v);

  return false;
}

void EqColoring::checkUpdateBackupGraphs(Vertex &v, int i){
  if(curr.rank >= curr.rankNC){
    updateBackupGraphs(v, i, false); 
  }else if(curr.rank < curr.rankNC){
    curr.createNewGraphs = true;
    curr.rankNC = curr.rank;
  }
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
