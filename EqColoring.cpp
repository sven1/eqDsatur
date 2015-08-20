#include "EqColoring.hpp"

EqColoring::EqColoring() : Coloring(){
}

EqColoring::EqColoring(const Parameters &parm) : Coloring(parm){
  setBounds(b.LB, calcUB());

  initPrevGraphsFF();

  dsaturClique();

  printAll();
}

bool EqColoring::initPrevGraphsFF(){
  prevGraphsFF.resize(b.UB);
  pmPrevGraphsFF.resize(b.UB);

  for(unsigned int i = 0; i < pmPrevGraphsFF.size(); i++){
    pmPrevGraphsFF[i].c = get(edge_capacity, prevGraphsFF[i].g);
    pmPrevGraphsFF[i].re = get(edge_reverse, prevGraphsFF[i].g);
    pmPrevGraphsFF[i].rc = get(edge_residual_capacity, prevGraphsFF[i].g);
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
  for(unsigned int i = b.LB; i <= curr.nColors; i++){
    if(!pruneFF(i)){
      return false;
    }
  }

  return true;
}

bool EqColoring::pruneFF(int color){
  gf = prevGraphsFF[color - 1].g;
  pmf = pmPrevGraphsFF[color - 1];

  Traits::vertex_descriptor s, t, sNew, tNew;
  EdgeFord eF1, eF2;
  VertexFordIter vIt1, vIt2;
  std::vector<VertexFord> vert;

  gf.clear();

  int n = 1 + curr.uncoloredVertices + cl.nCliques * color + color + 1;
  int nNew = n + 2;

  for(int i = 0; i < nNew; i++){
    vert.push_back(add_vertex(gf));
  }

  for(int i = 1; i <= curr.uncoloredVertices; i++){
    eF1 = add_edge(vert[0], vert[i], gf).first;
    eF2 = add_edge(vert[i], vert[0], gf).first;

    pmf.re[eF1] = eF2;
    pmf.re[eF2] = eF1;

    pmf.c[eF1] = 1;
    pmf.c[eF2] = 0;
  }

  int aUVPos = 1;
  int tmpIndex, colorPos;

  for(unsigned int i = 0; i < indClq.size(); i++){
    for(unsigned int j = 0; j < indClq[i].size(); j++){
      for(int k = 0; k < color; k++){
        if(pm.fbc[indClq[i][j]][k] == 0){
          tmpIndex = curr.uncoloredVertices + (k + 1) + i * (k + 1); 

          eF1 = add_edge(vert[aUVPos], vert[tmpIndex], gf).first;
          eF2 = add_edge(vert[tmpIndex], vert[aUVPos], gf).first;

          pmf.re[eF1] = eF2;
          pmf.re[eF2] = eF1;

          pmf.c[eF1] = 1;
          pmf.c[eF2] = 0;

          //zur Farbe
          colorPos = curr.uncoloredVertices + (k + 1) + (cl.nCliques + 1) * (k + 1); 

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
      for(int k = 0; k < color; k++){
        if(pm.fbc[*vIt1][k] == 0){
          tmpIndex = curr.uncoloredVertices + (k + 1) + cl.nCliques * (k + 1); 

          eF1 = add_edge(vert[aUVPos], vert[tmpIndex], gf).first;
          eF2 = add_edge(vert[tmpIndex], vert[aUVPos], gf).first;

          pmf.re[eF1] = eF2;
          pmf.re[eF2] = eF1;

          pmf.c[eF1] = 1;
          pmf.c[eF2] = 0;

          //zur Farbe
          colorPos = curr.uncoloredVertices + (k + 1) + (cl.nCliques + 1) * (k + 1); 

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

  int sumLB = 0, rU, rL;

  for(int i = 1; i <= color; i++){
    rU = parm.n / color + 1;
    rL = parm.n / color;

    rU -= cc.n[i-1];
    rL -= cc.n[i-1];

    if(rL < 0){
      rL = 0;
    }

    sumLB += rL;

    colorPos = curr.uncoloredVertices + i + (cl.nCliques + 1) * i; 

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

  std::vector<default_color_type> col(num_vertices(gf));
  std::vector<Traits::edge_descriptor> pred(num_vertices(gf));
  s = vert[nNew - 2], t = vert[nNew - 1];
  
  long flow = edmonds_karp_max_flow(gf, s, t, pmf.c, pmf.rc, pmf.re, &col[0], &pred[0]);

  if(!(flow == sumLB)){
    return true;
  }else{
    for(int i = 1; i <= color; i++){
      rU = parm.n / color + 1;
      rL = parm.n / color;

      rU -= cc.n[i-1];
      rL -= cc.n[i-1];

      if(rL < 0){
        rL = 0;
      }

      colorPos = curr.uncoloredVertices + i + (cl.nCliques + 1) * i; 

      eF1 = add_edge(vert[colorPos], vert[n - 1], gf).first;
      eF2 = add_edge(vert[n - 1], vert[colorPos], gf).first;

      pmf.re[eF1] = eF2;
      pmf.re[eF2] = eF1;

      pmf.c[eF1] = rU;
      pmf.c[eF2] = 0;

      pmf.rc[eF1] = rU - rL;
      pmf.rc[eF1] = 0;

      eF1 = add_edge(vert[colorPos], vert[nNew - 1], gf).first;
      eF2 = add_edge(vert[nNew - 1], vert[colorPos], gf).first;

      pmf.re[eF1] = eF2;
      pmf.re[eF2] = eF1;

      pmf.c[eF1] = 0;
      pmf.c[eF2] = 0;

      pmf.rc[eF1] = 0;
      pmf.rc[eF1] = 0;
    }
    
    eF1 = add_edge(vert[nNew - 2], vert[n - 1], gf).first;
    eF2 = add_edge(vert[n - 1], vert[nNew - 2], gf).first;

    pmf.re[eF1] = eF2;
    pmf.re[eF2] = eF1;

    pmf.c[eF1] = sumLB;
    pmf.c[eF2] = 0;

    pmf.rc[eF1] = 0;
    pmf.rc[eF1] = 0;
    
    eF1 = add_edge(vert[n - 1], vert[0], gf).first;
    eF2 = add_edge(vert[0], vert[n - 1], gf).first;

    pmf.re[eF1] = eF2;
    pmf.re[eF2] = eF1;

    pmf.c[eF1] = 0;
    pmf.c[eF2] = 0;

    pmf.rc[eF1] = 0;
    pmf.rc[eF1] = 0;

    std::vector<default_color_type> col2(num_vertices(gf));
    std::vector<Traits::edge_descriptor> pred2(num_vertices(gf));
    sNew = vert[nNew - 2], tNew = vert[nNew - 1];

    long flow2 = edmonds_karp_max_flow(gf, sNew, tNew, pmf.c, pmf.rc, pmf.re, &col2[0], &pred2[0]);

    if(flow2 == curr.uncoloredVertices){
      return false;
    }else{
      return true;
    }
  }
}

bool EqColoring::useNewIndepCliques(){
  std::vector<std::vector<Vertex> > indClqTmp = indClq;
  Cliques clTmp = cl;
  Graph tmpG;
            
  copy_graph(g, tmpG);
  setClique(0, 0, false);
            
  findIndepCliques(indClq, true, false);

  if(clTmp.nodesInClique > cl.nodesInClique){
    g = tmpG;
    cl = clTmp;
    indClq = indClqTmp;

    return false;
  }

  return true;
}

bool EqColoring::nodeClique(){
  if(curr.uncoloredVertices == 0){
    b.UB = curr.nColors;

    return true;
  }

  Vertex v = passVSS();


  for(int i = 1; i <= std::min(b.UB - 1, curr.nColors + 1); i++){
    if(pm.fbc[v][i - 1] == 0){
      if(parm.n >= (curr.M - 1) * std::max(curr.nColors, b.LB) + curr.T){
        if(!pruneFF()){
          colorVertex(v, i);

          if(pm.cl[v] != 0){
            useNewIndepCliques();
          }

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
