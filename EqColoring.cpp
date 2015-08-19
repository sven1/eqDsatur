#include "EqColoring.hpp"

EqColoring::EqColoring() : Coloring(){
}

EqColoring::EqColoring(const Parameters &parm) : Coloring(parm){
  setBounds(b.LB, calcUB());

  printAll();
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
  bool haveSwapped;
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
