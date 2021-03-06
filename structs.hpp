#ifndef STRUCTS_HPP
#define STRUCTS_HPP

#include "boost.hpp"
#include <string>

struct Current {
  long color;
  long nColors;
  long rank;
  Vertex node;
  long uncoloredVertices;
  int T;
  int M;
  bool createNewGraphs;
  int rankNC;
};

struct prevGraphs {
  GraphFord g;
  Vertex dlNode;
  long dlColor;
  int uncoloredVertices;
  std::vector<VertexFord> vert;
  int n;
  int nNew;
  int sumLB;
};

struct Bounds {
  long UB;
  long LB;
};

struct Backtracking {
  bool status;
  Vertex toNode;
  int toRank;
};

struct Count {
  long visitedNodes;
  long newCliques;
  long nFF;
  long backtracks;
};

struct Parameters {
  long nPruningRule;
  long timeLimit;
  long threshold;
  std::string ressource;
  long nRandomGraphs;
  long n;
  double p;
  char variant;
  char eqDsatur;
};

struct Time {
  long timeFF;
  long timeFindIndCliques;
  long timeTotal;
};

struct PropertyMap {
  NachbarMap n;
  CliqueMap cl;
  ColorMap c;
  IndexMap i;
  FBCMap fbc;
  GradMap g;
  RankMap r;
};

struct PropertyMapFF {
  CapacityMap c;
  ReverseEdgeMap re;
  ResidualCapacityMap rc;
  RefVertexMap rf;
};

struct Cliques {
  long nodesInClique;
  long nCliques;
  bool newClique;
};

struct Colors {
  std::vector<int> n;
};

#endif
