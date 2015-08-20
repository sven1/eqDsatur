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
};

struct prevGraphs {
  GraphFord g;
  Vertex dlNode;
  long dlColor;
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

struct Cliques {
  long nodesInClique;
  long nCliques;
  bool newClique;
};

struct Colors {
  std::vector<int> n;
};

#endif
