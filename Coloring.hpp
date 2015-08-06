#ifndef COLORING_HPP
#define COLORING_HPP

#include "structs.hpp"
#include "Input.hpp"
#include "Heuristic.hpp"

class Coloring {
  private:
    Graph g;
    Bounds b;
    PropertyMap pm;
    Time t;
    Count c;
    std::vector<Vertex> startClique;
    Current curr;
    Parameters parm;
    Backtracking bt;

    void printNeighbours(const Vertex &v) const;
    void printVertexInfo(const Vertex &v) const;
    void printFBC(const Vertex &v) const;

    bool checkColoring(const Vertex &v) const;

  protected:
    Colors cls;

  public:
    Coloring();
    Coloring(const Parameters &parm);
    Coloring(const Graph &g);
    Coloring(const Graph &g, const Parameters &parm);
  
    ~Coloring();

    // bounds
    long naiveUB(const Graph &g);
    long eqlLB(const Graph &g);

    bool node(Graph &g);

    void printVertexInfo() const;
    void printAdjMatrix() const;
    void printFBC() const;
    
    bool checkColoring() const;
};

#endif
