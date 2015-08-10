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
    Colors cc;

    bool initPropMap();
    bool initNeighbours();
    bool initResize();

    bool initVar();

    bool setBounds(int LB, int UB);
    bool setCurr(int c, int r, Vertex node, int uncoloredVertices);
    bool setBacktracking(bool status, Vertex node);

    bool findMaxClique(std::vector<Vertex> &clq, Vertex v, bool uncolored, bool inNoOtherClique);

    bool compareDegree(Vertex v, Vertex w);

    bool putNodeWithParm(Vertex v, std::vector<Vertex> &tmp, bool uncolored, bool inOtherClique);

    void printNeighbours(const Vertex &v) const;
    void printVertexInfo(const Vertex &v) const;
    void printFBC(const Vertex &v) const;

    bool checkColoring(const Vertex &v) const;
    bool checkClique(const std::vector<Vertex> &clq) const;

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
    void printClique(const std::vector<Vertex> &clq) const;
    
    bool findMaxClique(std::vector<Vertex> &clq, bool uncolored, bool inNoOtherClique);

    bool checkColoring() const;
};

#endif
