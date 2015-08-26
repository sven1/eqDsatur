#ifndef EQCOLORING_HPP
#define EQCOLORING_HPP

#include "structs.hpp"
#include "Coloring.hpp"

class EqColoring : Coloring{
  private:
    std::vector<prevGraphs> prevGraphsFF;
    std::vector<PropertyMapFF> pmPrevGraphsFF;

    GraphFord gf;
    PropertyMapFF pmf;

    bool initPrevGraphsFF();

    bool initA1(std::vector<VertexFord> &vert);
    bool initA2andA3(std::vector<VertexFord> &vert, int color);
    int initA4(std::vector<VertexFord> &vert, int color, std::pair<int, int> counts);

    bool initRespectLB(std::vector<VertexFord> &vert, std::pair<int, int> counts, int sumLB);
    bool removeRespectLB(std::vector<VertexFord> &vert, int color, std::pair<int, int> counts, int sumLB);
    
    long performEKMF(GraphFord &gf, VertexFord &vs, VertexFord &vt);

    bool updateIndepCliques(Vertex &v);
    void updateBackupGraphs(Vertex &v, bool removeVertex);
    void updateBackupGraphsHelp(Vertex &v, int i, bool removeVertex);

  public:
    EqColoring();
    EqColoring(const Parameters &parm);
    EqColoring(const Graph &g);
    EqColoring(const Graph &g, const Parameters &parm);

    ~EqColoring();

    bool node();
    bool nodeClique();

    bool dsatur();
    bool dsaturClique();

    bool pruningRulePaper();
    bool pruningRuleFF();

    bool checkEquitability() const;
    bool checkEqColoring() const;

    bool useNewIndepCliques(bool sBetterClique);

    bool pruneFF();
    bool pruneFF(int color);
    
    int calcUB();
    int naiveUB();

    bool swapNodeToColor(Vertex v, int color);
    std::pair<int, int> findMinMaxColorClass(int cMin, int cMax);
};

#endif
