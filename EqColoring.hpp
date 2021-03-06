#ifndef EQCOLORING_HPP
#define EQCOLORING_HPP

#include "structs.hpp"
#include "Coloring.hpp"

class EqColoring : Coloring{
  private:
    std::vector<prevGraphs> gf;
    std::vector<PropertyMapFF> pmf;

    bool initPrevGraphsFF();
    void initBackupGraphs();

    bool initA1(int j);
    bool initA2andA3(int l);
    int initA4(int l);

    bool initRespectLB(int l, int sumLB);
    bool removeRespectLB(int l, int sumLB);
    
    long performEKMF(int l, VertexFord &vs, VertexFord &vt);

    bool updateIndepCliques(Vertex &v);
    void updateBackupGraphs(Vertex &v, int k, bool removeVertex);
    void updateBackupGraphsHelp(Vertex &v, int i, int k, bool removeVertex);
    void checkUpdateBackupGraphs(Vertex &v, int i);

    void resetCap(int i);

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
