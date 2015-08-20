#ifndef EQCOLORING_HPP
#define EQCOLORING_HPP

#include "structs.hpp"
#include "Coloring.hpp"

class EqColoring : Coloring{
  private:
    std::vector<prevGraphs> bpGraphs;  

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
    
    int calcUB();
    int naiveUB();

    bool swapNodeToColor(Vertex v, int color);
    std::pair<int, int> findMinMaxColorClass(int cMin, int cMax);
};

#endif
