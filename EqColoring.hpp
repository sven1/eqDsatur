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

    bool node(Graph &g);

    bool pruningRulePaper();
    bool pruningRuleFF();

    bool checkEquitability() const;
    bool checkEqColoring() const;
};

#endif
