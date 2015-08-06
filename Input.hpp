#include "boost.hpp"

class Input{
  public:
    static bool readGraph(Graph &g);
    static bool readClique(Graph &g);

    static bool constructRandomGraph(Graph &g, int n, double p);
};
