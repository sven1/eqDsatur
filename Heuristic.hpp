#include "boost.hpp"
#include "structs.hpp"

class Heuristic{
  public:
    static bool findIndepCliques(Graph &g, PropertyMap &pm);
    static bool findMaxClique(Graph &g, PropertyMap &pm);
};
