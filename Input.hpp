#include "boost.hpp"
#include <string>

class Input{
  private:
    static bool readGraph(const std::string &filename, const std::string &sType , Graph &g);
    static bool readGraphCol(const std::string &filename, Graph &g);

  public:
    static bool readGraph(const std::string &filename, Graph &g);
    static bool readClique(Graph &g);
};
