#include "Coloring.hpp"
#include "EqColoring.hpp"

void usage(){
  std::cout << "Usage: ./eqDsatur [N | R] nPruningRule timeLimit threshold [FILE | nRandomGraphs n p]" << std::endl;
}

int main(int argc, char** argv){
  Parameters p;

  if(!Input::readInputArgs(argc, argv, p)){
    usage();

    exit(-1);
  }

  Coloring c(p);

  return 0;
}
