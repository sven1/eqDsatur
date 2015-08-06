#include "Input.hpp"
#include <fstream>

bool Input::readGraph(const std::string &filename, Graph &g){
  std::string sType;
  std::size_t index = filename.find(".");

  if(index != std::string::npos){
    sType = filename.substr(index+1);
        
    if(readGraph(filename, sType, g)){
      return true;
    }else{
      return false;
    }
  }else{
    std::cout << "Not a expected *.* file format!" << std::endl;

    return false;
  }
}

bool Input::readClique(Graph &g){
  return true;
}

bool Input::readGraph(const std::string &filename, const std::string &sType, Graph &g){
  if(sType.compare("col") == 0){
    if(readGraphCol(filename, g)){
      return true;
    }else{
      std::cout << "Error: Can not read graph." << std::endl;

      return false;
    }
  }else{
    std::cout << "No suitable parser avalaible." << std::endl;

    return false;
  }
}

bool Input::readGraphCol(const std::string &filename, Graph &g){
  std::ifstream file(filename.c_str(), std::ios::in);

  char control = 0;
  std::string tmp;
  int nVertices, nEdges, fromV, toV;
  bool error;

  if(file.good()){
    while(file >> control){
      switch(control){
        case 'c':{
          getline(file, tmp);
        } break;

        case 'p':{
          file >> tmp;

          if(tmp.compare("edge") != 0){
            std::cout << "Not \"edge\" format." << std::endl;

            return false;
          }

          file >> nVertices;
          file >> nEdges;
        } break;

        case 'e':{
          file >> fromV;
          file >> toV;

          error = add_edge(fromV, toV, g).second;

          if(!error){
            std::cout << "edge already exists" << std::endl;

            return false;
          }
        } break;

        default:{
          std::cout << "Unknown control character \"" << control << "\" in DIMACS format." << std::endl;

          file.close();
          
          return false;
        }break;
      }
    }
  }else{
    std::cout << "Can't open file." << std::endl;

    file.clear();
    file.close();
    
    return false;
  }

  // remove first Vertex for compatible reasons
  remove_vertex(*(vertices(g).first),g);

  file.close();
  
  return true;
}
