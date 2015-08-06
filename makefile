C = gcc
CXX = g++
CXXFLAGS = -std=c++11 -Wall
DEBUG = -g
OPT = -O3

eqDsatur : EqColoring.o Coloring.o main.cpp
	$(CXX) $(CXXFLAGS) -o eqDsatur main.cpp EqColoring.o Coloring.o

EqColoring.o : EqColoring.cpp EqColoring.hpp Coloring.hpp
	$(CXX) $(CXXFLAGS) -c EqColoring.cpp -o EqColoring.o

Coloring.o : Coloring.cpp Coloring.hpp
	$(CXX) $(CXXFLAGS) -c Coloring.cpp -o Coloring.o

clean:
	rm eqDsatur *.o

tar :
	tar -cvf eqDsatur.tar *.hpp *.cpp makefile
