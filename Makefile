CXX = g++ -m64
CXXFLAGS = -O3 -Wall -g -std=c++11


all: output
output: *.cpp *.h
	$(CXX) $(CXXFLAGS) *.cpp *.h
