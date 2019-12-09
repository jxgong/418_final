CXX = g++ -m64
CXXFLAGS = -O3 -Wall -g -std=c++11
NVCC = nvcc
NVCCFLAGS=-O3 -m64 --gpu-architecture compute_35

OBJDIR=objs
LIBS       :=
FRAMEWORKS :=

LDFLAGS=-L/usr/local/depot/cuda-8.0/lib64/ -lcudart
# LIBS += GL glut cudart
OBJS = $(OBJDIR)/main.o  $(OBJDIR)/cudaSimulator.o $(OBJDIR)/common.o\
	   $(OBJDIR)/fileio.o

LDLIBS  := $(addprefix -l, $(LIBS))
LDFRAMEWORKS := $(addprefix -framework , $(FRAMEWORKS))
all: output

dirs:
	mkdir -p $(OBJDIR)/

output: dirs $(OBJS) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS) $(LDFRAMEWORKS)

$(OBJDIR)/%.o: %.cu cudaSimulator.h
	$(NVCC) $< $(NVCCFLAGS) -c -o $@

$(OBJDIR)/%.o: %.cpp physics.h
	$(CXX) $< $(CXXFLAGS) -c -o $@

clean:
	rm -rf $(OBJDIR)
	rm output
