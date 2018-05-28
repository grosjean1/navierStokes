#  ubuntu 
UMFPACKINC = -I/home/elise/Documents/M2/Projet_NavierStokes_EliseGrosjean/SuiteSparse/UMFPACK/include
UMFPACKLIBS = -L/home/elise/Documents/M2/Projet_NavierStokes_EliseGrosjean/SuiteSparse/UMFPACK/Lib -lumfpack -lcholmod -lccolamd -lcolamd -lcamd -lamd -lsuitesparseconfig -lblas

CXXCHECK = -g -pg  -fno-optimize-sibling-calls -O0
CXXOPT =  -O3

CXXFLAGS =  $(CXXCHECK) -Wall -std=c++11 $(UMFPACKINC)
CXXFLAGS += -MMD -MP
PROGS =  NS
OBJS  = mesh.o mainNS.o
SRC = mesh.cpp mainNS.cpp
all: $(PROGS)

-include $(SRC:%.cpp=%.d)
	
NS: $(OBJS)
	$(CXX) -o $@ $^  $(CXXFLAGS) $(UMFPACKLIBS)

clean: 
	-rm $(PROGS) *.o *~  *.txt *.exe *.d 
