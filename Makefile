# Irrespective of the location of ./GeDi, automatically get the absolute path
# of the libdivsufsort dynamic library so GeDi can be run from any location.
libdivsufsort_lib = $(addprefix $(shell pwd), /libdivsufsort-master/build/lib/)
libdivsufsort_include = $(addprefix $(shell pwd), /libdivsufsort-master/build/include/)

OBJ=main.o util_funcs.o SNVIdentifier.o gsa.o GenomeMapper.o SamEntry.o CigarParser.o ssw.o ssw_cpp.o
LINKDIVSUF=-L $(libdivsufsort_lib) -ldivsufsort64 -Wl,-rpath $(libdivsufsort_lib)  
INCDIVSUF=-I $(libdivsufsort_include) -ldivsufsort64
EXE=GeDi
COMPFLAGS=-MMD -fopenmp -std=c++11 -O3 -DPWD=$(shell pwd)
CXX=g++

$(EXE):$(OBJ)
	$(CXX) $(COMPFLAGS) $(OBJ) -o $(EXE) $(LINKDIVSUF) -lz -lboost_regex -lboost_program_options

%.o: %.cpp
	$(CXX) $(COMPFLAGS) $(INCDIVSUF) -c $<
-include $(OBJ:.o=.d)	

.PHONY: clean

clean:
	rm ./*.o
	rm ./*.d

cleaner:
	rm ./$(EXE)
