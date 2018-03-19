OBJ=main.o util_funcs.o SNVIdentifier.o gsa.o GenomeMapper.o SamEntry.o CigarParser.o
LINKDIVSUF=-Wl,-R libdivsufsort-master/build/lib/  -L libdivsufsort-master/build/lib/ -ldivsufsort64
INCDIVSUF=-I libdivsufsort-master/build/include -ldivsufsort64
EXE=GeDi
CXX=g++
COMPFLAGS=-MMD -fopenmp -std=c++11 -O3 -rdynamic

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

