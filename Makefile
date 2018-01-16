OBJ=main.o util_funcs.o SuffixArray.o SNVIdentifier.o Reads.o GenomeMapper.o SamEntry.o
LINKDIVSUF=-Wl,-R libdivsufsort-master/build/lib/  -L libdivsufsort-master/build/lib/ -ldivsufsort64
INCDIVSUF=-I libdivsufsort-master/build/include -ldivsufsort64
EXE=GeDi
CXX=g++
COMPFLAGS=-MMD -fopenmp -std=c++11 -O2 -rdynamic

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

