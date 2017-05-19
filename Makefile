OBJ=main.o util_funcs.o SuffixArray.o BranchPointGroups.o Reads.o GenomeMapper.o string.o SamEntry.o
EXE=GeDi
CXX=g++
COMPFLAGS=-Wall -ggdb -MMD -pthread -std=c++11
OBJDIR=./objects/

$(EXE):$(OBJ)
	$(CXX) $(COMPFLAGS) $(OBJ) -o $(EXE) -lz -lboost_regex -lboost_program_options
#	mv *.o *.d ./obj

%.o: %.cpp
	$(CXX) $(COMPFLAGS) -c $<
-include $(OBJ:.o=.d)	

.PHONY: clean

clean:
	rm ./*.o
	rm ./*.d

cleaner:
	rm ./$(EXE)

