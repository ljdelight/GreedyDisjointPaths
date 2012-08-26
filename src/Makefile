CC=g++
FLAGS=-g -Wall
PROJECT=GreedyDisjointPaths
FILES=Graph.h Graph.cpp Makefile greedyDisointPaths.cpp


$(PROJECT): Graph.o greedyDisjointPaths.o
	$(CC) $(FLAGS) $^ -o $@

greedyDisjointPaths.o: greedyDisjointPaths.cpp
	$(CC) $(FLAGS) -c $^
Graph.o: Graph.cpp Graph.h
	$(CC) $(FLAGS) -c $^

tar: $(PROJECT)
	tar -cvzf proj_lucas_burson.tar $(FILES)
clean:
	rm -f *~ *.o *.gch $(PROJECT) 

