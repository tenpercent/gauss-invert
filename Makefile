CC = g++
CFLAGS = -Wall -O3
SOURCES = gauss.cpp multiplication.cpp simple_inversion.cpp tools_extended_cycles.cpp main-rewrite.cpp

all:	main

main:
	$(CC) $(CFLAGS) $(SOURCES) -o solve.out
clean:
	rm -v *.out
