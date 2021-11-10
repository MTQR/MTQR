default: App.o
	g++ App.o -o App path/to/QUASIMONT/build/libquasimont.a -lm -lgsl -lgslcblas -lgmp -lquadmath

App.o: App.cpp
	g++  -c App.cpp -o App.o -g -ansi -std=c++17 -Ipath/to/QUASIMONT/include
