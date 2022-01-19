default: App.o
	g++ App.o -o App libquasimont.a -lm -lquadmath -lgsl -lgslcblas
	./App

App.o: App.cpp
	g++  -c App.cpp -o App.o -ansi -std=c++17 -I../include
