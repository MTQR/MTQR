CXX = g++

CXXFLAGS = -g -Wall
DLIBS = -lm -lgsl -lgslcblas -lgmp -lquadmath

DEPS = ../../include
TEST1 = ../../output/test1
TEST2 = ../../output/test2
TEST3 = ../../output/test3

vpath %.h $(DEPS)

main: *.cpp $(DEPS)/*.h $(TEST1)/*.h
	g++ $(CXXFLAGS) -I $(DEPS) -I $(TEST1) *.cpp -o main $(DLIBS)

all: test1 test2 test3

test1: *.cpp $(DEPS)/*.h $(TEST1)/*.h
	g++ $(CXXFLAGS) -I $(DEPS) -I $(TEST1) *.cpp -o main $(DLIBS)
	./main test1

test2: *.cpp $(DEPS)/*.h $(TEST2)/*.h
	g++ $(CXXFLAGS) -I $(DEPS) -I $(TEST2) *.cpp -o main $(DLIBS)
	./main test2

test3: *.cpp $(DEPS)/*.h $(TEST3)/*.h
	g++ $(CXXFLAGS) -I $(DEPS) -I $(TEST3) *.cpp -o main $(DLIBS)
	./main test3

clean:
	rm -f main