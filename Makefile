CXX = icpc

SP = src

prog = CPA

flags = -Wall -O3

includes = -Iinclude

all : $(SP)/CPA.cpp
	$(CXX) $(flags) -o $(prog) $(SP)/CPA.cpp $(includes)
	
	
# clean all files
clean :
	rm -vf prog

