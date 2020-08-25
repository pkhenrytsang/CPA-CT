CXX = icpc

SP = src

prog = CPA

includes = -Iinclude

all : $(SP)/CPA.cpp
	$(CXX) -o $(prog) $(SP)/CPA.cpp $(includes)
	
	
# clean all files
clean :
	rm -vf prog

