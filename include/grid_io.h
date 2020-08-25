#include<string>
#include<complex>
#include"assert.h"
#include<cstdio>
#include<vector>

namespace grid_io {
	static void printfunc(std::string FN,std::vector<double> omega,std::vector<std::complex<double>> F)
	{
		size_t N = omega.size();
		assert(N==F.size());
		
		FILE * file;
		file = fopen(FN.c_str(), "w");
		if (file == NULL) { 
			char buff[100]; 
			snprintf(buff, sizeof(buff), "-- ERROR -- Failed to open file %s : ",FN.c_str()); 
			perror(buff); };
		for (size_t i=0;i<N;i++){
			fprintf(file,"%f %f %f\n",omega[i],real(F[i]),imag(F[i]));
		}
		fclose(file);
	};

	static void printfunc(std::string FN,std::vector<double> X,std::vector<double> Y)
	{
		size_t N = X.size();
		assert(N==Y.size());
		
		FILE * file;
		file = fopen(FN.c_str(), "w");
		if (file == NULL) { 
			char buff[100]; 
			snprintf(buff, sizeof(buff), "-- ERROR -- Failed to open file %s : ",FN.c_str()); 
			perror(buff); };
		for (size_t i=0;i<N;i++){
			fprintf(file,"%f %f\n",X[i],Y[i]);
		}
		fclose(file);
	};
}
