#include"inline_math.h"
#include"grid_io.h"

#include<complex>
#include<vector>
#include<cstdio>
#include<string>
#include <argp.h>


const char *argp_program_version =
  "CPA-CT Aug-2020";
const char *argp_program_bug_address =
  "<henrytsang222@gmail.com>";
  
/* Program documentation. */
static char doc[] =
  "CT model solver using CPA";
  
/* INPUT DESCRIPTION */
static char args_doc[] = "[t1 t2 ...]";

/* Options */
static struct argp_option options[] = {
  {"verbose",  'v', 0,      0,  "Produce verbose output" },
  {"quiet",    'q', 0,      0,  "Don't produce any output" },
  {"silent",   's', 0,      OPTION_ALIAS },
  {"U",    'U', "Value",      0,  "Interaction Strength" },
  {"V",    'V', "Value",      0,  "c-f hopping" },
  {"mu",    'u', "Value",      0,  "chemical potential" },
  {"ef",    'e', "Value",      0,  "CT gap" },
  {"iterations",    'i', "COUNT",      0,  "Number of DMFT iterations" },
  {"omegamax",    'm', "Value",      0,  "Cutoff frequency" },
  { 0 }
};


/* Used by main to communicate with parse_opt. */
struct arguments
{
  char **args;
  int silent, verbose,max_iter,N;
  double U,V,mu,ef,omega_max;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = (struct arguments*)(state)->input;


  switch (key)
    {
    case 'q': case 's':
      arguments->silent = 1;
      break;
    case 'v':
      arguments->verbose = 1;
      break;
    case 'U':
      arguments->U = arg ? atof(arg) : 2.0;
      break;
    case 'V':
      arguments->V = arg ? atof(arg) : 0.217;
      break;
    case 'u':
      arguments->mu = arg ? atof(arg) : 0.0;
      break;
    case 'e':
      arguments->ef = arg ? atof(arg) : -1.0;
      break;
    case 'i':
      arguments->max_iter = arg ? atoi(arg) : 10;
      break;
    case 'N':
      arguments->N = arg ? atoi(arg) : 1000;
      break;
    case 'm':
      arguments->omega_max = arg ? atof(arg) : 6.0;
      break;
		
    case ARGP_KEY_ARG:

      arguments->args[state->arg_num] = arg;

    break;

    case ARGP_KEY_END:
      if (state->arg_num < 1)
        /* Not enough arguments. */
      argp_usage (state);
      break;
		
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };


int main(int argc, char* argv[])
{
  struct arguments arguments;

  /* Default values. */
  arguments.silent = 0;
  arguments.verbose = 0;
  arguments.U = 2.0;
  arguments.V = 0.2;
  arguments.max_iter = 100;
  arguments.N = 2000;
  arguments.omega_max = 6.0;

	size_t count=0;
	while(argv[++count]);
	char ** args = new char * [count];
	arguments.args = args;
	
  argp_parse (&argp, argc, argv, 0, 0, &arguments);
  
	using namespace std;
	
	int verbose = arguments.verbose;
	int quiet = arguments.silent;
	
	size_t N = arguments.N;
	size_t maxit = arguments.max_iter;
	double U=arguments.U;
	
	double V=arguments.V;
	double mu=arguments.mu;
	double ef = arguments.ef;
	double ef1 = ef;
	double ef2 = ef1+U;
	
	if (!quiet) printf("max_iter=%lu N=%lu omega_max=%.02f\nU=%.02f V=%.02f mu=%.02f ef=%.02f\n",maxit,N,arguments.omega_max,U,V,mu,ef);
	
	//Read in number of t
	vector<double> ts;
	for (int i=0;i<count-1;++i) if(arguments.args[i]!=NULL) ts.push_back(atof(arguments.args[i]));
	
	for (auto & t:ts) {printf("Work on ");printf("%.02f ",t);printf("\n");}
	
	for (auto & t:ts){
	
		size_t iter = 0;
		
		if (!quiet) printf("Start CPA for t = %f\n",t);
		
		//Init grid
		double omega_max = arguments.omega_max;
		double domega = 2.0*omega_max/((double) N-1.0);
		vector<double> omega;
		for (int i=0;i<N;i++) omega.push_back(-omega_max+domega*i);
		
		//Initial Gc
		vector<complex<double>> Gc;
		for (int i=0;i<N;i++) Gc.push_back(complex<double>(0.0,-1.0));
		
		//Compute Deltac
		vector<complex<double>> Deltac;
		for (int i=0;i<N;i++) Deltac.push_back(sqr(t)*Gc[i]);
		
		while (iter<maxit){
		
			if (!quiet and verbose) printf("CPA: iter %d\n",iter);
			
			//Compute Deltaf
			vector<complex<double>> Deltaf;
			for (int i=0;i<N;i++) Deltaf.push_back(sqr(V)/(omega[i]+mu-Deltac[i]));
			 
			//Compute Gf from Gc
			vector<complex<double>> Gf;
			for (int i=0;i<N;i++){
				const complex<double> Gf1 = 1/(omega[i]+mu-ef1-Deltaf[i]);
				const complex<double> Gf2 = 1/(omega[i]+mu-ef2-Deltaf[i]);
				Gf.push_back(0.5*Gf1+0.5*Gf2);
			}
			
			//Now compute self energy Sigma
			vector<complex<double>> Sigma;
			for (int i=0;i<N;i++) Sigma.push_back(omega[i]+mu-ef-Deltaf[i]-1/Gf[i]);
			
			//Compute Phi
			vector<complex<double>> Phi;
			for (int i=0;i<N;i++) Phi.push_back(sqr(V)/(omega[i]+mu-Sigma[i]));
			
			//Compute new Gc
			for (int i=0;i<N;i++) Gc[i] = 1/(omega[i]+mu-Deltac[i]-Phi[i]);
			
			//Mixing
			for (int i=0;i<N;i++) Deltac[i] = 0.5*sqr(t)*Gc[i]+0.5*Deltac[i];
			
			
			if (iter==maxit-1){//Final iteration
				if (!quiet and verbose) printf("printing results to file\n");
				//Compute dos
				vector<double> dos;
				for (int i=0;i<N;i++) dos.push_back(-imag(Gf[i]+Gc[i])/M_PI);
				
				char buffer [100];
				snprintf ( buffer, 100, "V.%.02f.U.%.02f.t.%.02f",V,U, t );
				
				grid_io::printfunc("Sig.out."+string(buffer),omega,Sigma);
				grid_io::printfunc("Gf.out."+string(buffer),omega,Gf);
				grid_io::printfunc("Gc.out."+string(buffer),omega,Gc);
				grid_io::printfunc("dos.out."+string(buffer),omega,dos);
			}
			iter++;
		}
		
		if (!quiet) printf("End CPA\n");
	}
	
	delete [] args;
	return 0;
}
