#ifndef _random_h_
#define _random_h_

#include <math.h>	//for math functions
#include <time.h>	//for time()

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32										
#define NDIV (1+(IM-1)/NTAB)						//NDIV = 67108864
#define EPS 1.2e-30								    //for randomInterval
#define RNMX (1.0-EPS)								//for randomInterval()		
#define ISNAN(a) ((a) != (a))						//check for NaN
#define _PI_ 3.1415926535897932384626433832795		//PIE 
#define EPSILON 0.01				//epsilon for floating point comparison
#define EPS1 0.001					//epsilon for comparing two numbers
#define EPS2 1.0e-8					//for numerical stability
#define SQR(a) ((a)*(a))			//for using sqrt()
#define CUBE(a) ((a)*(a)*(a))		//for using cube()
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))			//accepts two values and return smaller
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))			//accepts two values and return bigger
#define MYHASH(a,b) ( (((uint_fast64_t) max(a,b)) << 32) | ((uint_fast64_t) min(a,b)) ) 
#define SIGN(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)	//to check +ve & -ve sign
#define XOR(a,b) ((!!a) != (!!b))					//takes two argument as operand & does XOR on everybit of two numbers

class randomClass
{
public:
	randomClass(long seed = 0);

	void init(long sd) { seed = sd; }
	//void init(long sd) { seed = (long) time(NULL);} hat nicht funktioniert

	long integer(void);
	double real(void);
/*	double gauss(void);
	double wls(void);*/
	double parkMiller( void );
	double randomInterval( double rmin, double rmax );

  private:
		long seed;
};

typedef randomClass *randomClassP;

#endif