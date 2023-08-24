#ifndef _polycrystal_h_
#define _polycrystal_h_

#include "world.h"
#include "mymath.h"
#include "grainCA.h"
#include "applic.h"

#define e_charge 1.602e-19
#define k_boltz 1.380662e-23

#define absMaxTemp (5273.15) //Units in Kelvin
#define TempFactor (absMaxTemp/(double)2097151) //The temperature will be encoded in 21 bits; the max possible value is 2097151

//growth machine velocity correction factors to let a polyeder grow closest to a sphere hull propagating with exactly the same normal velocity
//ideal island grain in monocrystal would also facet in twist and tilt boundaries thus in boundaries with most probable (slightly) intrinsic mobilities
//optimal: distorted polyeder with best fit to the sphere but minimal physical distortion of the node vertices
//how to get the compromise
//1. volume of convex polyeder hull closest to 1 --> (1 - Vpoly/Vsphere) < 0.002
//2. scaling constants as low as possible
//3. scaling constants as equal as possible
//compromise, because sphere volume is better approximated with larger and more anisotropic polyeder distortion
//compromise 26NN implementation 20130109 according to switch give quite a strong deviation // alternative // 26NNLoeckOriginal // 1 for face only

#define FGEOFACE 1 //1.0000 //1.0449		//1.0449 // 1 //1
#define FGEOEDGE 0.707106//0.7957 //1.1301 //0.8024 //1.1347 // 0.7957 //0.707106781186547
#define FGEODIAG 0.5773//0.6494 //1.1000 //0.5821 //1.0082 // 0.6494 //0.577350269189626 //###MK20121211, 26NN

#ifndef Real
typedef double Real;
#endif

typedef struct
{
	Real vMax;

	int output;   // #rg : not needed           
	int ElemNumb;      // #rg : not needed              
	int elemTyp;         // #rg : not needed            

	Real TwistAngleDeviation; // #rg : not needed
	Real AxisAngleDeviation; // #rg : not needed
	Real HighMobAngle; // #rg : not needed

        //rhos
	Real rhomax;
	Real rhoGND_av;
	Real rhoW_Average[25];     //0-9 ideal Oris + [11] overall rho average
	Real fgeo[27];

} physData;

typedef struct
{
   Real ASPIREDFILLPERSTEP;		// #rg : not needed
   Real MAXFILLPERSTEP;			
   char backgroundMode;
   int GGNumberOfDefGrainsLogged;	
   Real Tmelt;
   Real eta1;
   Real eta2;
   Real Uc;
   Real Tcrit_mean;
   Real Tcrit_sigma;
   Real delta ;  //as added ofr delta time
} genData;

typedef struct
{
	double mobil;
	double intern;
	double wall;
} allRho;

typedef class polyCrystal * polyCrystalP;

class polyCrystal : public world, public mathMethods
{
	friend class io;

public:
	polyCrystal( void * own );
	polyCrystal( void );
	~polyCrystal( void );

	caHdlP crystal( void ) { return (caHdlP) self; }
	Real GetMobilityWeight( grainP growGrain, grainP shrinkGrain);
	Real CalcMobilityWeight( grainP growGrain, grainP shrinkGrain );

public:
	Real time;
	Real deltaTime;
	Real timer;
	physData pd;
	genData gd;
	allRho aR;
	void * self;
};

#endif