#pragma once
#ifndef _grainCA_h_
#define _grainCA_h_

#include "io.h"
#include <stdint.h>

#define MINMIS 0.017 //0.0035

#define HASHCOLUMNMASK 0x000000f0
#define HASHCOLUMNSHIFT 4
#define HASHCOLUMNCOUNT ((HASHCOLUMNMASK>>HASHCOLUMNSHIFT)+1)
#define MAXLENGTH 40000
#define ASPIREDFILLPERSTEP_DEFAULT 0.2	//0.01  //default 0.2
#define MAXFILLPERSTEP_DEFAULT 1.0		//0.05  //default 1.0
#define TIMESLOTCOUNT_DEFAULT 1000		//Mohles work
#define TIMESLOTCOUNT 1000 //default 1000
#define NOACTIONOUTCOUNT 100000
#define _RHOFACUNIT 50
#define RHOFACUNIT (1.0/_RHOFACUNIT)

#define _MOBFACUNIT 0xFFFF
#define MOBFACUNIT (1.0/_MOBFACUNIT)

#define STDLEN 16

typedef double Real;

typedef class grainClass grain;
typedef grain* grainP;

typedef struct cellStruct cell;
typedef cell *cellP;

typedef struct cellPoolStruct cellPool;

class boundariesPool;
typedef boundariesPool * boundariesPoolP;

class cellsBoundary;
typedef cellsBoundary * cellsBoundaryP;

class caHdl;
typedef caHdl * caHdlP;

struct nuclei;
typedef nuclei * nucleiP;


inline long cellExists(cellP c) { return !((long)c & 0x00000001); }
inline long PorosityExists(cellP pd) { return !((long)pd & 0x00000001); } // #as:added for prorsoity check
inline long isCellSolid(cellP c) { return !((long)c & 0x00000006); }

//inline float getRhoFacOfCP(cellP c)		{ return RHOFACUNIT * (unsigned char)((long)(c) >> 8); }

inline float getRhoFacOfCP(cellP c)		{ return RHOFACUNIT * (unsigned char)((long)(c) >> 3); }
inline long grainArrayPosition(cellP c) { return ((long)0); }  //Because there is only one phase where the grains nucleate
//inline long grainArrayPosition(cellP c) { return ((long) (c) >> 11); }
inline long cellTemperature(cellP c) { return (((long)(c)) >> 11); }


struct coord
{
	short ix;
	short iy;
	short iz;
};

struct coordReal
{
	Real x;
	Real y;
	Real z;
};

struct gpData // grain pool data
{
	long oriCount;
	long count;                 //Number of grains

	long xPer;					//Number of grains in x direction
	long yPer;
	long zPer;
	long zFac;					//Number of grains in the xy plane - safe variable dereferencing

	long xCellCount;			//discrete grain size in cells is at least one
	long yCellCount;
	long zCellCount;
};

struct cpData // cell pool data
{
	long xPer;		//domain size in cells x >= y >= z
	long yPer;
	long zPer;
	long zFac;

	char syncRequired; // not needed.

	long count;                  //number of cells
	long retiredCount;

	Real CountCrit;             //MGS
	Real size;
	Real Volume;
	Real totalRXfrac;
};

struct cellPoolStruct
{
	cellP first[TIMESLOTCOUNT+1]; // this is the array used to strore the memory address of all the active cells.it is assumed that maximum number of active cells is equal to TIMESLOTCOUNT+1
	cellP next;
	cellP closeCell;

	cellP *all;

	cellP firstRecycled;
	long countRecycled;
	cpData data;
};

struct cellStruct		////////////// this struct should be kept small !! LB: This new cell struct is "only" 3 bytes larger than the original one but contains much more info
{
    short iz;// this is the global position of cell
	short iy;
	short ix;
	unsigned char rhoFac_RFU;  
	cellP next;
	Real rxFrac;

	Real T; //Let's save the temperature here for the moment.

	uint_fast64_t dataHolder; //Dataholder contains (right to left): infector (5 bits), P (16 bits), sign of P (1 bit), rx grain (21 bits), deformed grain (21 bits) = 64 bits

	//                                              sign of P|
	//                                 1: negative 0 positive|
	// |xxxx xxxx xxxx xxxx xxxx x|xxx xxxx xxxx xxxx xxxx xx|x|x xxxx xxxx xxxx xxx|x xxxx
	// |                          |                          | |                    |
	// |deformed grain            |rx grain                  | |P                   |infector


	//SOL: Deformed grain is not needed anymore because there is only one phase (the liquid). 21 bits allow for 2097151 different integers
	//that can be used to encode the temperature.
};

struct grainBoundingBox
{
	short xleft;
	short xright;
	short yleft;
	short yright;
	short zleft;
	short zright;
};

typedef struct{ //LB: For the microstructure generator use only the struct
	Real phi1;
	Real PHI;
	Real phi2;
	Real q0;
	Real q1;
	Real q2;
	Real q3;
	long count;
} exportOri;

class oriRepresentation
{
public:
	oriRepresentation( void );
	~oriRepresentation( void );

	Real phi1;
	Real PHI;
	Real phi2;
	Real q0;
	Real q1;
	Real q2;
	Real q3;
	long count;
	unsigned char RGBA[4];
	//Insert functions for ori operations

	void copyOriFrom(oriRepresentation* original);
};
typedef oriRepresentation * oriRepresentationP;

class grainPool : public io
{
public:
	grainPool( caHdl * own=NULL );
	~grainPool( void );

	caHdl * owner;

	int addOrientation( Real phi1, Real PHI, Real phi2 );
	long addGrain( grainP g );

	grainP first;
	grainP last;

	grainP *all;
	
	gpData data; 
    long Nuc; // this stores the total number of nucleus in the entire complete domain, not just one process.
	long defGrainsCount; // number of initial deformed grains. for AM simulation this remains 1 which signifies powder domain of orientation(-1,-1,-1)
	long NucMG;
	Real commonGrainDiameter;
	Real subGrainDiameter;
	oriRepresentationP * orientations;
	long noris;
	long ngrains;
};
typedef grainPool * grainPoolP;

class cellsBoundary //###LB: struct is better, owner and disori not needed
{
public:
	cellsBoundary( boundariesPoolP own = NULL, cellP up = NULL, cellP down = NULL );
	~cellsBoundary( void );

	coord upGrain;
	coord downGrain;

	cellsBoundaryP prev;
	cellsBoundaryP next;
	boundariesPool * owner;
	Real disori;
};

typedef struct
{
	long count;
	uint_fast32_t gposUp;
	uint_fast32_t gposDown;
	uint_fast64_t id;
	Real disori;
	cellsBoundaryP first;
	cellsBoundaryP last;
	long tempCount;
}boundaryEntry;
typedef boundaryEntry * boundaryEntryP;

typedef struct
{
	long len;
	long maxLen;
	boundaryEntry entry[];
}boundaryColumn;
typedef boundaryColumn * boundaryColumnP;

typedef struct
{
	boundariesPoolP owner;
	cellsBoundaryP first;
	boundaryColumnP gb;
}boundary;

typedef struct
{
    grainP identity;
	Real P;   // weight factors for mobility calculation

} hashTabEntry;

typedef struct
{
	long len;
    long maxLen;
    hashTabEntry entry[];
} hashColumn;
typedef hashColumn *hashColumnP;

class boundariesPool : public io
{
public:
	boundariesPool( caHdl * own = NULL ) : 
		owner(own), first(NULL), last(NULL), count(0), gbCount(0), offset(0)
		{
			memset( MODFCounts, 0, sizeof( MODFCounts ) );
			memset( MODF, 0, sizeof( MODF ) );
			memset( userMODF, 0, sizeof( userMODF ) );
			memset( userMODF_Cumulative, 0, sizeof( userMODF_Cumulative ) );
		}

	~boundariesPool( void );
	void emptyBoundaries( void );
	caHdl * owner;

	boundaryColumnP * boundaryBucket;
	boundaryEntryP * boundaryEntries;

	cellsBoundaryP first;
	cellsBoundaryP last;
	long count;
	long gbCount;
	long bucketSize;
	long offset;

	long MODFCounts[32];
	Real MODF[32];
	Real MODF_Cumulative[32];

	Real userMODF[32];
	Real userMODF_Cumulative[32];

	void calculateMODF( void );
	void calculateFastMODF( void );
	void initializeBoundaries( void );
	void addBoundaryPiece( long pUp, long pDown, Real disori );
	void recalculateDisorientation( long posga, long posgb );
	Real calculateDisorientation( grainP up, grainP down );
	void sortGrainBoundaryCells( void );
};

struct nuclei
{
	short ix;
	short iy;
	short iz;
	long grainIdx;
	nucleiP next;
	nucleiP prev;
};

struct cellsToNucleate
{
	nucleiP first;
	nucleiP last;
	long count;
	char nucleateFlag;
};

typedef struct 
{
	short x;
	short y;
	short z;
	long ig;
} exportNucleus;

typedef struct
{
	long xl;              //outermost left cell of particle in x
	long yl;              //outermost left cell of particle in y
	long zl;              //outermost left cell of particle in z
	long xr;             //outermost left cell of particle in x
    long yr;             //outermost left cell of particle in y
    long zr;             //outermost left cell of particle in z
    long count;          //counts no of cells per particle
	long pos;
} exportParticle;

struct nucleiPacket
{
	long count;
	exportNucleus ** nucList;
};

class grainClass
{
public:
	grainClass( grainPoolP own,  Real p1, Real P, Real p2, double rhoGND, double rm0, double ri0, double rw0,int n1,int n2,int n3, short Oriindex, Real NGLS, long DeformGrainIndex );
	grainClass(  grainPoolP own, Real p1, Real P, Real p2, int index1, Real scatter );
	grainClass( void );
	~grainClass( void );

	caHdlP automaton( void ) { return owner->owner; }

	grainPoolP owner;
              
	Real RVFac3;                 /*1/NGLS*/
	double rhoGND;               /*geometrical necessary dislocations*/
	double rhoMobil0;            /*mobile disloc.dens. initial*/
	double rhoIntern0;           /*interior disloc.dens. initial*/
	double rhoWall0;             /*wall disloc.dens. initial*/
    double rhoM;                 /*mobile disloc.dens. time-updated*/
	double rhoI;                 /*interior disloc.dens. time-updated*/
	double rhoW;                 /*wall disloc.dens. time-updated*/
	double NGLS;                 /*number of active slip systems*/

	double rhoTotal;

	int nuc1, nuc2, nuc3;        /*saves Nuc-Modus 1: SB 2:GB 3: TB*/ 

	long alreadyNuc;             /*flag to distinguish if or if not nucleated*/
	short Oriindex;              /*Index represents the affiliation to an ideal ori*/
	long DeformGrainIndex;       /*Index of grains correspondig to grain input list*/
	Real scatter;               /*Oriscatter*/

	long CellCount;              /*counts no of 100% Rx cells, respectively Particlecells per grain*/
	long initialCellNumber;

	grainP next;
	grainP prev;

	cellP zombie;

	int colR;
	int colG;
	int colB;

	hashColumnP cp[HASHCOLUMNCOUNT];
	int oriIndex;
	long arrayPos;
	grainBoundingBox mySpace;

	cellsToNucleate GBnuclei;					//list of nuclei to be set in the grain
	cellsToNucleate TBnuclei;
	cellsToNucleate SBnuclei;
	cellsToNucleate PSNnuclei;

	char somethingToNucleate;
};

// #as: for storing the data of porosity
struct porosity {
	vector<Real> x;
	vector<Real> y;
	vector<Real> z;
	Real* total;
	vector<Real> tmax;

};



#endif
 