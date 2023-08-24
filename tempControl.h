#pragma once
#ifndef _tempControl_h_
#define _tempControl_h_

#include "io.h"
#include <stdint.h>
#include "caHdl.h"
#include "grainCA.h"

#define MIN_TCELLS_PERSUBDOMAIN (10)
#define NO_MELT (12345)
#define LIQUID (6789)
#define TIMEDISCRETIZATION (16383)

//Corners of the temperature cell
#define c000 (0)
#define c001 (1)
#define c010 (2)
#define c011 (3)
#define c100 (4)
#define c101 (5)
#define c110 (6)
#define c111 (7)

typedef double Real;

class caHdl;//caHdl will check after increasing the time in which time interval the simulation is an increment only the counter
typedef caHdl* caHdlP;

typedef class temperatureNode TNode;
typedef temperatureNode* TNodeP;

typedef class tempNodePool TPool;
typedef TPool* tempNodePoolP;

typedef class tempCellClass tempCell;
typedef tempCell* tempCellP;

typedef class tempCellClass TCell;
typedef TCell* TCellP;

class tempCellPool;
typedef tempCellPool* tempCellPoolP;

/*typedef struct cellStruct cell;
typedef cell* cellP;
typedef struct cellPoolStruct cellPool;*/

/*struct nuclei;
typedef nuclei* nucleiP;*/

struct temperatureData
{
	TNodeP vertexTemp[8]; // pointer to temperatureNode class 
};

struct temperatureTimer 			// timer for temperature calculation
{
	long key;
	long size;
	Real* timeSeries;
};

struct TNodeRawData
{
	coordReal pos; // stores the x, y, z coordiantes of temperature node.
	long id; // Id of temp node.
};

struct cellNodesRawData
{
	long cellnodes[8]; // Stores the id of neighbouring 8 temperature Nodes.
};

class temperatureNode
{
public:
	temperatureNode(Real x, Real y, Real z, int identity); // constructor of temperatureNode class
	~temperatureNode(void);								   // destructor of temperatureNode class

	void allocateTemperatureTimeSeries(Real* Ttseries, long entries);
	Real interpolateTemperatureTime(Real prevTime, Real currentTime, Real nextTime, long pos, long nextPos);

	int id; // Id of the temperature node.
	coordReal pos; //x, y, z coordiante of the temperature node.
	Real * Temperature; // this array is used to store the temperature of the given temperature node for all the time steps.
};

class tempNodePool : public io
{
public:
	tempNodePool(tempCellPool* own = NULL);
	~tempNodePool(void);

	void readNodalTemperatureTimeSeries(const char* filename);

	TNodeP* all;

	tempCellPool* owner; 

	TNodeRawData* tnodesraw;
	Real* temperatureSeriesRawData;
	long size; // number of temperature nodes
	char* nodalTemperatureFilename;
};

class tempCellPool : public io
{
public:
	tempCellPool(caHdl* own = NULL);
	~tempCellPool(void);
	Real shortesttNucleationtime; //#as added to store the shortest nulceation time
	void setTimeEncodeFactor(void);
	long addTempCell(tempCellP Tc); // does nothing for now.
	void addTempCelltoBucket(tempCellP Tcell);
	Real getCellNucleationTime(uint_fast64_t data);
	Real getCellNucleationTime(cellP cinactive);
	void setTtimer(Real minTime);
	
	caHdl* owner;
	cellNodesRawData* tcellraw; 
	tempNodePoolP nodes; // pointer to tempNodePool class. tempNodePool class is owned by tempCellPool class.

	tempCellP first;
	tempCellP last;
	tempCellP* all;

	struct parallelNodeBucket 
	{
		tempCellP* TCparallelNode; // this is an array that stores the tempCellP of all the tenperature cell lying in my process.
		short size;
		short maxLen;

	} PNBucket;

	Real timeEncodeFactor; // What exactly is this..??

	Real xmin, xmax, ymin, ymax, zmin, zmax; // gives the min and man of the entire doamin.

	short xPerT, yPerT, zPerT;  //number of temperature cells in respective directions		
	long ntcells; // total number of temperature cells in the entire domain.

	Real* xdiscr;  //Temperature cells discretization
	Real* ydiscr;
	Real* zdiscr;	

	temperatureTimer T_timer;
};

class tempCellClass 
{
public:
	tempCellClass(tempCellPoolP own, TNodeP T000, TNodeP T001, TNodeP T010, TNodeP T011, TNodeP T100, TNodeP T101, TNodeP T110, TNodeP T111);
	tempCellClass(void);
	tempCellClass(tempCellPoolP own, short xmin, short xmax, short ymin, short ymax, short zmin, short zmax);
	~tempCellClass(void);

	tempCellPoolP owner(void) { return own; }
	caHdlP automaton(void) { return own->owner; }
	void allocateTemperatureTimeLocalNodes(int size);
	Real interpolateTemperatureSpace(short ix, short iy, short iz, Real currentTime, long timeKey, long nextKey);
	
	Real determineNucleationTime(short ix, short iy, short iz);
	Real TempParticle(short ix, short iy, short iz); // defined in tempControl.cpp
	tempCellPoolP own;
	temperatureData data;
	long initialCellNumber;
	TCellP next;
	TCellP prev;
	grainBoundingBox mySpace; // this is basically temperature bounding box. Directly copied from original code and never changed.
	short id;
	//cellP zombie;

};

#endif