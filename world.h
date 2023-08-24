#pragma once
#ifndef _world_h_
#define _world_h_

#include "interface.h"
#include "applic.h"

#ifdef CORe_MPI
#include <mpi.h>
#endif

#define BOXPER			(0x07)
#define NONBOXPER		(0x0F)  //same identifiers? (0x0F)

#define PLOTRGB			(0x01)
#define PLOTINFPROGRESS (0x02)

#define NODE_PER		(0x0E)
#define NON_NODE_PER	(0x0F) //same identifiers?

#define NOCOMM			(0x09)
#define DETERMINISTIC	(0x0D)
#define STATISTICAL		(0x0B)
#define COMM_ENABLED	(0x06)

#define STEPS			(0x01)

#define	INTERFACE		(0x01)
#define BULK			(0x02)

//###current nucleation model
#define NODF_USINGLE	(0x01)
#define NODF_UMULT		(0x02)
#define NODF_UPRNG		(0x03)

#define ONEPROCESS (1)

//###current nuclei placement model
#define NSITE_UHOMOFULL			(0x11)
#define NSITE_UHOMOEIGHT		(0x12)
#define NSITE_UCENTRE			(0x13)
#define NSITE_UYZFULL			(0x14)
#define NSITE_UYZLINE			(0x15)
#define NSITE_UPLACEEXPLICITLY	(0x16)
#define NSITE_UCLUSTERING		(0x17)
#define NSITE_UCLUSTERINGSINSQ	(0x18)
#define NSITE_UREHER			(0x19)


class world
{
public:
	world( void ); 			//constructor 
	~world( void );			//destructor

	virtual int sendToInterface( long ix, long iy, long iz, int mynode, int nextNode );
	virtual void distributeReceivedCells( void );
	int myNode( void ) { return rank; }
	virtual bool isMyRankInTheList(short* list, int size);		//checking rank is in the list or not

	virtual void sendReceive( void );

	virtual void receiveSend( int sender );
	void receiveSend( void );
	
	int neighbourNode( int who ) { return nodeFace[who]->nextNode; }
	
public:

	Real xsize; // length in x direction of total domain
	Real ysize; // length in y direction of total domain
	Real zsize; // length in z direction of total domain

	char boxPer;
	char nodePer;
	char protocol;
	char commType;
	char insertionSite;

	char plotmode;//rg: not needed.

	char nucodfmethod;//rg: not needed.
	char nucsitemethod;//rg: not needed.

	//###20121228done

public:
	// these set of values are set for each processes differently.
	long x0;					//discretized starting offset where the subdomain is can be abstracted as a this0 vector
	long y0;					
	long z0;

	long xPer;                  //number of cells in direction x within the subdomain automat
	long yPer;					//number of cells in direction y within the subdomain automat
	long zPer;					//number of cells in direction z within the subdomain automat
	long zFac;

	int rank;
	   
	worldInterfaceP nodeFace[6]; //###20121228done was 26
	char somethingToSend;			//###20121228done was 26
	char somethingReceived; 	//###20121228done was 26
};

#endif