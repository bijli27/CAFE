#pragma once
#ifndef _parallel_h_

#define _parallel_h_

#include <stdlib.h>
#include "applic.h"

#define MASTER (0)
#define EVEN false
#define ODD true

#define SEND 98
#define RECV 98
#define SENDRECV 99
#define STATUSSEND 100

class caHdl;

#define CORe_MPI

#ifdef CORe_MPI
#include <mpi.h>
#else
//insert here dummy MPI functions

#define MPI_COMM_WORLD 1

#define MPI_CHAR 0
#define MPI_INT 1
#define MPI_DOUBLE 2
#define MPI_SHORT 3

#define MPI_MAX 1
#define MPI_BOR 1
#define MPI_SUCCESS 0

typedef int MPI_Status;

int MPI_Init( int * argc, char ** argv[]);
int MPI_Comm_size( int dummy , int * nodes );
int MPI_Comm_rank(int dummy, int * mynode );
int MPI_Finalize( void );
int MPI_Barrier( int world );
int MPI_Bcast( void * buffer, int count, int type, int root, int comm );
int MPI_Sendrecv( void * buffer, int countsd, int typesd, int to, int tagsd, void * recv, int count, int typercv, int countrcv, int tagrcv, int comm, MPI_Status * st );
int MPI_Send( void * buffer, int count, int type, int source, int tag, int comm );
int MPI_Recv( void * buffer, int count, int type, int source, int tag, int comm, MPI_Status * st );
int MPI_Allreduce( void * in, void * out, int count, int type, int  op, int comm );

#endif

typedef class processHdl * processHdlP; 
typedef class parallelData * parallelDataP;

class parallelData   // defining parallel data class
{
public:
	parallelData( long x0, long x1, long y0, long y1, long z0, long z1, int r, int xn, int yn, int zn ); 
	~parallelData( void ); //destructor - it will be invoked when scope of the object gets over

	inline int myNode( void ) { return rank; }

	int rank;
	long minx;
	long maxx;
	long miny;
	long maxy;
	long minz;
	long maxz;
	int ixn;
	int iyn;
	int izn;
	int neighbourNodes[6]; //###20121228done was 26
	//###20121228done
};


class processHdl   // classdefinition processHdl
{
public:
	processHdl( caHdl * own, int totalNodes, long xPer, long yPer, long zPer, long xbox, long ybox, long zbox );	//constructor
	~processHdl( void );		//destructor - it will be invoked when scope of the object gets over

	parallelDataP * all;		//pointer to parallelDataP
	caHdl * owner;				//pointer to caHdl

	int cellNode( long ix, long iy, long iz ) { return xPer * yPer * znodes[iz] + xPer * ynodes[iy] + xnodes[ix]; }
	int getNextNode( long nix, long niy, long niz, long dx,  long dy, long dz );

	unsigned char * xnodes; // array whose index is cell x-cordiante and value is process number along x direction.
	unsigned char * ynodes; // array whose index is cell y-cordiante and value is process number along y direction.
	unsigned char * znodes; // array whose index is cell z-cordiante and value is process number along z direction.

	int nodes; // total number of processes.

	long xPer;  //this is the number of nodes in each direction!
	long yPer;
	long zPer;

	long xbox; // Number of automaton cells in the process in x direction.
	long ybox; // Number of automaton cells in the process in y direction.
	long zbox; // Number of automaton cells in the process in z direction.
};

#endif