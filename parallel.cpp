#include "parallel.h"
#include "mymath.h"
#include "caHdl.h"

#ifndef CORe_MPI

//Insert here dummy functions for MPI

int MPI_Init( int * argc, char ** argv[])
{
	return 0;
}

int MPI_Comm_size( int dummy , int * nodes )
{
	*nodes = 1;
	return 1;
}

int MPI_Finalize( void )
{
	return 1;
}

int MPI_Comm_rank(int dummy, int * mynode )
{
	*mynode = 0;
	return 1;
}

int MPI_Barrier( int world )
{
	return 0;
}

int MPI_Bcast( void * buffer, int count, int type, int root, int comm ) //broadcast
{
	return 1;
}

int MPI_Sendrecv( void * buffer, int countsd, int typesd, int to, int tagsd, void * recv, int count, int typercv, int countrcv, int tagrcv, int comm, MPI_Status * st )
{
	return 1;
}

int MPI_Send( void * buffer, int count, int type, int to, int tag, int comm )
{
	return 1;
}

int MPI_Recv( void * buffer, int count, int type, int source, int tag, int comm, MPI_Status * st )
{
	return 1;
}

int MPI_Allreduce( void * in, void * out, int count, int type, int  op, int comm )
{
	if( type == MPI_CHAR )
		*(char *)out = *(char *)in;
	if( type == MPI_DOUBLE )
		*(Real *)out = *(Real *) in;
	if( type == MPI_INT )
		*(int *)out = *(int *) in;
	if( type == MPI_SHORT )
		*(short *)out = *(short *) in;
	return 1;
}

#endif

//prototype
parallelData::parallelData( long x0=0, long x1=0, long y0=0, long y1=0, long z0=0, long z1=0, int r=0, int xn=0, int yn=0, int zn=0 )
{
	minx = x0;
	maxx = x1;
	miny = y0;
	maxy = y1;
	minz = z0;
	maxz = z1;
	rank = r;
	ixn = xn;
	iyn = yn;
	izn = zn;
	memset( neighbourNodes, 0, sizeof( neighbourNodes ) );
}

parallelData::~parallelData( void )  //destructor for parallelData
{

}

processHdl::processHdl( caHdl * own=NULL, int totalNodes=1, long x=1, long y=1, long z=1, long xb=0, long yb=0, long zb=0 )
{
	owner = own; // caHdl is the owner.

	nodes = totalNodes; // Total number of processes.
	xPer = x; // number of processes in x- direction
	yPer = y;
	zPer = z;

	xbox = xb; // Number of automaton cells in the process in x direction.
	ybox = yb;
	zbox = zb;

	all = NULL;

	xnodes = (unsigned char *) calloc( xbox, sizeof(unsigned char) );
	ynodes = (unsigned char *) calloc( ybox, sizeof(unsigned char) );
	znodes = (unsigned char *) calloc( zbox, sizeof(unsigned char) );


	if( !xnodes || !ynodes || !znodes ) owner->err->reportError( ERRTXT( "Cannot allocate more memory" ) );

}

processHdl::~processHdl( void )
{
	if( all )
	{
		for( int i=0; i<nodes; i++ ) 
			delete all[i];

		free(all);
	}

	free( xnodes );
	free( ynodes );
	free( znodes );
}

int processHdl::getNextNode( long nix, long niy, long niz, long dx,  long dy, long dz )
{
	const char bper = owner->boxPer;
	const char nper = owner->nodePer;

	long xPern = xPer; //this is the number of nodes in each direction!
	long yPern = yPer;
	long zPern = zPer;

	long ix = nix + dx;
	long iy = niy + dy;
	long iz = niz + dz;
	//###20121228 how does this system behave if totalNodes > 1 holds? are the nodes then a closed 3D surface but not utilized
	//###20121228 or is nothing routed in them because protocols are explicitly asked during each potential send after updateFullRX?
	if( (bper & nper) == (BOXPER & NON_NODE_PER) || (bper & nper) == (NONBOXPER & NON_NODE_PER) )
	{
		if( iz < 0 )		iz += zPern;
		if( iz >= zPern )	iz -= zPern;
		if( iy < 0 )		iy += yPern;
		if( iy >= yPern )	iy -= yPern;
		if( ix < 0 )		ix += xPern;
		if( ix >= xPern )	ix -= xPern;

		return xPern * yPern * iz + xPern * iy + ix;
	}	
	return xPern * yPern * niz + xPern * niy + nix;
}