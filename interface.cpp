
#include "interface.h"
#include <string.h>
#include "world.h"
#include "grainCA.h"
#include "caHdl.h"


oriBucket::oriBucket( worldInterfaceP own ) :
	owner(own)
{
	for( int i = 0; i < BUCKETSIZE; i++ )	oriContainer[i] = NULL;
	allOri = 0;
}
	
oriBucket::~oriBucket( void )
{
	for( int i = 0; i < BUCKETSIZE; i++ )	free( oriContainer[i] );
}

void oriBucket::addOri( long oriIndex, long maxOri )  //maxOri is an estimator how one expects the oriBin to be###
{
	int pos = 0;

	if( maxOri < BUCKETSIZE )
		pos = oriIndex;
	else
		pos = (int) ( ( BUCKETSIZE - 1 ) * ( (Real) oriIndex / maxOri ) );

	oriBinP bin = oriContainer[pos];

	if( !bin ) //create the bin if it not exists
	{
		bin = (oriBinP) malloc( sizeof(oriBin) + STDLEN * sizeof(oriCounter) );
		if( !bin ) owner->err->reportError( ERRTXT( "Cannot allocate more memory" ) );

		bin->len = 0;
		bin->maxLen = STDLEN;
		oriContainer[pos] = bin;
	}

	int i;
	for ( i = 0; i<bin->len; i++ )
	{
		if( bin->entry[i].oriIndex == oriIndex ) //if the current oriIndex is the one in the oribucket, it is increased
		{
			bin->entry[i].count++;
			return;
		}
	}

	//and account for this also for the whole OriBucket, this might make a re-malloc necessary
	if( i < bin->maxLen )
	{                                                                       
		bin->entry[i].oriIndex = oriIndex;                                            
		bin->entry[i].count = 1;
		bin->len++;
		allOri++;
		return;
	}

	//*re-malloc procedure
	if( bin->len != bin->maxLen ) owner->err->reportError( ERRTXT("bin->len != bin->maxLen") );

	oriBinP binlong = (oriBinP) malloc( sizeof(oriBin) + bin->maxLen * 2 * sizeof(oriCounter) );
	if( !binlong ) owner->err->reportError( ERRTXT( "Cannot allocate more memory" ) );

	binlong->len = bin->maxLen;
	binlong->maxLen = 2 * bin->maxLen;

	for( i=0; i<bin->maxLen; i++)
	{
		binlong->entry[i] = bin->entry[i];
	}

	binlong->entry[bin->maxLen].oriIndex = oriIndex;
	binlong->entry[bin->maxLen].count = 1;
	binlong->len++;
	allOri++;

	oriContainer[pos] = binlong;
	free(bin);
}

long * oriBucket::countOrientations(long * size)
{
	int i;
	int k = 0;       //temporary counter
	int counter = 0; //counts how many orientation in total in this bucket

	long * temp = (long *)calloc(2 * allOri, sizeof(long));

	for (i = 0; i < BUCKETSIZE; i++)
	{
		if (!oriContainer[i]) continue;

		oriBinP bin = oriContainer[i];

		int length = bin->len;

		for (int j = 0; j < length; j++)
		{
			temp[k] = bin->entry[j].oriIndex;
			temp[k + 1] = bin->entry[j].count;
			k += 2;
			counter += bin->entry[j].count;
		}

	}
	*size = 2 * allOri;

	QUICKASSERT(k / 2 == allOri);

	emptyBucket();

	return temp;
}

void oriBucket::emptyBucket(void)
{
	for (int i = 0; i < BUCKETSIZE; i++)
	{
		free(oriContainer[i]);
		oriContainer[i] = NULL;
	}
	allOri = 0;
}

worldInterface::worldInterface(caHdlP own, long ix, long iy, short me) :
	owner(own), aPer(ix), bPer(iy), self(me)
{
	count = ix * iy;

	recvCells = NULL;

	interCells = (long *)calloc(count, sizeof(long));

	if (!interCells) err->reportError(ERRTXT("Cannot allocate more memory"));

	sendRecv = 0;

	somethingToSend = NOTHING;
	somethingReceived = NOTHING;
	bucket = NULL;

	if (owner)
	{
		err = owner->err;
		protocol = owner->protocol;
		if (protocol == STATISTICAL) bucket = new oriBucket(this);
	}
	else
		err = NULL;

	active = 0;
	retired = 0;
	recvSize = 0;

}

worldInterface::~worldInterface(void)		//destructor
{
	free( interCells );
	delete bucket;
}

char worldInterface::addCell100(cellP infector, int a, int b)
{
	long a0 = 0; //a and b are global addresses in the whole CA domain!
	long b0 = 0; //a0 and b0 are global addresses where the interface is modeled sectioned from the domain
	long aPer = 0;
	long bPer = 0;
	//###why are z0 and y0 switched in their order
	if (self == C_LEFT || self == C_RIGHT) { a0 = owner->z0; b0 = owner->y0; aPer = owner->zPer; bPer = owner->yPer; }
	if (self == C_BOTTOM || self == C_TOP) { a0 = owner->x0; b0 = owner->z0; aPer = owner->xPer; bPer = owner->zPer; }
	if (self == C_REAR || self == C_FRONT) { a0 = owner->x0; b0 = owner->y0; aPer = owner->xPer; bPer = owner->yPer; }
	//###20121228done

	int aface = a - a0; //face is a coordinate vector in the interface!
	int bface = b - b0;

	int n = bface * aPer + aface; //implicit array coordinates utilizing the length of the edge of the rectangular interface aPer

	grainP g = owner->grainOfCell(infector);

	long arrayPos = g->arrayPos;

	if (interCells[n] != 0) return 0;

	interCells[n] = arrayPos;

	if (this->protocol == STATISTICAL)
		bucket->addOri(arrayPos, owner->grains.ngrains);

	somethingToSend = SOMETHING;

	return 1;
}

long * worldInterface::RLEencoding( long * size )
{
	
	long i = 0;
	long j = 0;
	long k = 0;

	long * tempG = (long *) malloc(2 * count * sizeof(long));
	if( !tempG ) err->reportError(ERRTXT( "Cannot allocate more memory" ) );
	
	i = 0;
	j = 0;
	k = 0;

	tempG[0] = interCells[0];

	if( tempG[0] < 0 ) tempG[0] = -1;

	while( i < count )
	{
		int tmp = tempG[k];
		int v = interCells[i];

		if( tmp < 0 ) tmp = -1;
		if( v < 0 ) v = -1;

		if( v == tmp ){
			j++;
			i++;
		}
		else
		{
			k++;
			tempG[k] = j; k++;
			tempG[k] = v;
			j=0;
		}

		if( i == count )
		{
			k++;
			tempG[k] = j;
		}
	}
	*size = k+1;

	long * resultG = (long *) calloc( *size, sizeof(long) );

	if( !resultG ) err->reportError( ERRTXT( "Cannot allocate more memory" ) );

	memcpy( resultG, tempG, *size * sizeof(long) );

	free( tempG );

	return resultG;
}

int * worldInterface::prepareDataToSend( int * n )
{
	return NULL;
}

#ifdef RLE

void worldInterface::receiveDataSendData( short source, long sizeToRecv )
{
	MPI_Status st;
	long size = 0;

	long * encodedCells = NULL;

	if( sizeToRecv > 0 )
	{
		recvCells = (long *) calloc( sizeToRecv, sizeof(long) );

		if( !recvCells ) err->reportError( ERRTXT( "Cannot allocate more memory" ) );

		MPI_Recv( recvCells, sizeToRecv, MPI_LONG, source, SENDRECV, MPI_COMM_WORLD, &st );

		owner->somethingReceived = SOMETHING;
		somethingReceived = SOMETHING;
		recvSize = sizeToRecv;
	}

	if( somethingToSend )
	{
		if( protocol == DETERMINISTIC )	encodedCells = RLEencoding( &size );
		if( protocol == STATISTICAL )	encodedCells = bucket->countOrientations( &size );
	}

	MPI_Send( &size, 1, MPI_LONG, source, SEND, MPI_COMM_WORLD );

	if( somethingToSend )
	{
		MPI_Send( encodedCells, size, MPI_LONG, source, SEND, MPI_COMM_WORLD );

		for( long i = 0; i < count; i++ ) if( interCells[i] > 0 ) interCells[i] *= -1;

		somethingToSend = NOTHING;
		free( encodedCells );
	}
}

void worldInterface::sendDataRecvData( void )
{
	QUICKASSERT( owner );

	if( nextNode == owner->myNode() )	return;
	
	MPI_Status sta;
	long sizeToRecv = 0;

	long size = 0;
	long * encodedCells = NULL;

	if( somethingToSend )
	{
		if( protocol == DETERMINISTIC ) encodedCells = RLEencoding( &size );
		if( protocol == STATISTICAL )	encodedCells = bucket->countOrientations( &size );
	}

	long message[2] = { (long) self, size };

	MPI_Send( message, 2, MPI_LONG, nextNode, STATUSSEND, MPI_COMM_WORLD );

	if( somethingToSend )
	{
		MPI_Sendrecv( encodedCells, size, MPI_LONG, nextNode, SENDRECV, &sizeToRecv, 1, MPI_LONG, nextNode, RECV, MPI_COMM_WORLD, &sta );

		for( long i = 0; i < count; i++ ) if( interCells[i] > 0 ) interCells[i] *= -1;

		somethingToSend = NOTHING;

		free( encodedCells );
	}else
		MPI_Recv( &sizeToRecv, 1, MPI_LONG, nextNode, RECV, MPI_COMM_WORLD, &sta );

	MPI_Status stb;

	if( sizeToRecv > 0 )
	{
		recvCells = (long *) calloc( sizeToRecv, sizeof(long) );

		if( !recvCells ) err->reportError( ERRTXT( "Cannot allocate more memory" ) );

		MPI_Recv( recvCells, sizeToRecv, MPI_LONG, nextNode, RECV, MPI_COMM_WORLD, &stb );
		owner->somethingReceived = SOMETHING;
		somethingReceived = SOMETHING;
		recvSize = sizeToRecv;
	}
}

/*void worldInterface::distributeImportedCells100( void )
{
	long c0 = 0;
	long cPer = 0; //###only for readability?
	long aPer = 0;
	long bPer = 0;

	long xPer = owner->xPer; //dimensions of the "owner" subdomain box 
	long yPer = owner->yPer;
	long zPer = owner->zPer;
	long x0 = owner->x0; //global edge coordinates of the "owner" subdomain box
	long y0 = owner->y0;
	long z0 = owner->z0;

	if( self == C_LEFT )		{ aPer = zPer; bPer = yPer; cPer = xPer; c0 = x0; }

	if( self == C_RIGHT )		{ aPer = zPer; bPer = yPer; cPer = xPer; c0 = x0 + xPer - 1; }

	if( self == C_BOTTOM )		{ aPer = xPer; bPer = zPer; cPer = yPer; c0 = y0; }

	if( self == C_TOP )			{ aPer = xPer; bPer = zPer; cPer = yPer; c0 = y0 + yPer - 1; }

	if( self == C_FRONT )		{ aPer = xPer; bPer = yPer; cPer = zPer; c0 = z0; }

	if( self == C_REAR )		{ aPer = xPer; bPer = yPer; cPer = zPer; c0 = z0 + zPer - 1; }

	caHdlP automaton = owner;

	QUICKASSERT( recvCells );

	QUICKASSERT( recvSize > 0 );

	int isum = 0;

	for( int i = 0; i < recvSize; i+=2 )
	{
		isum += recvCells[i+1];

		QUICKASSERT( recvCells[i+1]>0 );

		if( recvCells[i] <= 0 ) continue;

		int rep = recvCells[i+1];

		while( rep > 0 )
		{
			int oriIndex = recvCells[i];

			int n = isum - rep ;

			long ix = 0; 
			long iy = 0;
			long iz = 0;

			int b0 = n / aPer;
			int a0 = n % aPer; //conversion of implicit array coordinates in explicit interface planar coordinates

			if( self == C_LEFT || self == C_RIGHT ) { ix = c0; iy = b0 + y0; iz = a0 + z0; }
			if( self == C_TOP || self == C_BOTTOM ) { ix = a0 + x0; iy = c0; iz = b0 + z0; }
			if( self == C_FRONT || self == C_REAR ) { ix = a0 + x0; iy = b0 + y0; iz = c0; }
			//###20121228done assigns domain coordinates

			cellP c = automaton->getCellPer( ix, iy, iz );

			if( c ) //growth increment into already infected cells "on the other site" gets lost
			{
				rep--;
				continue;
			}

			Real phi1 = automaton->grains.orientations[oriIndex]->phi1;
			Real PHI = automaton->grains.orientations[oriIndex]->PHI;
			Real phi2 = automaton->grains.orientations[oriIndex]->phi2;

			cellP newCell = automaton->newNucleus( phi1, PHI, phi2, ix, iy, iz, 0 );

			rep--;
		}
	}
	free( recvCells );
	somethingReceived = NOTHING;
	recvCells = NULL;
	recvSize = 0;
}*/

void worldInterface::distributeImportedCells100( void )
{
	long c0 = 0;
	long cPer = 0; 
	long aPer = 0;
	long bPer = 0;

	long xPer = owner->xPer;  //dimensions of the "owner" subdomain box 
	long yPer = owner->yPer;
	long zPer = owner->zPer;
	long x0 = owner->x0;//global edge coordinates of the "owner" subdomain box
	long y0 = owner->y0;
	long z0 = owner->z0;

	if( self == C_LEFT )		{ aPer = zPer; bPer = yPer; cPer = xPer; c0 = x0; }

	if( self == C_RIGHT )		{ aPer = zPer; bPer = yPer; cPer = xPer; c0 = x0 + xPer - 1; }

	if( self == C_BOTTOM )		{ aPer = xPer; bPer = zPer; cPer = yPer; c0 = y0; }

	if( self == C_TOP )			{ aPer = xPer; bPer = zPer; cPer = yPer; c0 = y0 + yPer - 1; }

	if( self == C_FRONT )		{ aPer = xPer; bPer = yPer; cPer = zPer; c0 = z0; }

	if( self == C_REAR )		{ aPer = xPer; bPer = yPer; cPer = zPer; c0 = z0 + zPer - 1; }

	caHdlP automaton = owner;

	QUICKASSERT( recvCells );

	QUICKASSERT( recvSize > 0 );

	long isum = 0;

	for( long i = 0; i < recvSize; i+=2 )
	{
		isum += recvCells[i+1];

		QUICKASSERT( recvCells[i+1]>0 );

		if( recvCells[i] <= 0 ) continue;

		long rep = recvCells[i+1];

		while( rep > 0 )
		{
			long arrayPos = recvCells[i];

			long n = isum - rep ;

			long ix = 0; 
			long iy = 0;
			long iz = 0;

			long b0 = n / aPer;
			long a0 = n % aPer; //conversion of implicit array coordinates in explicit interface planar coordinates

			if( self == C_LEFT || self == C_RIGHT ) { ix = c0; iy = b0 + y0; iz = a0 + z0; }
			if( self == C_TOP || self == C_BOTTOM ) { ix = a0 + x0; iy = c0; iz = b0 + z0; }
			if( self == C_FRONT || self == C_REAR ) { ix = a0 + x0; iy = b0 + y0; iz = c0; }
			//###20121228done assigns domain coordinates

			cellP c = automaton->getCellPer( ix, iy, iz );

			if( c ) //growth increment into already infected cells "on the other site" gets lost
			{
				rep--;
				continue;
			}

			grainP g = automaton->grains.all[arrayPos];

			cellP newCell = automaton->addTransferredCell( g, ix, iy, iz );

			rep--;
		}
	}
	free( recvCells );
	somethingReceived = NOTHING;
	recvCells = NULL;
	recvSize = 0;
}

void worldInterface::distributeImportedCells( char nucleationSite ) //###LB: Not yet modified
{
	long c0 = 0;
	long cPer = 0;
	long aPer = 0;
	long bPer = 0;

	long xPer = owner->xPer;
	long yPer = owner->yPer;
	long zPer = owner->zPer;
	long x0 = owner->x0;
	long y0 = owner->y0;
	long z0 = owner->z0;

	if( self == C_LEFT )		{ aPer = zPer; bPer = yPer; cPer = xPer; c0 = x0; }

	if( self == C_RIGHT )		{ aPer = zPer; bPer = yPer; cPer = xPer; c0 = x0 + xPer - 1; }

	if( self == C_BOTTOM )		{ aPer = xPer; bPer = zPer; cPer = yPer; c0 = y0; }

	if( self == C_TOP )			{ aPer = xPer; bPer = zPer; cPer = yPer; c0 = y0 + yPer - 1; }

	if( self == C_FRONT )		{ aPer = xPer; bPer = yPer; cPer = zPer; c0 = z0; }

	if( self == C_REAR )		{ aPer = xPer; bPer = yPer; cPer = zPer; c0 = z0 + zPer - 1; }

	caHdlP automaton = owner;

	QUICKASSERT( recvCells );

	QUICKASSERT( recvSize > 0 );

	long isum = 0;

	for( long i = 0; i < recvSize; i+=2 )
	{
		isum += recvCells[i+1];

		QUICKASSERT( recvCells[i+1]>0 );

		if( recvCells[i] <= 0 ) continue;

		long rep = recvCells[i+1];

		while( rep > 0 )
		{
			int oriIndex = recvCells[i];

			long ix = 0; 
			long iy = 0;
			long iz = 0;

			long a0 = (long) ( ( aPer - 1 ) * owner->myRandom.parkMiller() );
			long b0 = (long) ( ( bPer - 1 ) * owner->myRandom.parkMiller() );

			if( nucleationSite == BULK ) c0 = (long) ( ( cPer - 1 ) * owner->myRandom.parkMiller() );

			if( self == C_LEFT || self == C_RIGHT ) 
			{ 
				ix = c0;

				if( nucleationSite == BULK ) ix += x0;

				iy = b0 + y0; 
				iz = a0 + z0; 
			}

			if( self == C_TOP || self == C_BOTTOM ) 
			{ 
				ix = a0 + x0; 
				iy = c0;

				if( nucleationSite == BULK ) iy += y0;

				iz = b0 + z0; 
			}

			if( self == C_FRONT || self == C_REAR ) 
			{ 
				ix = a0 + x0; 
				iy = b0 + y0; 
				iz = c0;

				if( nucleationSite == BULK ) iz += z0;
			}

			cellP c = automaton->getCellPer( ix, iy, iz );

			if( c )
			{
				rep--;
				continue;
			}

			/*Real phi1 = automaton->grains.orientations[oriIndex]->phi1;
			Real PHI = automaton->grains.orientations[oriIndex]->PHI;
			Real phi2 = automaton->grains.orientations[oriIndex]->phi2;

			cellP newCell = automaton->newNucleus( phi1, PHI, phi2, ix, iy, iz, 0 );*/


			grainP g = automaton->grains.all[oriIndex];

			cellP newCell = automaton->addTransferredCell( g, ix, iy, iz );

			rep--;
		}
	}
	free( recvCells );
	somethingReceived = NOTHING;
	recvCells = NULL;
	recvSize = 0;
}

#else

void worldInterface::receiveDataSendData( int source, int message )
{
	MPI_Status st;
	char somethingToReceive = (char) message;

	if( somethingToReceive )
	{
		recvCells = (int *) calloc( count, sizeof(int) );
		if( !recvCells ) err->reportError( ERRTXT( "Cannot allocate more memory" ) );

		MPI_Recv( recvCells, count, MPI_INT, source, SENDRECV, MPI_COMM_WORLD, &st );
		owner->somethingReceived = SOMETHING;
		somethingReceived = SOMETHING;
	}

	MPI_Send( &somethingToSend, 1, MPI_CHAR, source, SEND, MPI_COMM_WORLD );

	if( somethingToSend )
	{
		MPI_Send( interCells, count, MPI_INT, source, SEND, MPI_COMM_WORLD );
		for( int i = 0; i < count; i++ ) interCells[i] *= -1;
		somethingToSend = NOTHING;
	}
}

void worldInterface::sendDataRecvData( void )
{
	QUICKASSERT( owner );
	if( nextNode == owner->myNode() )
		return;
	
	MPI_Status sta;
	char somethingToRecv = NOTHING;
	int message[2] = { self, (int) somethingToSend };	

	MPI_Send( message, 2, MPI_INT, nextNode, STATUSSEND, MPI_COMM_WORLD );

	if( somethingToSend )
	{
		MPI_Sendrecv( interCells, count, MPI_INT, nextNode, SENDRECV, &somethingToRecv, 1, MPI_CHAR, nextNode, RECV, MPI_COMM_WORLD, &sta );

		for( int i = 0; i < count; i++ ) if( interCells[i] > 0 ) interCells[i] *= -1;

		somethingToSend = NOTHING;
	}else
		MPI_Recv( &somethingToRecv, 1, MPI_CHAR, nextNode, RECV, MPI_COMM_WORLD, &sta );

	MPI_Status stb;

	if( somethingToRecv )
	{
		recvCells = (int *) calloc( count, sizeof(int) );
		if( !recvCells ) err->reportError( ERRTXT( "Cannot allocate more memory" ) );

		MPI_Recv( recvCells, count, MPI_INT, nextNode, RECV, MPI_COMM_WORLD, &stb );
		owner->somethingReceived = SOMETHING;
		somethingReceived = SOMETHING;
	}
}

void worldInterface::distributeImportedCells( void )
{
	long c0 = 0;
	long cPer = 0;
	long aPer = 0;
	long bPer = 0;

	long xPer = owner->xPer;
	long yPer = owner->yPer;
	long zPer = owner->zPer;
	long x0 = owner->x0;
	long y0 = owner->y0;
	long z0 = owner->z0;

	if( self == C_LEFT )		{ aPer = zPer; bPer = yPer; cPer = xPer; c0 = x0; }

	if( self == C_RIGHT )		{ aPer = zPer; bPer = yPer; cPer = xPer; c0 = x0 + xPer - 1; }

	if( self == C_BOTTOM )	{ aPer = xPer; bPer = zPer; cPer = yPer; c0 = y0; }

	if( self == C_TOP )		{ aPer = xPer; bPer = zPer; cPer = yPer; c0 = y0 + yPer - 1; }

	if( self == C_FRONT )		{ aPer = xPer; bPer = yPer; cPer = zPer; c0 = z0; }

	if( self == C_REAR )		{ aPer = xPer; bPer = yPer; cPer = zPer; c0 = z0 + zPer - 1; }

	caHdlP automaton = (caHdlP) owner->self;

	QUICKASSERT( recvCells );

	for( int i = 0; i < count; i++ )
	{
		int oriIndex = recvCells[i];

		if( oriIndex <= 0 ) continue;

		long ix = 0; 
		long iy = 0;
		long iz = 0;

		int b0 = i / aPer;
		int a0 = i % aPer;

		if( self == C_LEFT || self == C_RIGHT ) { ix = c0; iy = b0 + y0; iz = a0 + z0; }
		if( self == C_TOP || self == C_BOTTOM ) { ix = a0 + x0; iy = c0; iz = b0 + z0; }
		if( self == C_FRONT || self == C_REAR ) { ix = a0 + x0; iy = b0 + y0; iz = c0; }

		cellP c = automaton->getCellPer( ix, iy, iz );

		if( c ) continue;

		Real phi1 = automaton->grains.orientations[oriIndex]->phi1;
		Real PHI = automaton->grains.orientations[oriIndex]->PHI;
		Real phi2 = automaton->grains.orientations[oriIndex]->phi2;

		cellP newCell = automaton->newNucleus( phi1, PHI, phi2, ix, iy, iz, 0 );
	}
	free( recvCells );
	somethingReceived = NOTHING;
	recvCells = NULL;
}

#endif