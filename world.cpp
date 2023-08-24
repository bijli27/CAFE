#include "world.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

world::world(void)
{
	xPer=0;
	yPer=0;
	zPer=0;
	zFac=0;
	x0=0;
	y0=0;
	z0=0;

	boxPer = BOXPER;
	nodePer = NON_NODE_PER;
	protocol = DETERMINISTIC;
	commType = STEPS;

	memset( nodeFace, NULL, sizeof(nodeFace) );
	somethingToSend = NOTHING;
	somethingReceived = NOTHING;

	rank=0;
}

world::~world(void)
{
	//###20121228done
	for( int i=0; i< SIXFACES; i++ )
		if( nodeFace[i] ) delete nodeFace[i];
}

void world::sendReceive( void )
{
	//###20121228done add additional sendDataRecvData (//was 6)
	for( int direction=0; direction < 6; direction++ )	
		nodeFace[direction]->sendDataRecvData();
}

void world::distributeReceivedCells( void )
{
}

void world::receiveSend( void )
{

}

int world::sendToInterface( long ix, long iy, long iz, int mynode, int nextNode )
{
	return 0;
}

bool world::isMyRankInTheList(short* list, int size)
{
	return false;
}

#ifdef RLE

void world::receiveSend( int source )
{
	long envelope[2] = { -1, -1 };

	MPI_Status st;

	MPI_Recv( envelope, 2, MPI_LONG, source, STATUSSEND, MPI_COMM_WORLD, &st );

	long sender = envelope[0];
	long sizeToRecv = envelope[1];

	QUICKASSERT( sender >= 0 && sizeToRecv >= 0 );	//A message has to be received always

	//Opposite counterparts receive

	//###20121228done
	//###20120108###consider early branching strategies
	if( sender == C_RIGHT )		nodeFace[C_LEFT]->receiveDataSendData( (short) source, sizeToRecv );

	if( sender == C_FRONT )		nodeFace[C_REAR]->receiveDataSendData( (short) source, sizeToRecv );

	if( sender == C_TOP )		nodeFace[C_BOTTOM]->receiveDataSendData( (short) source, sizeToRecv );

	if( sender == C_BOTTOM )	nodeFace[C_TOP]->receiveDataSendData( (short) source, sizeToRecv );

	if( sender == C_REAR )		nodeFace[C_FRONT]->receiveDataSendData( (short) source, sizeToRecv );

	if( sender == C_LEFT )		nodeFace[C_RIGHT]->receiveDataSendData( (short) source, sizeToRecv );


	if( sender == C_FRONTTOP)		nodeFace[C_REARBOTTOM]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_FRONTBOTTOM)	nodeFace[C_REARTOP]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_FRONTLEFT)		nodeFace[C_REARRIGHT]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_FRONTRIGHT)		nodeFace[C_REARLEFT]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_REARTOP)		nodeFace[C_FRONTBOTTOM]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_REARBOTTOM)		nodeFace[C_FRONTTOP]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_REARLEFT)		nodeFace[C_FRONTRIGHT]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_REARRIGHT)		nodeFace[C_FRONTLEFT]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_RIGHTTOP)		nodeFace[C_LEFTBOTTOM]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_RIGHTBOTTOM)	nodeFace[C_LEFTTOP]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_LEFTTOP)		nodeFace[C_RIGHTBOTTOM]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_LEFTBOTTOM)		nodeFace[C_RIGHTTOP]->receiveDataSendData( (short) source, sizeToRecv );


	if( sender == C_FRONTTOPLEFT)		nodeFace[C_REARBOTTOMRIGHT]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_FRONTTOPRIGHT)		nodeFace[C_REARBOTTOMLEFT]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_FRONTBOTTOMLEFT)	nodeFace[C_REARTOPRIGHT]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_FRONTBOTTOMRIGHT)	nodeFace[C_REARTOPLEFT]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_REARTOPLEFT)		nodeFace[C_FRONTBOTTOMRIGHT]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_REARTOPRIGHT)		nodeFace[C_FRONTBOTTOMLEFT]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_REARBOTTOMLEFT)		nodeFace[C_FRONTTOPRIGHT]->receiveDataSendData( (short) source, sizeToRecv );
	if( sender == C_REARBOTTOMRIGHT)	nodeFace[C_FRONTTOPLEFT]->receiveDataSendData( (short) source, sizeToRecv );
}

#else

void world::receiveSend( int source ) //was 6 
{
	int envelope[2] = { -1, -1 };
	MPI_Status st;
	MPI_Recv( envelope, 2, MPI_INT, source, STATUSSEND, MPI_COMM_WORLD, &st );

	int sender = envelope[0];
	int message = envelope[1];

	QUICKASSERT( sender >= 0 && message >= 0 );	//A message has to be received always

	//Opposite counterparts receive

	if( sender == C_RIGHT )		nodeFace[C_LEFT]->receiveDataSendData( source, message );

	if( sender == C_FRONT )		nodeFace[C_REAR]->receiveDataSendData( source, message );

	if( sender == C_TOP )		nodeFace[C_BOTTOM]->receiveDataSendData( source, message );

	if( sender == C_BOTTOM )	nodeFace[C_TOP]->receiveDataSendData( source, message );

	if( sender == C_REAR )		nodeFace[C_FRONT]->receiveDataSendData( source, message );

	if( sender == C_LEFT )		nodeFace[C_RIGHT]->receiveDataSendData( source, message );


	if( sender == C_FRONTTOP)		nodeFace[C_REARBOTTOM]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_FRONTBOTTOM)	nodeFace[C_REARTOP]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_FRONTLEFT)		nodeFace[C_REARRIGHT]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_FRONTRIGHT)		nodeFace[C_REARLEFT]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_REARTOP)		nodeFace[C_FRONTBOTTOM]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_REARBOTTOM)		nodeFace[C_FRONTTOP]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_REARLEFT)		nodeFace[C_FRONTRIGHT]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_REARRIGHT)		nodeFace[C_FRONTLEFT]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_RIGHTTOP)		nodeFace[C_LEFTBOTTOM]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_RIGHTBOTTOM)	nodeFace[C_LEFTTOP]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_LEFTTOP)		nodeFace[C_RIGHTBOTTOM]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_LEFTBOTTOM)		nodeFace[C_RIGHTTOP]->receiveDataSendData( source, sizeToRecv );


	if( sender == C_FRONTTOPLEFT)		nodeFace[C_REARBOTTOMRIGHT]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_FRONTTOPRIGHT)		nodeFace[C_REARBOTTOMLEFT]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_FRONTBOTTOMLEFT)	nodeFace[C_REARTOPRIGHT]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_FRONTBOTTOMRIGHT)	nodeFace[C_REARTOPLEFT]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_REARTOPLEFT)		nodeFace[C_FRONTBOTTOMRIGHT]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_REARTOPRIGHT)		nodeFace[C_FRONTBOTTOMLEFT]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_REARBOTTOMLEFT)		nodeFace[C_FRONTTOPRIGHT]->receiveDataSendData( source, sizeToRecv );
	if( sender == C_REARBOTTOMRIGHT)	nodeFace[C_FRONTTOPLEFT]->receiveDataSendData( source, sizeToRecv );

}

#endif