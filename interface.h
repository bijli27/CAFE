#pragma once
#ifndef _interface_h_
#define _interface_h_

#include <stdlib.h>
#include <stdint.h>
#define RLE

//Changed

/*
 this struct should be kept small!! 
 LB: This new cell struct is "only" 3 bytes larger than the original one but contains much more info
*/
struct cellStruct;
typedef cellStruct cell;
typedef cell * cellP;

#define C_RIGHT				0
#define C_FRONT				1
#define C_TOP				2
#define C_BOTTOM			5
#define C_REAR				4
#define C_LEFT				3

#define C_FRONTTOP			6
#define C_FRONTBOTTOM		7
#define C_FRONTLEFT			8
#define C_FRONTRIGHT		9
#define C_REARTOP			10
#define C_REARBOTTOM		11
#define C_REARLEFT			12
#define C_REARRIGHT			13
#define C_RIGHTTOP			14
#define C_RIGHTBOTTOM		15
#define C_LEFTTOP			16
#define C_LEFTBOTTOM		17 //6-<18

#define C_FRONTTOPLEFT		18
#define C_FRONTTOPRIGHT		19
#define C_FRONTBOTTOMLEFT	20
#define C_FRONTBOTTOMRIGHT	21
#define C_REARTOPLEFT		22
#define C_REARTOPRIGHT		23
#define C_REARBOTTOMLEFT	24
#define C_REARBOTTOMRIGHT	25 //18-<26


//###20121228done #defines 
#define SIXFACES 6
#define TWENTYSIXFACES 26
//#define TWELVEEDGES 18
//#define EIGHTCORNERS 26

#define NOTHING (0x00)
#define SOMETHING (0x01)
#define NOTHINGTOINFECT (0)

#define _RXFRACUNIT 0xFF

#define NOTASSIGNED 0x0101
#define RETIRED 0x0

#define BUCKETSIZE 100
#define NO_ORI -1

class caHdl;
typedef caHdl * caHdlP; // caHdlP is a caHdl datatype

class interfaceHdl;
typedef interfaceHdl * interfaceHdlP; //interfaceHdlP is a interfaceHdl datatype

class worldInterface;
typedef worldInterface * worldInterfaceP; //worldInterfaceP is a worldInterface datatype

class error;
typedef error * errorP; //errorP is a error datatype

struct oriCounter{
	long count;
	long oriIndex;
};
typedef oriCounter * oriCounterP;		//oriCounterP is a oriCounter datatype

struct oriBin{
	long len;
	long maxLen;
	oriCounter entry[];
};
typedef oriBin * oriBinP;

class oriBucket
{
public :
	oriBucket( worldInterfaceP own = NULL );		//constructor
	~oriBucket( void );							//destructor
	
	void addOri( long oriIndex, long maxOri );
	long * countOrientations( long * size );
	void emptyBucket(void);
	
	long allOri;
	oriBinP oriContainer[BUCKETSIZE];
	worldInterfaceP owner;
};
typedef oriBucket * oriBucketP; //oriBucketP is a oriBucket datatype

class worldInterface
{
public:

	worldInterface( caHdlP own = NULL, long ix = 0, long iy = 0, short me = 0 );
	~worldInterface(void);

	void sendDataRecvData( void );
	char addCell100( cellP infector, int ix, int iy );
	int * prepareDataToSend( int * n );
	void distributeImportedCells100( void );
	void distributeImportedCells( char nucleiInsertionSite );

	//int * RLEencoding( short * size );
	long * RLEencoding( long * size );

	oriBucketP bucket;
	char protocol;

	errorP err;

#ifdef RLE
	void receiveDataSendData( short source, long message );
#else
	void receiveDataSendData( int source, int message );
#endif

	long aPer;
	long bPer;
	short self;
	long count;
	char somethingToSend;
	char somethingReceived;

	caHdlP owner;

	//int * recvCells;
	//int * interCells;

	long * interCells;
	long * recvCells;

	long sendRecv;
	int nextNode;
	long active;
	long retired;
	long recvSize;
};

#endif