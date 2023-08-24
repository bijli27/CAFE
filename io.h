#ifndef _io_h_		//header guard
#define _io_h_	

#include <iostream> 		//header file for input output stream
#include <string>			//header file for string
#include <fstream> 			//header file for file stream
#include <stdio.h>			//header file for standard input output
#include <stdlib.h>			//header file for standard library
#include <string.h>			//header file for string
#include "applic.h"			//header file for applic
#include <stddef.h>			//header file for standard definition

#define BUFSIZE 1024 						//buffer size for reading data from file
#define GGAPPROACH_MAXCASIZE 64 			//maximum size of CA

#define IDTAG0 0x7769						//ID tag1 for file format
#define IDTAG1 0x6e746572					//ID tag2 for file format
#define IDTAG2 0x6973636f					//ID tag3 for file format
#define IDTAG3 0x6d696e67					//ID tag4 for file format

#define _NO_ 0x00							//no
#define _YES_ 0x01							//yes

//#define INCLUDECELLSIZE true  //option to allow for cell size commitment to perform digital line interception with mstoPNG
//(set true for mstoGrainSize compatible binaries, in which case you also have to set it to true in mstoPNG's mstoPNG.cpp!!)

#define xWinLoc (30)				//x location of window
#define yWinLoc (50)				//y location of window

#define xWinSize (1000)				//x size of window
#define yWinSize (700)				//y size of window

#define xWinCentre (xWinSize/2)		//x centre of window
#define yWinCentre (yWinSize/2)		//y centre of window

#define stringsNotEqual(a,b) strcmp(a,b)		//compare two strings
#define stringsEqual(a,b) !stringsNotEqual(a,b)	//compare two strings

#define READ   0			//read
#define WRITE  1			//write
#define APPEND 2			//append
#define FAILURE 0			//failure
#define SUCCESS 1			//success

class dataLine;
class dataBlock; 
typedef dataLine * dataLineP;
typedef dataBlock * dataBlockP;
class io;
typedef io * ioP;
class error;
typedef error * errorP;

using namespace std;

typedef struct			//structure for data line
{
	union
	{
		float f;
		long i;
		char *s;
	} d;
	char type;

} univData;
typedef univData * univDataP;		//pointer to universal data

class dataLine 			//class for data line
{
public:
	dataLine( void );		//constructor of dataLine class
	dataLine( int size );	//constructor overloading of dataLine class
	~dataLine( void );		//destructor of dataLine class

	univDataP dat;
	long dataCount;

	dataLineP next;
	dataLineP prev;
};

class dataBlock 			//class for data block
{
public:
	dataBlock( void ); 			//constructor of dataBlock class
	~dataBlock( void );			//destructor of dataBlock class
	long columnCount;
	long lineCount;
	
	dataLineP first;		//pointer to first data line
	dataLineP last;			//pointer to last data line

	char head[BUFSIZE];		
	char name[128];
};

// MPI_IO structs

typedef struct					//calculation not platform independent!
{
    int OriIndex;               // 4 byte
    float phi1,PHI,phi2;        // 3*4 byte  //Ori Output only single precision ###!
} MPI_IO_OriIndex  ;            //// 16 byte

typedef struct
{
    int MPIRank;              // 4 byte
    int x0, xmax, y0, ymax, z0, zmax;   // 6*4 byte
} MPI_IO_Node;                  //// 28 byte

typedef struct
{
    int OriIndex;               // 4 byte
    char state;                 // 1 byte
} MPI_IO_Cell;                  //// 5 byte

typedef struct
{
	float InfProgress;			// 4 byte  ###coded as 255*(1- c->rxFraction/1) (black is deformed, white is 100% recrystallized)
	char state;					// 1 byte
} MPI_IO_CellInfState;

typedef struct
{
	double phi1_TB;                	// 8 byte
	double PHI_TB;                 		
	double phi2_TB;                         	
	double RVFac3;					// 3*8 byte
	double scatter;					// 1*8 byte
	double rhoGND;               
	double rhoMobil0;            	// 2*8 byte
	double rhoIntern0;           
	double rhoWall0;            
    double rhoM;                 
	double rhoI;                 			
	double rhoW;                		// 3*8 byte
	double NGLS;                
	double rhoTotal;  				// 6*8 byte
	double phi1;
	double PHI;
	double phi2;                   // 8 bytes * 17
	

	short Oriindex; 
	short xleft;
	short xright;
	short yleft;
	short yright;
	short zleft;
	short zright; // 2 bytes * 7

	long alreadyNuc;
	long DeformGrainIndex;
	long CellCount;		// 4 bytes * 3
	long arrayPos;
	long initialCellNumber;


	int oriIndex;   // 4 bytes * 7
	int nuc1; 
	int nuc2; 
	int nuc3;
	int colR;
	int colG;
	int colB;
} MPI_IO_GrainInfo;

// ---

class io
{

public:
	io( void  );		//constructor of io class
	~io( void );		//destructor of io class

	short open( short rw, const char * filename, FILE ** file );		//open file for reading or writing
	short open( const char * filename, ofstream * file );				//open file for writing (overloading)
	short open( const char * filename, ifstream file );					//open file for reading (overloading)
	short write( char * message, FILE * file );							//write message to file
	short write( string message, ofstream * file );

	dataBlockP readDataBlock( const char *name, const char *fileName );		//read data block from file (overloading)
	dataBlockP readDataBlock( const char *name, FILE * file );				//read data block from file (overloading)

	char *newString(const char *s);						//create new string
	Real getReal( dataLineP line, long column );		//get real number from data line
	Real getReal2( dataLineP line, long column );		
	long getInt( dataLineP line, long column );			//get integer from data line
	char * getString( dataLineP line, long column );	//get string from data line

	Real geTReal( const char *s, dataBlockP db );   		//get real number from data block
	long geTInt( const char *s, dataBlockP db );			//get integer from data block
	char * geTString( const char *s, dataBlockP db );		//get string from data block

};

class error
{
public:
	
	error( ofstream * log = NULL, int rank = 0 ) :			//constructor of error class
	  logFile(log), myRank(rank)							//set log file
	  {
	  }

	~error( void );							//destructor of error class

	void reportError( const char * message, int rank = 0 );			//report error
	void reportWarning( const char * message, int rank = 0 );		//report warning

	ofstream * logFile;		
	io ioHdl;		//io handler
	int myRank;		//rank of process
};

#endif

