#ifndef _caHdl_h_
#define _caHdl_h_

#include "simhdl.h"
#include "polycrystal.h"
#include <stdio.h>
#include <list>
#include <math.h>
#include "applic.h"
#include "tempControl.h"
#include "grainCA.h"
#include <time.h>
#include <exception>
#include <iomanip>
#include "gaussian.h"

/*#ifdef CORe_MPI
#include <mpi.h>
#endif*/

#ifdef _MEMLEAK_
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
	#include <crtdbg.h>
#endif

inline time_t currentTime(time_t *t) {return time(t);}

#ifndef Real
typedef double Real;
#endif

#define DRX 1e-33
#define NOTNEWCELLCREATED NULL
#define GRAINSDEFINED 0
#define VOLUMEDEFINED 1
#define MINIMUM(a,b,c) (a < b ? (a < c ? a : c) : (b < c ? b : c))
#define MAXIMUM(a,b,c) (a > b ? (a > c ? a : c) : (b > c ? b : c))
#define NUCLEUS 0x200020
#define ORICLASSMAX			20 
#define SOLIDMASK 0xFFFFFFF9
#define LIQUIDMASK 0xFFFFFFFA 
//#define MELTINGPOINT (1375) //as: reading from input file now

class caHdl : public simHdl, public polyCrystal
{
public:
	caHdl( void );
	~caHdl( void );

	Real initialDeltaTime(void);  //used

	void createMPItypes( void );
	char handlePeriodicity( long * ix, long * iy, long * iz );
	cellP getCellPer( long ix, long iy, long iz );
	cellP * getCellPerP( long ix, long iy, long iz );

	grainP getOrigGrain( long ix, long iy, long iz );
	grainP * getGrainPos( long ix, long iy, long iz );

	cellP infectCell( cellP infector, long dx, long dy, long dz, int infectorvictimindex, Real overshoot ); //###MK20121211, 26NN
	cellP addTransferredCell( grainP g, long ix, long iy, long iz );
	cellP newNucleusSol(Real p1, Real ps, Real p2, long ix, long iy, long iz, TCellP Tcell, long DeformGrainIndex, int Oriindex, Real TimeToMelt);
	
	void retireCell( cellP c );
	int sendToInterface( cellP infector, long dx, long dy, long dz, int nextNode, Real overshoot );
	void sendReceive( void );
	void receiveSend( int sender );
	void receiveSend( void );
	bool isMyRankInTheList(short* list, int size);

	inline grainP grainOfCell( cellP c ) { return grains.all[(long) ((c->dataHolder & 0x7FFFFFFFFFF)>>22)]; }//rx grain
	inline grainP defGrainOfCell( cellP c ) {return grains.all[(long) (c->dataHolder>>43)];}// change this to store tempearature or something else.
	inline long posOfGrainOfCell( cellP c ) { return (long) ((c->dataHolder & 0x7FFFFFFFFFF)>>22); }
	inline int cellInfector( cellP c ) { return (int) (c->dataHolder & 0x1F); }
	inline grainP getOrigGrain( cellP c ) { return grains.all[grainArrayPosition(c)]; }
	
	inline char isCellRecrystallized( cellP c ) { if( c->rxFrac >= 1000.0 ) return 1; return 0; }
	inline Real returnPositiveAngle(Real angle) { while (angle < 0) { angle += (2 * _PI_); } return angle; };
	void freeCellMemory(cellP c);
	long calcGrowthStepMG(void);
	void updateFullRX(void);

	void getGeneralInput( Real * Ncells );
	void initializeSolCA( const char * inputfile, int iautomaton );
	void initializeGrainsSolCA( void );
	void prepareAutomata( void );
	
	void determineCAsizeS(Real Ncells);
	void determineNodeTopology( void );
	void startParallelProcesses( int argc, char *argv[] );
	void endParallelProcesses( void );
	void startSimulation( void );

	void initializeCells( void );
	void insertNucleiSolidification(void);
	void insertBulkNucleiSolidification(void);
	void insertSurfaceNucleiSolidification(void);
	void insertEpitaxialNucleationSolidification(void); //for epitaial nucleation
	void addOrientation( Real phi1, Real PHI, Real phi2 );
	void distributeReceivedCells( void );
	void determineInitialTexture( void );
	grainP addGrain( MPI_IO_GrainInfo * grainData );
	void prepareGrainForBroadcast( cellP newCell, MPI_IO_GrainInfo * grainData );
	void prepareGrainForBroadcast( long arrPos, MPI_IO_GrainInfo * grainData );
	
	char isParticle( cellP c );
	
	void addNucleus( cellsToNucleate * cn, short ix, short iy, short iz, grainP rxg );//It could be useful in the future
	//void prepareNucleation( void );//Not used but it could be used in the future
	
	Real currentCellTemperature(cellP c, Real time, long timeKey);

	//////////////USED BUT NOT NEEDED///////////////////////////
	cellP rhoFacCP(float x); //Used but not needed
	cellP updateFacCP( float x, cellP c ); //Used but not needed
	Real nucleationProbability( Real rhoW );//Used but not needed
	cellP newNucleusMG(Real p1, Real ps, Real p2, long ix, long iy, long iz, Real Temperature, Real rT, Real rGND, double rm0, double ri0, double rw0, int n1, int n2, int n3, Real NGLS, long DeformGrainIndex, int Oriindex); //used but not used
	
	//////////////OUTPUT///////////////////

    int outCount, oriOutCount;;
    MPI_Datatype MPI_IO_OriIndex_Type;
    MPI_Datatype MPI_IO_Node_Type;
    MPI_Datatype MPI_IO_Cell_Type;
	MPI_Datatype MPI_IO_CellInfState_Type;
	MPI_Datatype MPI_IO_GrainInfo_Type;
	MPI_Datatype MPI_IO_Nucleus_Type;
	MPI_Datatype MPI_IO_Particle_Type;
	MPI_Datatype MPI_IO_Ori_Type;

    bool userDefinedMSOutput;
	bool generateMSBinary;
	char MSBinaryPlotMode;

    list<double> userDefinedFrequencies;    //MS Frequencies!
    list<double> userDefinedOriFrequencies;
    
	void Output(void);
	void OriOutput( void );
    void WriteOri( void );
	void RxFracOutput( void );

    void MicrostructureOutput( void );
    void WriteMicrostructurePRNGRGB( void );
	void WriteMicrostructureInfGreyscale( void );

	void write_voxeldata_coloring_ipfz(void);
	void write_xsection_coloring_ipfz(void);
	void write_voxeldata_oris(void);
	char InsufficientMelting(cellP c); // #as: for porosity
	void write_prorosity(void);
	//void write_prorosity(vector<porosity> nomelt); //#as: for porosity writing
	//private:
	//	vector<porosity> nomelt;

	public:
	
	long * compressNodeInformation( long * size );

	void saveMicrostructure( const char * filename );
	void writeHeader( MPI_Offset * offset );
	nucleiPacket **  writeGrains( MPI_Offset * offset );
	void writeOrientations( MPI_Offset * offset );
	void writeNuclei( MPI_Offset * offset, nucleiPacket *** packets );
	void writeCells( MPI_Offset * offset );

	////////TEMPERATURE MANAGEMENT//////////////  //Maybe it would be a good idea to define a temperature handler

	void readTemperatureTimeDiscretization(dataBlockP tTData);
	TNodeRawData* readTemperatureNodes(dataBlockP Tnodes);
	cellNodesRawData* readTemperatureCells(dataBlockP Tcells);
	void readTemperatureCellsDiscretization(dataBlockP xdiscr, dataBlockP ydiscr, dataBlockP zdiscr);
	void initializeTemperatureCells(void);
	void establishCellOwnershipToTemperatureCell(void);
	inline tempCellP getTCellofCell(cellP c) {
		long index = (((long)(c)) >> 15);
			return Tcells->all[index];
	};
	inline Real getNucTimeInactiveCell(cellP c) { return Tcells->getCellNucleationTime(c); }
	inline TCellP TCellOfActiveCell(cellP c) { return Tcells->all[(long)(c->dataHolder >> 43)]; }

	////////MEMORY MANAGEMENT////////////////

	cellP getCellMemory(void);

	////SEPI DATA FUNCTIONS.
	void orientations_SEPI(void); // write id , phi_1, phi. phi_2.
	void write_state_SEPI(void);// writes whether crystallized, powder or liquid.
	void write_final_state_SEPI(void);
	void write_voxeldata_oris_and_oris_id_SEPI(void);
	void write_nucleus_data_SEPI(void); // write the x, y, z, phi_1,phi, phi_2, id
	void writingVoxelOrientationID_SEPI();
	void draw_pole_figure(void);
	
public:
	Real Vmax = 0; // maximum cryst growth speed
	int totalNodes;
	MPI_File theFile;
	ofstream logFile;		//Log file; a real log file.
	FILE* core8;		//Is read by many instances
	double* grainOriNuc;
	double numberInputNucOris;
	//double numberInputOris; //#rg: Not needed.
	// long tempCounter; // #rg :No use
	
	cellPool cells;				// struct									
	physData parm;				//struct				
	genData gen;				//struct
	cpData box;					//struct
	//allRho rhoDot;
	porosity PT;
	randomClass myRandom;// class
	guassianclass nucprobability;  // gaussian class
	errorP err; //class
	tempCellPoolP Tcells; //class tempCellPool
	grainPool grains;// class grainPool
	processHdlP parallelHdl;// class
	oriBucketP texture; // class

	//grainPool tempgrains; //#rg: Not needed. Never used
	//grainPool Idealgrains; //#rg: Not needed. Never used.	
	//boundariesPool boundaries;

	trigger transfer;
	trigger txtOutput;
	trigger rxOutput;
	trigger mstOutput;
	trigger rvOutput;



};

typedef caHdl * caHdlP;

#endif