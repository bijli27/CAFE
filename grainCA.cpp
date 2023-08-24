#include "io.h"
#include "caHdl.h"
//#include "nucleation.h"


grainPool::grainPool( caHdl * own )
{
	owner = own;
	first = NULL;
	last = NULL;
	data.count = 0;
	defGrainsCount = 0;
	NucMG = 0; // Not needed.

	data.xPer = 0;
	data.yPer = 0;
	data.zPer = 0;
	data.zFac = 0;

	Nuc = 0;

	commonGrainDiameter = -1.0; // Not needed
	subGrainDiameter = 1.0;// Not needed
	
	noris = 0;
	ngrains = 0;
	orientations = NULL;
	all = NULL;
}

long grainPool::addGrain( grainP g )
{
	if( !all )
	{
		all = new grainP[1];
		all[0] = g;
		ngrains++;
		return 0;
	}
	size_t newSize = (size_t) ngrains + 1;
	grainP * newArr = new grainP[newSize];

	if( !newArr )
		owner->err->reportError( ERRTXT("Too many grains, not enough memory") );

	memcpy( newArr, all, ngrains * sizeof(oriRepresentationP) );

	delete [] all;
	all = newArr;

	all[ngrains] = g;

	ngrains= (long) newSize;
	return ngrains - 1;
}

int grainPool::addOrientation( Real phi1, Real PHI, Real phi2 )
{
	for( int i=0; i<noris; i++ )
	{
		Real mis = owner->misorientationCubic( phi1, PHI, phi2, orientations[i]->phi1, orientations[i]->PHI, orientations[i]->phi2 );
		if( mis <= MINMIS ) 
			return i;
	}

	Real q[4] = {0};
	Real angles[3] = { phi1, PHI, phi2 };
	
	//eulerToQuaternion( phi1, PHI, phi2, q );

	unsigned char rgb[3] = { UCHAR_RANGE_MIN, UCHAR_RANGE_MIN, UCHAR_RANGE_MIN };
	double locationinsst[2] = { 0.0, 0.0 };

	owner->bunge2ipfz(phi1, PHI, phi2, rgb, locationinsst);

	owner->euler2quaternion( angles, q ); //###LB: In degrees
	
	if( !orientations ) 
	{
		orientations = new oriRepresentationP[1];
		orientations[0] = new oriRepresentation;
		orientations[0]->phi1 = phi1;
		orientations[0]->PHI = PHI;
		orientations[0]->phi2 = phi2;
		orientations[0]->q0 = q[0];
		orientations[0]->q1 = q[1];
		orientations[0]->q2 = q[2];
		orientations[0]->q3 = q[3];
		orientations[0]->count = 0;
		orientations[0]->RGBA[RED] = rgb[RED];
		orientations[0]->RGBA[GREEN] = rgb[GREEN];
		orientations[0]->RGBA[BLUE] = rgb[BLUE];
		orientations[0]->RGBA[ALPHA] = UCHAR_RANGE_MAX;

		noris++;
		return 0;
	}

	size_t newSize = (size_t)noris + 1;
	oriRepresentationP * newArr = new oriRepresentationP[newSize];

	memcpy( newArr, orientations, noris * sizeof(oriRepresentationP) );

	delete [] orientations;
	orientations = newArr;

	orientations[noris] = new oriRepresentation;
	orientations[noris]->phi1 = phi1;
	orientations[noris]->PHI  = PHI;
	orientations[noris]->phi2 = phi2;
	orientations[noris]->q0 = q[0];
	orientations[noris]->q1 = q[1];
	orientations[noris]->q2 = q[2];
	orientations[noris]->q3 = q[3];
	orientations[noris]->RGBA[RED] = rgb[RED];
	orientations[noris]->RGBA[GREEN] = rgb[GREEN];
	orientations[noris]->RGBA[BLUE] = rgb[BLUE];
	orientations[noris]->RGBA[ALPHA] = UCHAR_RANGE_MAX;
	orientations[noris]->count = 0;

	noris = newSize;
	return noris - 1;
}

grainPool::~grainPool(void)
{
	for( long i=0; i<noris; i++ ) if( orientations[i] ) delete orientations[i];
	
	delete [] orientations;

	grainP gnext = NULL;
	for( grainP g=this->first; g; g=gnext )
	{
		gnext = g->next;
		delete g;
	}

	free( all );
}

oriRepresentation::oriRepresentation( void )
{
	PHI = 0;
	phi1 = 0;
	phi2 = 0;
	q0 = 0;
	q1 = 0;
	q2 = 0;
	q3 = 0;
	count = 0;
}

oriRepresentation::~oriRepresentation( void )
{

}

void oriRepresentation::copyOriFrom(oriRepresentation* original)
{
	phi1 = original->phi1;
	PHI = original->PHI;
	phi2 = original->phi2;
	q0 = original->q0;
	q1 = original->q1;
	q2 = original->q2;
	q3 = original->q3;
	RGBA[RED] = original->RGBA[RED];
	RGBA[GREEN] = original->RGBA[GREEN];
	RGBA[BLUE] = original->RGBA[BLUE];
	RGBA[ALPHA] = original->RGBA[ALPHA];
}

grainClass::~grainClass(void)
{
	for( int i=0; i<STDLEN; i++ )
	{
		if( this->cp[i] ) 
			free(cp[i]);
	}
}

grainClass::grainClass( void )
{
	owner = NULL;
	this->Oriindex         = 0;
	this->DeformGrainIndex = 0;
	initialCellNumber = 0;

	nuc1   = 0;                                                        //SB
	nuc2   = 0;                                                        //GB
	nuc3   = 0;									                    //TB

	this->NGLS   = 0;                                                      //DZ

	RVFac3     = 0;                                                      //Recovery-Factor for cross slip, dependence on 1/NGLS

	this->rhoGND = 0;

	rhoMobil0  = 0;
	rhoIntern0 = 0;
	rhoWall0   = 0;

	rhoM       = 0;
	rhoI       = 0;
	rhoW       = 0;
	rhoTotal = 0;

	alreadyNuc   = 0;

	CellCount = 0;
	initialCellNumber = 0;

	colR = 30000 + rand()%40000;
	colG = 30000 + rand()%40000;
	colB = 30000 + rand()%40000;

	mySpace.xleft = 0;
	mySpace.xright = 0;
	mySpace.yleft = 0;
	mySpace.yright = 0;
	mySpace.zleft = 0;
	mySpace.zright = 0;

	somethingToNucleate = NOTHING;

	GBnuclei.count = 0;
	GBnuclei.first = NULL;
	GBnuclei.last= NULL;
	GBnuclei.nucleateFlag = NOTHING;

	TBnuclei.count = 0;
	TBnuclei.first = NULL;
	TBnuclei.last = NULL;
	TBnuclei.nucleateFlag = NOTHING;

	SBnuclei.count = 0;
	SBnuclei.first = NULL;
	SBnuclei.last = NULL;
	SBnuclei.nucleateFlag = NOTHING;

	PSNnuclei.count = 0;
	PSNnuclei.first = NULL;
	PSNnuclei.last = NULL;
	PSNnuclei.nucleateFlag = NOTHING;

	this->zombie = NULL;

	long i;
	for (i=0; i<HASHCOLUMNCOUNT; i++)
	{
		this->cp[i]= NULL;
	}
}

grainClass::grainClass( grainPoolP own,  Real p1, Real P, Real p2, double rhoGND, double rm0, double ri0, double rw0,int n1,int n2,int n3, short Oriindex, Real NGLS, long DeformGrainIndex )
{
	owner = own;
	this->Oriindex         = Oriindex; // this correspomds to ideal orientaion thing.
	this->DeformGrainIndex = DeformGrainIndex;
	oriIndex = owner->addOrientation( p1, P, p2 );

	nuc1   = n1;                                                        //SB
	nuc2   = n2;                                                        //GB
	nuc3   = n3;									                    //TB

	this->NGLS   = NGLS;                                                      //DZ

	RVFac3     = 1;    //NGLS;                                                   //Recovery-Factor for cross slip, dependence on 1/NGLS

	this->rhoGND = rhoGND;

	rhoMobil0  = rm0;
	rhoIntern0 = ri0;
	rhoWall0   = rw0;

	rhoM       = 0;
	rhoI       = ri0;
	rhoW       = rw0;
	rhoTotal = 0;

	alreadyNuc   = 0;

	CellCount = 0;
	initialCellNumber = 0;

	colR = 30000 + rand() % 40000;    //RND-colors
	colG = 30000 + rand() % 40000;
	colB = 30000 + rand() % 40000;

	mySpace.xleft = 0;
	mySpace.xright = 0;
	mySpace.yleft = 0;
	mySpace.yright = 0;
	mySpace.zleft = 0;
	mySpace.zright = 0;

	somethingToNucleate = NOTHING;

	GBnuclei.count = 0;
	GBnuclei.first = NULL;
	GBnuclei.last= NULL;
	GBnuclei.nucleateFlag = NOTHING;

	TBnuclei.count = 0;
	TBnuclei.first = NULL;
	TBnuclei.last = NULL;
	TBnuclei.nucleateFlag = NOTHING;

	SBnuclei.count = 0;
	SBnuclei.first = NULL;
	SBnuclei.last = NULL;
	SBnuclei.nucleateFlag = NOTHING;

	PSNnuclei.count = 0;
	PSNnuclei.first = NULL;
	PSNnuclei.last = NULL;
	PSNnuclei.nucleateFlag = NOTHING;

	if( owner->first )	owner->first->prev = this;
	this->next = owner->first;
	this->prev = NULL;
	owner->first = this;
	if( !(owner->last) )	owner->last = this;
	owner->data.count++;

	this->zombie = NULL;

	long i;
	for (i=0; i<HASHCOLUMNCOUNT; i++)
	{
		this->cp[i]= NULL;
	}
	arrayPos = own->addGrain( this );
}

grainClass::grainClass(  grainPoolP own, Real p1, Real P, Real p2, int index1, Real scatter )
{
	oriIndex = owner->addOrientation( p1, P, p2 );
	this->Oriindex  = (short) index1;
	this->scatter = scatter;
	arrayPos = own->addGrain( this );
	this->initialCellNumber = 0;
	this->CellCount = 0;

	mySpace.xleft = 0;
	mySpace.xright = 0;
	mySpace.yleft = 0;
	mySpace.yright = 0;
	mySpace.zleft = 0;
	mySpace.zright = 0;

	somethingToNucleate = NOTHING;

	GBnuclei.count = 0;
	GBnuclei.first = NULL;
	GBnuclei.last= NULL;
	GBnuclei.nucleateFlag = NOTHING;

	TBnuclei.count = 0;
	TBnuclei.first = NULL;
	TBnuclei.last = NULL;
	TBnuclei.nucleateFlag = NOTHING;

	SBnuclei.count = 0;
	SBnuclei.first = NULL;
	SBnuclei.last = NULL;
	SBnuclei.nucleateFlag = NOTHING;

	PSNnuclei.count = 0;
	PSNnuclei.first = NULL;
	PSNnuclei.last = NULL;
	PSNnuclei.nucleateFlag = NOTHING;
}

void boundariesPool::initializeBoundaries( void )
{
	if( owner )
	{
		offset = owner->grains.defGrainsCount;
		bucketSize = owner->grains.Nuc;
		boundaryBucket = (boundaryColumnP *) calloc( bucketSize, sizeof(boundaryColumnP) );

		if( !boundaryBucket ) owner->err->reportError( ERRTXT( "Cannot allocate more memory" ) );

		//memset( boundaryBucket, 0, sizeof(boundaryBucket) );

		for( int i = 0; i < bucketSize; i++ )	boundaryBucket[i] = NULL;
	}
	else boundaryBucket = NULL;
}

void boundariesPool::calculateFastMODF( void ) //###LB: Bin size for MODF 2.03°
{
	MPI_Barrier( MPI_COMM_WORLD );
	Real binSize = 2.03 * _PI_ / 180;

	memset( MODFCounts, 0, sizeof(MODFCounts) );
	memset( MODF_Cumulative, 0, sizeof(MODF_Cumulative) );

	long i;
	long sum = 0;

	for( i=0; i < bucketSize; i++ )
	{
		boundaryColumnP bin = boundaryBucket[i];

		if( !bin ) continue;

		long len = bin->len;
		
		for( long j = 0; j < len; j++ )
		{
			Real mis = bin->entry[j].disori;
			long idx = (long) ( mis / binSize );
			MODFCounts[idx] += bin->entry[j].count;
			sum += bin->entry[j].count;
		}
	}

	long allCounts[32] = {0};
	long totalCounts = 0;

	MPI_Allreduce( MODFCounts, allCounts, 32, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &(this->count), &totalCounts, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

	long suma = 0;

	for( i = 0; i < 32; i++ )
	{
		MODF[i] = (Real) allCounts[i] / (totalCounts * binSize);
		MODF_Cumulative[i] = ( i > 0 ) ? ( (Real) allCounts[i] / totalCounts ) + MODF_Cumulative[i-1] : ( (Real) allCounts[i] / totalCounts);
		suma += allCounts[i];
	}
	if( suma != totalCounts )
		owner->err->reportError(ERRTXT("suma != totalCounts"));

	QUICKASSERT( suma == totalCounts );
}

void boundariesPool::calculateMODF( void ) //###LB: Bin size for MODF 2.03°
{
	MPI_Barrier( MPI_COMM_WORLD );
	Real binSize = 2.03 * _PI_ / 180;

	memset( MODFCounts, 0, sizeof(MODFCounts) );
	memset( MODF_Cumulative, 0, sizeof(MODF_Cumulative) );

	for( cellsBoundaryP b = first; b; b = b->next )
	{
		long idx = (long) (b->disori / binSize);
		MODFCounts[idx]++;
	}
	long allCounts[32] = {0};
	long totalCounts = 0;

	MPI_Allreduce( MODFCounts, allCounts, 32, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &(this->count), &totalCounts, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

	long sum = 0;

	for( int i = 0; i < 32; i++ )
	{
		MODF[i] = (Real) allCounts[i] / (totalCounts * binSize);
		MODF_Cumulative[i] = ( i > 0 ) ? ( (Real) allCounts[i] / totalCounts ) + MODF_Cumulative[i-1] : ( (Real) allCounts[i] / totalCounts);
		sum += allCounts[i];
	}
}

void boundariesPool::addBoundaryPiece( long pUp, long pDown, Real disori )
{
	uint_fast64_t id = MYHASH( pUp, pDown );

	long pmax = MAX( pUp, pDown );
	
	long pos = pmax - offset;

	boundaryColumnP bin = boundaryBucket[pos];

	if( !bin ) //create the bin if it not exists
	{
		bin = (boundaryColumnP) malloc( sizeof(boundaryColumn) + STDLEN * sizeof(boundaryEntry) );
		      
		if( !bin ) owner->err->reportError( ERRTXT( "Cannot allocate more memory" ) );

		bin->len = 0;
		bin->maxLen = STDLEN;
		boundaryBucket[pos] = bin;
	}

	int i;
	for ( i = 0; i<bin->len; i++ )
	{
		if( bin->entry[i].id == id )
		{
			QUICKASSERT( SQR(bin->entry[i].disori - disori) < SQR(1e-3) );  //###LB: test; delete
			bin->entry[i].count++;
			return;
		}
	}

	//and account for this also for the whole OriBucket, this might make a re-malloc necessary
	if( i < bin->maxLen )
	{                                                                       
		bin->entry[i].id = id;                                            
		bin->entry[i].count = 1;
		bin->entry[i].disori = disori;
		bin->entry[i].gposDown = pDown;
		bin->entry[i].gposUp = pUp;
		bin->entry[i].first = NULL;
		bin->entry[i].last = NULL;
		bin->entry[i].tempCount = 0;
		bin->len++;
		gbCount++;
		return;
	}

	if( bin->len != bin->maxLen ) owner->err->reportError( ERRTXT("bin->len != bin->maxLen") );

	boundaryColumnP binlong = (boundaryColumnP) malloc( sizeof(boundaryColumn) + bin->maxLen * 2 * sizeof(boundaryEntry) );

	if( !binlong ) 
		owner->err->reportError( ERRTXT( "Cannot allocate more memory" ) );

	binlong->len = bin->maxLen;
	binlong->maxLen = 2 * bin->maxLen;

	for( i=0; i<bin->maxLen; i++)
	{
		binlong->entry[i] = bin->entry[i];
	}
	
	gbCount++;

	binlong->entry[bin->maxLen].id = id;                                            
	binlong->entry[bin->maxLen].count = 1;
	binlong->entry[bin->maxLen].disori = disori;
	binlong->entry[bin->maxLen].gposDown = pDown;
	binlong->entry[bin->maxLen].gposUp = pUp;
	binlong->entry[bin->maxLen].first = NULL;
	binlong->entry[bin->maxLen].last = NULL;
	binlong->entry[bin->maxLen].tempCount = 0;
	binlong->len++;

	boundaryBucket[pos] = binlong;
	free(bin);
	//re-malloc procedure if oribin is too small
}

boundariesPool::~boundariesPool( void )
{
	cellsBoundaryP bnext = NULL;

	for( long i = 0; i < gbCount; i++ )
	{
		for( cellsBoundaryP b = boundaryEntries[i]->first; b; b=bnext )
		{
			bnext = b->next;
			delete b;
		}
		boundaryEntries[i]->first = NULL;
	}
	for( long i = 0; i < bucketSize; i++ )		free( boundaryBucket[i] );
	free( boundaryEntries );
	free( boundaryBucket );
}

void boundariesPool::emptyBoundaries( void )
{
/*	cellsBoundaryP bnext = NULL;

	for( long i = 0; i < gbCount; i++ )
	{
		for( cellsBoundaryP b = boundaryEntries[i]->first; b; b=bnext )
		{
			bnext = b->next;
			delete b;
		}
		boundaryEntries[i]->first = NULL;
	}
	for( long i = 0; i < bucketSize; i++ )
	{
		free( boundaryBucket[i] );
		boundaryBucket[i] = NULL;
	}*/
}

cellsBoundary::cellsBoundary( boundariesPool * own, cellP up, cellP down )
{
	owner = own;
	upGrain.ix = up->ix;
	upGrain.iy = up->iy;
	upGrain.iz = up->iz;

	downGrain.ix = down->ix;
	downGrain.iy = down->iy;
	downGrain.iz = down->iz;

	if( owner->first )	owner->first->prev = this;

	this->next = owner->first;
	this->prev = NULL;
	owner->first = this;

	if( !(owner->last) )	owner->last = this;
	owner->count++;

	caHdlP masterOwner = owner->owner;

	long pu = masterOwner->posOfGrainOfCell( up );
	long pd = masterOwner->posOfGrainOfCell( down );

	grainP upGrain = masterOwner->grains.all[pu];
	grainP downGrain = masterOwner->grains.all[pd];

	disori = owner->calculateDisorientation( upGrain, downGrain );
	owner->addBoundaryPiece( pu, pd, disori );
}

cellsBoundary::~cellsBoundary( void )
{
}

Real boundariesPool::calculateDisorientation( grainP up, grainP down )
{
	long ip = up->oriIndex;
	long iq = down->oriIndex;

	oriRepresentationP p = up->owner->orientations[ip];
	oriRepresentationP q = down->owner->orientations[iq];

	Real qq[4] = { q->q0, q->q1, q->q2, q->q3 };
	Real pq[4] = { p->q0, p->q1, p->q2, p->q3 };

	Real misq[4];

	owner->misorientationQuaternionCubic( pq, qq, misq );

	Real theta = misq[0];

	QUICKASSERT( theta < 1.01 );

	if( theta > 1.0 )
		theta = (Real) (int) theta;

	return 2*acos(theta);
}

//###LB: What this function does is to take the 2D grain boundary array, which is sparsely filled and allocate it in a linear array for subsequent use of fast
//###LB: operations
void boundariesPool::sortGrainBoundaryCells( void )
{
	MPI_Barrier(MPI_COMM_WORLD);

	boundaryEntries = (boundaryEntryP *) calloc( gbCount, sizeof(boundaryEntryP) );
	if( !boundaryEntries ) owner->err->reportError(ERRTXT("Cannot allocate more memory"));

	cellsBoundaryP bnext = NULL;
	for( cellsBoundaryP b = first; b; b = bnext )
	{
		bnext = b->next;

		long ixu = b->upGrain.ix;
		long iyu = b->upGrain.iy;
		long izu = b->upGrain.iz;
		long ixd = b->downGrain.ix;
		long iyd = b->downGrain.iy;
		long izd = b->downGrain.iz;

		cellP * cup = owner->getCellPerP( ixu, iyu, izu );
		cellP * cdo = owner->getCellPerP( ixd, iyd, izd );

		//cellP cup = owner->getCellPer( ixu, iyu, izu );
		//cellP cdo = owner->getCellPer( ixd, iyd, izd );

		long pup = grainArrayPosition( *cup );
		long pdo = grainArrayPosition( *cdo );

		//long pup = owner->posOfGrainOfCell( cup );
		//long pdo = owner->posOfGrainOfCell( cdo );

		long bucketColumn = MAX( pup, pdo ) - offset;

		boundaryColumnP bin = boundaryBucket[bucketColumn];

		if( !bin )
			owner->err->reportError(ERRTXT("!bin"));
		
		long len = bin->len;

		for( long j = 0; j < len; j++ )
		{
			if( MAX( bin->entry[j].gposDown, bin->entry[j].gposUp ) != ( bucketColumn + offset ) )
				owner->err->reportError(ERRTXT("MAX( bin->entry[j].gposDown, bin->entry[j].gposUp ) != ( bucketColumn + offset )"));

			if( MYHASH(pup,pdo) == bin->entry[j].id )
			{
				if( bin->entry[j].first )	bin->entry[j].first->prev = b;

				b->next = bin->entry[j].first;
				b->prev = NULL;
				bin->entry[j].first = b;

				if( !(bin->entry[j].last) )	bin->entry[j].last = b;
				bin->entry[j].tempCount ++;
			}
		}

	}
	
	long k = 0;
	long i;

	long bucketStart = offset;
	
	for( i = bucketStart ; i < bucketSize; i++ )
	{
		boundaryColumnP bin = boundaryBucket[i];

		if( !bin ) continue;

		long len = bin->len;

		for( long j = 0; j < len; j++ )
		{
			boundaryEntries[k] = &(bin->entry[j]);

			if( bin->entry[j].count != bin->entry[j].tempCount )
				owner->err->reportError(ERRTXT("bin->entry[j].count != bin->entry[j].tempCount"));

			k++;
		}
	}

	if( gbCount != k )
		owner->err->reportError(ERRTXT("gbCount != k"));
}

void boundariesPool::recalculateDisorientation( long posga, long posgb ) //posa and posb were passed already with the offset
{
	long i;

	long minPos = MIN( posga, posgb );

	long bucketStart = minPos - offset;
	
	for( i = bucketStart ; i < bucketSize; i++ )
	{
		boundaryColumnP bin = boundaryBucket[i];

		long grainPos = i + offset;

		if( !bin ) continue;

		long len = bin->len;

		for( long j = 0; j < len; j++ )
		{
			long ga = bin->entry[j].gposDown;
			long gb = bin->entry[j].gposUp;

			if( ( ga != posga && ga != posgb ) && ( gb != posga && gb != posgb ) ) 
				continue;

			if(!(ga == posga || ga == posgb || gb == posga || gb == posgb)) 
				owner->err->reportError(ERRTXT("ga == posga || ga == posgb || gb == posga || gb == posgb"));

			if(!(grainPos == ga || grainPos == gb)) 
				owner->err->reportError(ERRTXT("grainPos == ga || grainPos == gb"));

			grainP upGrain = owner->grains.all[gb];
			grainP downGrain = owner->grains.all[ga];
			
			bin->entry[j].disori = calculateDisorientation( upGrain, downGrain );
			
		}
	}
}
