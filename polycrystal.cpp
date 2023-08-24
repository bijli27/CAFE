#include "polycrystal.h"
#include "caHdl.h"

polyCrystal::polyCrystal( void * own )
{
	self = own;
}

polyCrystal::polyCrystal( void )
{
}

polyCrystal::~polyCrystal(void)
{

}

////////////////////////////misorientation//////////////////////////////////////
Real polyCrystal::GetMobilityWeight( grainP growGrain, grainP shrinkGrain )
{
	caHdlP k = crystal();
	errorP err = k->err;

	long hashIndex = (((long)growGrain) & HASHCOLUMNMASK) >> HASHCOLUMNSHIFT;       //calculates the hash index from the address of the grain, i.e.the tab structure is for germs and seeds
	hashColumnP col = shrinkGrain->cp[hashIndex];									//0..15 selects one of the 16 column pointers (-> initialize the pointer: new grain) the shrinking grain stores which germs grow in it	if( col == NULL )
	{
		col = (hashColumnP) malloc(sizeof(hashColumn)+STDLEN*sizeof(hashTabEntry));   // Storage for columns

		if( !col ) err->reportError( ERRTXT( "Cannot allocate more memory" ) );

		col->len = 0;                                                                        //Initialisieren
		col->maxLen = STDLEN;                                                                //Specifying the max
		shrinkGrain->cp[hashIndex] = col;                                                    //cPointer grain on column
	}

	long i;
	for ( i = 0; i<col->len; i++ )                                            //if mobility already exists
	{
		if( col->entry[i].identity == growGrain )                       //identity == address grain//identity==Adresse Korn
		{
			return col->entry[i].P;                                //Return of the table value//Rueckgabe des Tabellen-Wertes
		}
	}

	Real P = CalcMobilityWeight( growGrain, shrinkGrain);       //if mobility is not available, simply calculate//falls Mobilitaet nicht vorhanden, einfach berechnen

	if( i < col->maxLen )
	{                                                                       //if there is still space in the table//falls in der Tabelle noch Platz ist
		col->entry[i].P = P;                                            //Write the calculated m in// berechnetes m reinschreiben
		col->entry[i].identity = growGrain;
		col->len++;
		return P;                                                       //return m
	}

	ASSERT(col->len==col->maxLen);                                          //Extend table if space is insufficient // Tabelle verl�ngern, falls platz nicht ausreichend

	hashColumnP collong = (hashColumnP) malloc(sizeof(hashColumn)+col->maxLen*2*sizeof(hashTabEntry));   //Allocate memory
	if( !collong ) err->reportError( ERRTXT( "Cannot allocate more memory" ) );

	collong->len = col->maxLen;   //Endwert der einen gleich Anfangswert der neuen Spalte
	collong->maxLen = col->maxLen*2;  //neues Maximum

	for( i=0; i<col->maxLen; i++)              //Transfer values //Werte uebertragen
	{
		collong->entry[i] = col->entry[i];
	}
	collong->entry[col->maxLen].P = P;            //m Write to new memory locations // m in neue Speicherpl�tze schreiben
	collong->entry[col->maxLen].identity = growGrain;    //also the grain addresse //benso die Kornadresse
	collong->len++;

	shrinkGrain->cp[hashIndex] = collong;		//new column pointer//neuer Spaltenzeiger
	free(col);									//free old column//alte Spalte freigeben
	return (P);									//return m//m zurueckgeben
}

Real polyCrystal::CalcMobilityWeight( grainP growGrain, grainP shrinkGrain )		//calculate mobility//Mobilitaet berechnen
{
	Real weightedMob = 0;

	Real misq[4] = {0};
	Real maxDev40_111 = 0.174533;	// Has to be user defined; here 10�. 
	 								//It's different from the previously defined TwistAngleDeviation
	
	//caHdlP own = this->owner;		//###LB: deviation can be taken from owner and be defined in core8.uds


	Real _sqrt3 = 1 / sqrt(3.0); 
	Real oneNinth = 1.0 / 9.0;
	Real m40_111[4] = { cos( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ) };

	int ip = growGrain->oriIndex;
	int iq = shrinkGrain->oriIndex;

	oriRepresentationP p = growGrain->owner->orientations[ip];
	oriRepresentationP q = shrinkGrain->owner->orientations[iq];

	Real qq[4] = { q->q0, q->q1, q->q2, q->q3 };
	Real pq[4] = { p->q0, p->q1, p->q2, p->q3 };

	misorientationQuaternionCubic( pq, qq, misq );

	Real theta = misq[0];

	QUICKASSERT( theta < 1.01 );

	if( theta > 1.0 )
		theta = (Real) (int) theta;

	theta = 2*acos(theta);

	QUICKASSERT( theta <= 1.099 );

	if( theta <= 0.2618 ) //15�
		weightedMob = -1.0;
	else
	{
		weightedMob = 0.0;

		Real dev_40_111 = misorientationCubicQxQ( misq[0], misq[1], misq[2], misq[3], m40_111[0], m40_111[1], m40_111[2], m40_111[3] );
		if( dev_40_111 <  maxDev40_111 )
			weightedMob = SQR( cos( 0.5 * _PI_ * dev_40_111 / maxDev40_111 ) );
	}
	return weightedMob;
}















