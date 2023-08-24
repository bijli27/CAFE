#include "mymath.h"
#include "applic.h"
#include "random.h"

//###LB: Unless indicated, all functions in radians.

mathMethods::mathMethods( void )
{
	r.init( -4256 );
}

mathMethods::~mathMethods( void )
{
}

char mathMethods::isPrime( long n )					//Prime function
{
	long n_root = (long) ( sqrt((double)n) + 0.5 );

	for( long i = 2; i <= n_root; i++ )
		if( n % i == 0 ) return 0;

	return 1;
}

void mathMethods::factorizeIn3( long n, long * factors )		//Factor
{
	long f[3]={0};

	if( isPrime( n ) )
	{
		factors[0]=n;
		factors[1]=1;
		factors[2]=1;
		return;
	}

	long nhalf = (long) (0.5*n);
	long *Fac = (long *) calloc( nhalf, sizeof( long ) );

	if( !Fac ) exitus("ERROR: Cannot allocate more memory :(");

	int count=0;

	for( long i = nhalf; i > 0; i-- )
	{
		if( n % i==0 )
		{
			Fac[count]=i;
			count++;
		}
	}

	long ft = count;
	count=0;
	long sum = 3*n;

	for( long i=0;i<ft;i++ )
		for( long j=0;j<ft;j++ )
		{
			long p = Fac[i]*Fac[i]*Fac[j];
			long s = Fac[i]+Fac[i]+Fac[j];

			if( p == n && s < sum )
			{
				f[0]=Fac[i];
				f[1]=Fac[i];
				f[2]=Fac[j];
				sum = s;
			}
		}

		for( long i=0;i<ft;i++ )
			for( long j=i+1;j<ft;j++ )
				for( long k=j+1;k<ft;k++ )
				{
					long p = Fac[i]*Fac[j]*Fac[k];
					long s = Fac[i]+Fac[j]+Fac[k];
					if( p == n && s < sum )
					{
						f[0]=Fac[i];
						f[1]=Fac[j];
						f[2]=Fac[k];
						sum = s;
					}
				}

		sortInt( f,3 );
		factors[0]=f[2];
		factors[1]=f[1];
		factors[2]=f[0];
		free(Fac);
}

double mathMethods::areaPolygon( int *poly, int counter )
{
        double sum=0;

        for( int i=0;i<(2*counter-2);i+=2 )
        {
                sum += (poly[i]+poly[i+2])*(poly[i+3]-poly[i+1]);
        }
        if( sum < 0 ) sum *= -1;
        return 0.5 * sum;
}

int mathMethods::pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
        int i, j, c = 0;
        for (i = 0, j = nvert-1; i < nvert; j = i++)
        {
                if ( ((verty[i]>testy) != (verty[j]>testy)) && (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
                        c = !c;
        }
        return c;
}

void mathMethods::bubbleSort ( Real arr [ ], int size ) // Sort components of a quaternion
 { 
    int last = size - 2; 
    int isChanged = 1; 

    while ( last >= 0 && isChanged ) 
    { 
            isChanged = 0; 
            for ( int k = 0; k <= last; k++ ) 
                if ( arr[k] > arr[k+1] ) 
                { 
                    swap ( arr[k], arr[k+1] ); 
                    isChanged = 1; 
                } 
            last--; 
    } 
 }

void mathMethods::bubbleSortIndex( Real arr [], long idx[], int size ) // Sort components of a quaternion
 { 
    int last = size - 2; 
    int isChanged = 1; 

    while ( last >= 0 && isChanged ) 
    {
            isChanged = 0; 
            for ( int k = 0; k <= last; k++ ) 
                if ( arr[k] > arr[k+1] ) 
                { 
                    swap( arr[k], arr[k+1] );
					swapInt( idx[k], idx[k+1] );
                    isChanged = 1; 
                } 
            last--; 
    } 
 }

void mathMethods::sortInt( long arr [ ], long size ) // Sort integers
 { 
    long last = size - 2; 
    long isChanged = 1; 

    while ( last >= 0 && isChanged ) 
    { 
            isChanged = 0; 
            for ( long k = 0; k <= last; k++ ) 
                if ( arr[k] > arr[k+1] ) 
                { 
                    swapInt ( arr[k], arr[k+1] ); 
                    isChanged = 1; 
                } 
            last--; 
    } 
 }

void mathMethods::swapInt( long& x, long& y ) //Required for the sorting
{
    long temp;
    temp = x;
    x = y;
    y = temp;
}
 
void mathMethods::swap( Real& x, Real& y ) //Required for the sorting
{
    Real temp;
    temp = x;
    x = y;
    y = temp;
}

void mathMethods::counting_sort(int *A, int size, int range) 	//Sorting counted numbers
{
	int *B = new int[size];
	int *C = new int[range + 1];

	for(int i = 0; i <= range; ++i)
		C[i] = 0;

	for(int i = 0; i < size; ++i)
		++C[A[i]];

	for(int i = 1; i <= range; ++i)
		C[i] += C[i - 1];

	for (int i = size - 1; i >= 0; --i)
	{
		B[C[A[i]] - 1] = A[i];
		--C[A[i]];
	}

	for (int i = 0; i < size; ++i)
		A[i] = B[i];

  delete[] B;
  delete[] C;
}

////////////////////////////- Convert Weightpercent to Atompercent/100 -////////////////////////////
Real mathMethods::convertConc(Real conc,int elemTyp)
{
/*     int k,kMax;        //Al elemTyp= Cr-1 Cu-2 Fe-3 Mg-4 Mn-5 Si-6 Ti-7 Zn-8 (elemTyp=8)
     kMax = 9;         //number of chemical elements considered in ClaNG
     Real ElemC[ORICLASSMAX];
     Real SumC=0.0;
     for (k = 0 ;k<kMax ;k++ )
     {
            if(k == (elemTyp-1)) ElemC[k] = conc/100.0;
            else ElemC[k]=parm.ElemC[k]/100.0;
            if(k == 0) ElemC[k]= 1 - 1/100.0 * (parm.ElemC[k+1]+parm.ElemC[k+2]+parm.ElemC[k+3]+parm.ElemC[k+4]+parm.ElemC[k+5]+parm.ElemC[k+6]+parm.ElemC[k+7]+parm.ElemC[k+8]);
//            printf("%f\n", ElemC[k]);
            SumC+= ElemC[k]/parm.ElemMolM[k];
     }
     QUICKASSERT(SumC > 0.0);
     conc = ElemC[elemTyp-1]/(parm.ElemMolM[elemTyp-1]*SumC);

     return conc;*/
	return 0;
}

void mathMethods::K_S_Test( double * data1, unsigned long n1, double * data2, unsigned long n2, double * d, double * prob )
{
	unsigned long j1 = 1;
	unsigned long j2 = 1;
	double d1, d2, dt, en1, en2, en;
	double fn1 = 0.0;
	double fn2 = 0.0;

	sort( n1, data1 );
	sort( n2, data2 );

	en1 = n1;
	en2 = n2;
	*d = 0.0;

	while( j1 <= n1 && j2 <= n2 ) {
		if( (d1=data1[j1]) <= (d2=data2[j2]) ) fn1 = j1++ / en1;
		if( d2 <= d1 ) fn2 = j2++/en2;
		if( (dt=fabs(fn2-fn1)) > *d ) *d = dt;
	}
	en = sqrt(en1*en2/(en1+en2));
	*prob = probks( (en+0.12+0.11/en)*(*d) );
}

double mathMethods::probks( double alam )
{
	int j;
	double a2, term;
	double fac = 2.0;
	double sum = 0.0;
	double termbf = 0.0;

	a2 = -2.0 * SQR(alam);
	for( j=1; j<=100; j++ )
	{
		term = fac * exp( a2 * SQR( j ) );
		sum += term;
		if( fabs(term) <= EPS1 * termbf || fabs(term) <= EPS2 * sum ) return sum;
		fac = -fac;
		termbf = fabs( fac );
	}
	return 1.0;
}

void mathMethods::sort(int n, double *ra)
{
	int l, j, ir, i;
	float rra;

	l = (n >> 1) + 1;
	ir = n;

	for ( ; ; )
	{
		if (l > 1)					/* still in hiring phase */
			rra = ra[--l];
		else						/* in retirement-and-promotion phase */
		{
			rra = ra[ir];           /* clear a space at end of array */
			ra[ir]=ra[1];			/* retire the top of the heap into it */
			if (--ir == 1) 			/* done with last promotion */
			{
				ra[1] = rra;
				return;
			}						/* end if */
		}							/* end else */
		i = l;						/* whether we are in the hiring phase */
		j = l << 1;					/* or promotion phase, we here set up */
		while ( j <= ir )
		{
			if ( (j < ir) && (ra[j] < ra[j + 1]) )
				++j;				/* compare to the better underling */
				if ( rra < ra[j] )	/* demote rra */
				{
					ra[i] = ra[j];
					j += (i = j);
				}
				else
					j = ir + 1;		/* this is rra's level; set j to */
		}                           /* terminate the sift-down */
		ra[i] = rra;				/* put rra into its slot */
	}
}

double mathMethods::misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 )
{
        int i;

        Real p[4] = {q01,q11,q21,q31};
	Real q[4] = {q02,q12,q22,q32};

	Real qm1[4];    //Inverse of quaternion q

	for(i=0;i<4;i++)       //Inverting unit quaternion
        {
	        qm1[i]=q[i];
                if( i>0 ) qm1[i]*=-1;
        }

	Real r[4];     //Resulting quaternion, rotation of the two previous quaternions pq-1

        multiplyQuaternions( p, qm1, r );

        //Now, we have to determine the smallest angle.

	Real r0[6][4];    //There are 12 possible angles

        //Note: The notation r0 is due to the definition of the quaternion which lie
        //in the fundamental zone, this vector possesses the smallest angle, in such a way
        //that r0 is actually the scalar part of this quaternion

	double a,b,c,d;
	Real rt3=sqrt(3.0);

	a=r[0]; b=r[1]; c=r[2]; d=r[3];

	Real fac=0.70710678;

	r0[0][0]=(r[0]+r[1])*fac; r0[0][1]=(r[0]-r[1])*fac; r0[0][2]=(r[2]+r[3])*fac; r0[0][3]=(r[2]-r[3])*fac;
	r0[1][0]=(r[0]+r[2])*fac; r0[1][1]=(r[0]-r[2])*fac; r0[1][2]=(r[1]+r[3])*fac; r0[1][3]=(r[1]-r[3])*fac;
	r0[2][0]=(r[0]+r[3])*fac; r0[2][1]=(r[0]-r[3])*fac; r0[2][2]=(r[1]+r[2])*fac; r0[2][3]=(r[1]-r[2])*fac;
	r0[3][0]=(r[0]+r[1]+r[2]+r[3])*0.5; r0[3][1]=(r[0]+r[1]-r[2]-r[3])*0.5; r0[3][2]=(r[0]-r[1]+r[2]-r[3])*0.5; r0[3][3]=(r[0]-r[1]-r[2]+r[3])*0.5;
	r0[4][0]=(r[0]+r[1]+r[2]-r[3])*0.5; r0[4][1]=(r[0]+r[1]-r[2]+r[3])*0.5; r0[4][2]=(r[0]-r[1]+r[2]+r[3])*0.5; r0[4][3]=(r[0]-r[1]-r[2]-r[3])*0.5;
	r0[5][0]=r[0];r0[5][1]=r[1];r0[5][2]=r[2];r0[5][3]=r[3];


	Real omega=0.0;

	for(i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega )
				omega=fabs(r0[i][j]);

	QUICKASSERT( omega < 1.01 );

	if( omega > 1.0 )
		omega = (Real) (int) omega;

	omega=2*acos(omega);
	QUICKASSERT( omega <= 1.099 );
	return omega;
}

void mathMethods::randomMisorientationAxisConsidered(  double * qref, double * qr, double maxDev  )
{
        Real theta = cos( 0.5 * _PI_ );

        double q[4] = {0};

        double dev = cos(0.5 * maxDev);

        double refNorm = sqrt( SQR(qref[0]) + SQR(qref[1]) + SQR(qref[2]) + SQR(qref[3]) );

        qref[0] /= refNorm;
        qref[1] /= refNorm;
        qref[2] /= refNorm;
        qref[3] /= refNorm;

        while( theta < dev  )
        {
                double s = r.parkMiller();
                double sigma1 = sqrt(1-s);
                double sigma2 = sqrt(s);
                double theta1 = 2 * _PI_ * r.parkMiller();
                double theta2 = 2 * _PI_ * r.parkMiller();

                q[0]=sigma2*cos(theta2);
                q[1]=sigma1*sin(theta1);
                q[2]=sigma1*cos(theta1);
                q[3]=sigma2*sin(theta2);

				
                double norm = sqrt(SQR(q[0])+SQR(q[1])+SQR(q[2])+SQR(q[3]));

                q[0] /= norm;
                q[1] /= norm;
                q[2] /= norm;
                q[3] /= norm;

				long idx[4]={0,1,2,3};
				int signQ[4] = {SIGN(q[0]),SIGN(q[1]),SIGN(q[2]),SIGN(q[3])};

				for( int i = 0; i<4; i++ ) q[i] *= signQ[idx[i]];

                bubbleSortIndex( q, idx, 4 );

				for( int i = 0; i<4; i++ ) q[i]*=signQ[idx[i]];

                theta = q[3]*qref[0] + q[2]*qref[1] + q[1]*qref[2] + q[0]*qref[3];
        }
        qr[0]=q[3];
        qr[1]=q[2];
        qr[2]=q[1];
        qr[3]=q[0];

}

double mathMethods::angleBetweenQuaternions( double * q, double * p )
{
        double _qNorm = 1 / sqrt(SQR(q[1])+SQR(q[2])+SQR(q[3])+SQR(q[0]));
        double _pNorm = 1 / sqrt(SQR(p[1])+SQR(p[2])+SQR(p[3])+SQR(p[0]));
        return acos( _qNorm * _pNorm * ( q[0]*p[0] + q[1]*p[1] + q[2]*p[2] + q[3]*p[3] ) );
}

void mathMethods::randomOrientation( double * result )
{
	double q[4]={0,0,0,0};

	double s = r.parkMiller();
	double sigma1 = sqrt(1-s);
	double sigma2 = sqrt(s);
	double theta1 = 2 * _PI_ * r.parkMiller();
	double theta2 = 2 * _PI_ * r.parkMiller();

	q[0]=sigma2*cos(theta2);
	q[1]=sigma1*sin(theta1);
	q[2]=sigma1*cos(theta1);
	q[3]=sigma2*sin(theta2);

	long idx[4]={0,1,2,3};
	int signQ[4] = {SIGN(q[0]),SIGN(q[1]),SIGN(q[2]),SIGN(q[3])};

	for( int i = 0; i<4; i++ ) q[i] *= signQ[idx[i]];

	bubbleSortIndex( q, idx, 4 );

	for( int i = 0; i<4; i++ ) q[i]*=signQ[idx[i]];

	Real qr[4];

	qr[0]=q[3];
	qr[1]=q[2];
	qr[2]=q[1];
	qr[3]=q[0];

	Real angles[3] = {0};

	quaternion2Euler( qr, angles );

	result[0] = angles[0];
	result[1] = angles[1];
	result[2] = angles[2];
}

void mathMethods::randomMisorientation( double theta, double* qr  ) // theta in radians
{
        double q[4]={0,0,0,0};
        double qcrit = cos(0.5*theta);

        while( q[3]<qcrit )
        {
                double s = r.parkMiller();
                double sigma1 = sqrt(1-s);
                double sigma2 = sqrt(s);
                double theta1 = 2 * _PI_ * r.parkMiller();
                double theta2 = 2 * _PI_ * r.parkMiller();

                q[0]=sigma2*cos(theta2);
                q[1]=sigma1*sin(theta1);
                q[2]=sigma1*cos(theta1);
                q[3]=sigma2*sin(theta2);

                long idx[4]={0,1,2,3};
				int signQ[4] = {SIGN(q[0]),SIGN(q[1]),SIGN(q[2]),SIGN(q[3])};
				
				for( int i = 0; i<4; i++ ) q[i] *= signQ[idx[i]];

                bubbleSortIndex( q, idx, 4 );

				for( int i = 0; i<4; i++ ) q[i]*=signQ[idx[i]];
        }
        qr[0]=q[3];
        qr[1]=q[2];
        qr[2]=q[1];
        qr[3]=q[0];
}

void mathMethods::multiplyQuaternions( double *q, double* p, double* r )
{
        r[0]=q[0]*p[0]-q[1]*p[1]-q[2]*p[2]-q[3]*p[3];
        r[1]=q[1]*p[0]+q[0]*p[1]-q[3]*p[2]+q[2]*p[3];
        r[2]=q[2]*p[0]+q[3]*p[1]+q[0]*p[2]-q[1]*p[3];
        r[3]=q[3]*p[0]-q[2]*p[1]+q[1]*p[2]+q[0]*p[3];
}

void mathMethods::misorientationQuaternionCubic( double* p, double* q, double* quat  )
{

	Real qm1[4];
	int i;

	for(i=0;i<4;i++)                 //inverting unit quaternion q
        {				//Copy quaternion; not really necessary
		qm1[i]=q[i];
                if( i>0 ) qm1[i]*=-1;
        }

	Real r1[4];

        multiplyQuaternions( p, qm1, r1 );

	Real sqrt2=1/sqrt(2.0);
	Real rot[6][4];
	Real a, b, c ,d;

	//The six fundamental quaternions describing the misorientation

	a=r1[0]; b=r1[1]; c=r1[2]; d=r1[3];

	rot[0][0] = a; rot[0][1] = b; rot[0][2] = c; rot[0][3] = d;

	rot[1][0] = sqrt2 * ( a + b ); rot[1][1] = sqrt2 * ( a - b ); rot[1][2] = sqrt2 * ( c + d ); rot[1][3] = sqrt2 * ( c - d );

	rot[2][0] = sqrt2 * ( a + c ); rot[2][1] = sqrt2 * ( a - c ); rot[2][2] = sqrt2 * ( b + d ); rot[2][3] = sqrt2 * ( b - d );

	rot[3][0] = sqrt2 * ( a + d ); rot[3][1] = sqrt2 * ( a - d ); rot[3][2] = sqrt2 * ( b + c ); rot[3][3] = sqrt2 * ( b - c );

	rot[4][0] = 0.5 * ( a + b + c + d ); rot[4][1] = 0.5 * ( a + b - c - d ); rot[4][2] = 0.5 * ( a - b + c - d ); rot[4][3] = 0.5 * ( a - b - c + d );

	rot[5][0] = 0.5 * ( a + b + c - d ); rot[5][1] = 0.5 * ( a + b - c + d ); rot[5][2] = 0.5 * ( a - b + c + d ); rot[5][3]= 0.5 * ( a - b - c - d );

	Real rq[4];

	int mi;
	Real max=0.0;
	int j=0;

	for( i=0;i<6;i++ )						//Determing the quaternion with the maximal component and the component itself
		for( j=0;j<4;j++ )
		{
			if( fabs(rot[i][j]) > max )
			{
				max=fabs(rot[i][j]);
				mi=i;
				//mj=j;
			}
		}

	rq[0] = fabs( rot[mi][0] );					//Desorientation requires all components positive
	rq[1] = fabs( rot[mi][1] );
	rq[2] = fabs( rot[mi][2] );
	rq[3] = fabs( rot[mi][3] );

	bubbleSort( rq,4 );						//Sorting into ascending order, because a desorientation in the SST
									//requires a quaternion with q0>q1>q2>q3 which represents a minimal 
	quat[0] = rq[3];						//rotation angle and an axis fulfilling h>k>l
	quat[1] = rq[2];
	quat[2] = rq[1];
	quat[3] = rq[0];

}

double mathMethods::misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 )
{
	int i;
	Real oria[3] = { pa1, Pa, pa2 };
	Real orib[3] = { pb1, Pb, pb2 };

	Real p[4];
	Real q[4];

	euler2quaternion( oria, p );
	euler2quaternion( orib, q );

	Real qm1[4];    //Inverse of quaternion q

	for(i=0;i<4;i++)               //Inverting unit quaternion
        {
	        qm1[i]=q[i];
                if( i>0 ) qm1[i]*=-1;
        }

	Real r[4];     //Resulting quaternion, rotation of the two previous quaternions pq-1

        multiplyQuaternions( p, qm1, r );

        //Now, we have to determine the smallest angle.

	Real r0[6][4];    //There are 12 possible angles

        //Note: The notation r0 is due to the definition of the quaternion which lie
        //in the fundamental zone, this vector possesses the smallest angle, in such a way
        //that r0 is actually the scalar part of this quaternion

	double a,b,c,d;
	Real rt3=sqrt(3.0);

	a=r[0]; b=r[1]; c=r[2]; d=r[3];

	Real fac=0.70710678;

	r0[0][0]=(r[0]+r[1])*fac; r0[0][1]=(r[0]-r[1])*fac; r0[0][2]=(r[2]+r[3])*fac; r0[0][3]=(r[2]-r[3])*fac;
	r0[1][0]=(r[0]+r[2])*fac; r0[1][1]=(r[0]-r[2])*fac; r0[1][2]=(r[1]+r[3])*fac; r0[1][3]=(r[1]-r[3])*fac;
	r0[2][0]=(r[0]+r[3])*fac; r0[2][1]=(r[0]-r[3])*fac; r0[2][2]=(r[1]+r[2])*fac; r0[2][3]=(r[1]-r[2])*fac;
	r0[3][0]=(r[0]+r[1]+r[2]+r[3])*0.5; r0[3][1]=(r[0]+r[1]-r[2]-r[3])*0.5; r0[3][2]=(r[0]-r[1]+r[2]-r[3])*0.5; r0[3][3]=(r[0]-r[1]-r[2]+r[3])*0.5;
	r0[4][0]=(r[0]+r[1]+r[2]-r[3])*0.5; r0[4][1]=(r[0]+r[1]-r[2]+r[3])*0.5; r0[4][2]=(r[0]-r[1]+r[2]+r[3])*0.5; r0[4][3]=(r[0]-r[1]-r[2]-r[3])*0.5;
	r0[5][0]=r[0];r0[5][1]=r[1];r0[5][2]=r[2];r0[5][3]=r[3];

	Real omega=0.0;

	for(i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega )
				omega=fabs(r0[i][j]);

	QUICKASSERT( omega < 1.01 );

	if( omega > 1.0 )
		omega = (Real) (int) omega;

	omega=2*acos(omega);
	QUICKASSERT( omega <= 1.099 );
	return omega;
}

void mathMethods::newOrientationFromReferenceFixedAngularCone(double * oriOri, double maxDev,double angle, double u, double v, double w,double * newOri)
{
        Real qr[4];
        Real ori[4], qideal[4];

        Real _norm = 1.0 / sqrt(SQR(u)+SQR(v)+SQR(w));

        euler2quaternion( oriOri, qideal );
        Real qref[4] = { cos( 0.5* angle ), u * _norm * sin( 0.5 * angle ), v * _norm * sin( 0.5 * angle ), w * _norm * sin( 0.5 * angle ) };

        randomMisorientationAxisConsidered(  qref, qr, maxDev  );
        multiplyQuaternions( qr,qideal,ori );

        Real euler[3];
        quaternion2Euler( ori, euler );

        newOri[0] = euler[0];
        newOri[1] = euler[1];
        newOri[2] = euler[2];
}

void mathMethods::rotateOrientation( double *oriOri, double angle, double u, double v, double w, double *newOri ) //Angle in radians
{
        Real qr[4];
        Real ori[4], qideal[4];
        Real _norm = 1.0 / sqrt(SQR(u)+SQR(v)+SQR(w));

        euler2quaternion( oriOri, qideal );

        Real qref[4] = { cos( 0.5* angle ), u * _norm * sin( 0.5 * angle ), v * _norm * sin( 0.5 * angle ), w * _norm * sin( 0.5 * angle ) };

        multiplyQuaternions( qref, qideal, ori );

        double euler[3] = {0};

        quaternion2Euler( ori, euler );

        newOri[0] = euler[0];
        newOri[1] = euler[1];
        newOri[2] = euler[2];
}

void mathMethods::euler2quaternion( double * euler, double * q )
{
	double p1 = euler[0];
	double t  = euler[1];
	double p2 = euler[2];

	Real co1 = cos(t/2);
	Real s1 = sin(t/2);

	Real p[4] = {co1*cos((p1+p2)/2),s1*cos((p1-p2)/2),s1*sin((p1-p2)/2),co1*sin((p1+p2)/2)};

	double test[4]={0};

	quaternion2Euler( p,test);

	q[0] = p[0];
	q[1] = p[1];
	q[2] = p[2];
	q[3] = p[3];
}

void mathMethods::linearRegression(Real* x, Real* y, long size, Real* slope, Real* intercept, Real* dev)
{
	Real x2sum = 0;
	Real xysum = 0;
	Real xsum = 0;
	Real ysum = 0;

	// y = mx + c0

	for (int i = 0; i < size; i++)
	{
		x2sum += SQR(x[i]);
		xysum += (x[i] * y[i]);
		xsum += x[i];
		ysum += y[i];
	}
	Real num = (size * xysum - xsum * ysum);
	Real den = (size * x2sum - SQR(xsum));

	Real m =  1e30;

	if (den != 0) m = num / den;

	Real c0 = (ysum - m * xsum) / size;

	Real yavg = ysum / size;
	Real rss = 0;
	Real tss = 0;

	for (int i = 0; i < size; i++)
	{
		rss += SQR(y[i] - (m * x[i] + c0));
		tss += SQR(y[i] - yavg);
	}

	*slope = m;
	*intercept = c0;
	*dev = (1 - (rss / tss));
}

void mathMethods::quaternion2Euler( double * quat, double * euler )
{
	Real q0 = quat[0];
	Real q1 = quat[1];
	Real q2 = quat[2];
	Real q3 = quat[3];
	Real PHI, sP, phi1, phi2;

	Real cosPHI = SQR(q3)-SQR(q2)-SQR(q1)+SQR(q0);

	Real y0 = 2*q1*q3-2*q0*q2;
	Real x0 = 2*q2*q3+2*q0*q1;
	Real y1 = 2*q1*q3+2*q0*q2;
	Real x1 = -2*q2*q3+2*q0*q1;

	if( cosPHI > 1. ) cosPHI = 1.;

	if( SQR( 1. - cosPHI ) <= 1e-20 )
		PHI = 0.;
	else
		PHI = acos( cosPHI );

	sP = sin(PHI);

	if( sP != 0 )
	{
		phi2 = atan2( y0 / sP, x0 / sP );
		phi1 = atan2( y1 / sP, x1 / sP );
	}else
	{
		phi1 = atan2( (2*q1*q2+2*q0*q3),SQR(q0)+SQR(q1)-SQR(q2)-SQR(q3) );
		phi2 = 0.;
	}

	euler[0]=phi1;
	euler[1]=PHI;
	euler[2]=phi2;
}

void mathMethods::newOrientationFromReference( double *oriOri, double deviation, double *newOri )
{
	Real qr[4];
	Real ori[4], qideal[4];

	euler2quaternion( oriOri, qideal );

	randomMisorientation( deviation, qr  );
	multiplyQuaternions( qr,qideal,ori );

	Real euler[3];

	quaternion2Euler( ori, euler );

	newOri[0] = euler[0];
	newOri[1] = euler[1];
	newOri[2] = euler[2];
}

void mathMethods::devtorefEuler2RGB(double* bunge, double* ideal, double maxDev, unsigned char* rgb)
{
	//blue channel stretch from 0.0 to maxDev in radians, all other orientations white
	double dev = misorientationCubic(bunge[0], bunge[1], bunge[2], ideal[0], ideal[1], ideal[2]);

	if (dev < 0.0 || dev > MAXIMUM_MISORI_FCC) { //black as error mark
		rgb[0] = 0;
		rgb[1] = 0;
		rgb[2] = 0;
		return;
	}

	//valid and in range
	if (dev <= maxDev) { //more than maxDev_ref larger than maxDev are white
		double colorscale = 255.0 * pow((dev / maxDev), 1.0);
		int col = colorscale;

		rgb[0] = col;
		rgb[1] = col;
		rgb[2] = 255;
		return;
	}
	//otherwise white
	rgb[0] = RGBRANGE;
	rgb[1] = RGBRANGE;
	rgb[2] = RGBRANGE;
}

void mathMethods::QuatOnVector3D(double* q, double* v, double* r)
{
	//in accordance with Spie�2009 and Morawiec Pospiech 1989
	double a0 = q[0];
	double a1 = q[1];
	double a2 = q[2];
	double a3 = q[3];

	//get rotation matrix components
	double r11 = SQR(a0) + SQR(a1) - SQR(a2) - SQR(a3);
	double r12 = 2 * (a1 * a2 - a0 * a3);
	double r13 = 2 * (a1 * a3 + a0 * a2);

	double r21 = 2 * (a1 * a2 + a0 * a3);
	double r22 = SQR(a0) - SQR(a1) + SQR(a2) - SQR(a3);
	double r23 = 2 * (a2 * a3 - a0 * a1);

	double r31 = 2 * (a1 * a3 - a0 * a2);
	double r32 = 2 * (a2 * a3 + a0 * a1);
	double r33 = SQR(a0) - SQR(a1) - SQR(a2) + SQR(a3);

	//x,y,z vector3d v
	r[0] = (r11 * v[0]) + (r12 * v[1]) + (r13 * v[2]);
	r[1] = (r21 * v[0]) + (r22 * v[1]) + (r23 * v[2]);
	r[2] = (r31 * v[0]) + (r32 * v[1]) + (r33 * v[2]);
}

void mathMethods::project2fundamentalregion_ipfz(double* qtest, double* xy)
{
	//MK::ok, define crystallographic m-3m symmetry operators in quaternion representation
	double qsymm[SYMMETRIES_IN_FCC][4];

	double _sqrt2 = 1.0 / pow(2.0, 0.5);
	double half = 0.5;

	//<100>90*degrees four-fold symmetries
	qsymm[0][0] = 1.0;				qsymm[0][1] = 0.0;				qsymm[0][2] = 0.0;				qsymm[0][3] = 0.0;			/*identity*/
	qsymm[1][0] = _sqrt2;			qsymm[1][1] = _sqrt2;			qsymm[1][2] = 0.0;				qsymm[1][3] = 0.0;			/*[100]90*/
	qsymm[2][0] = 0.0;				qsymm[2][1] = 1.0;				qsymm[2][2] = 0.0;				qsymm[2][3] = 0.0;			/*[100]180*/
	qsymm[3][0] = -1.0 * _sqrt2;	qsymm[3][1] = _sqrt2;			qsymm[3][2] = 0.0;				qsymm[3][3] = 0.0;			/*[100]270*/
	qsymm[4][0] = _sqrt2;			qsymm[4][1] = 0.0;				qsymm[4][2] = _sqrt2;			qsymm[4][3] = 0.0;			/*[010]90*/ //MTex -sq2 0 -sq2 0
	qsymm[5][0] = 0.0;				qsymm[5][1] = 0.0;				qsymm[5][2] = 1.0;				qsymm[5][3] = 0.0;			/*[010]180*/ //##MK::be careful, this is consistent with Pomana it is inconsistent with Mtex stating 0,0,-1,0 but q = -q are the same quaternions
	qsymm[6][0] = -1.0 * _sqrt2;	qsymm[6][1] = 0.0;				qsymm[6][2] = _sqrt2;			qsymm[6][3] = 0.0;			/*[010]270*/
	qsymm[7][0] = _sqrt2;			qsymm[7][1] = 0.0;				qsymm[7][2] = 0.0;				qsymm[7][3] = _sqrt2;		/*[001]90*/
	qsymm[8][0] = 0.0;				qsymm[8][1] = 0.0;				qsymm[8][2] = 0.0;				qsymm[8][3] = 1.0;			/*[001]180*/
	qsymm[9][0] = -1.0 * _sqrt2;	qsymm[9][1] = 0.0;				qsymm[9][2] = 0.0;				qsymm[9][3] = _sqrt2;		/*[001]270*/

	//<110>180*degrees two-fold symmetries
	qsymm[10][0] = 0.0;				qsymm[10][1] = _sqrt2;			qsymm[10][2] = _sqrt2;			qsymm[10][3] = 0.0;			/*[110]180*/
	qsymm[11][0] = 0.0;				qsymm[11][1] = _sqrt2;			qsymm[11][2] = -1.0 * _sqrt2;	qsymm[11][3] = 0.0;			/*[1-10]180*/
	qsymm[12][0] = 0.0;				qsymm[12][1] = _sqrt2;			qsymm[12][2] = 0.0;				qsymm[12][3] = _sqrt2;		/*[101]180*/
	qsymm[13][0] = 0.0;				qsymm[13][1] = -1.0 * _sqrt2;	qsymm[13][2] = 0.0;				qsymm[13][3] = _sqrt2;		/*[-101]180*/ //mtex332 0 sq2 0 -sq2
	qsymm[14][0] = 0.0;				qsymm[14][1] = 0.0;				qsymm[14][2] = _sqrt2;			qsymm[14][3] = _sqrt2;		/*[011]180*/ //mtex 0 0 -sq -sq
	qsymm[15][0] = 0.0;				qsymm[15][1] = 0.0;				qsymm[15][2] = -1.0 * _sqrt2;	qsymm[15][3] = _sqrt2;		/*[110]180*/

	//<111>120*degrees, three-fold symmetries
	qsymm[16][0] = half;			qsymm[16][1] = half;			qsymm[16][2] = half;			qsymm[16][3] = half;		/*[111]120*/
	qsymm[17][0] = -1.0 * half;		qsymm[17][1] = half;			qsymm[17][2] = half;			qsymm[17][3] = half;		/*[111]240*/
	qsymm[18][0] = half;			qsymm[18][1] = half;			qsymm[18][2] = -1.0 * half;		qsymm[18][3] = half;		/*[1-11]240*/
	qsymm[19][0] = -1.0 * half; 	qsymm[19][1] = half;			qsymm[19][2] = -1.0 * half;		qsymm[19][3] = half;		/*[1-11]240*/
	qsymm[20][0] = half;			qsymm[20][1] = -1.0 * half;		qsymm[20][2] = half;			qsymm[20][3] = half;		/*[-111]120*/
	qsymm[21][0] = -1.0 * half;		qsymm[21][1] = -1.0 * half;		qsymm[21][2] = half;			qsymm[21][3] = half;		/*[-111]240*/ //mtex h h -h -h
	qsymm[22][0] = half;			qsymm[22][1] = -1.0 * half;		qsymm[22][2] = -1.0 * half;		qsymm[22][3] = half;		/*[-1-11]120*/ //Mtex332-h h h -h
	qsymm[23][0] = -1.0 * half;		qsymm[23][1] = -1.0 * half;		qsymm[23][2] = -1.0 * half;		qsymm[23][3] = half;		/*[-1-11]240*/
//cout << "All m-3m fcc crystal symmetry quaternion operators loaded successfully." << endl;

	double nd[3] = { 0.0, 0.0, 1.0 }; //z-direction normal vector

	//project to standard triangle by calculating symmetric variants, first ndp
	double triposp[SYMMETRIES_IN_FCC][2];

	for (int s = 0; s < SYMMETRIES_IN_FCC; s++) {
		double qtestqsym[4];

		double qq[4] = { qsymm[s][0], qsymm[s][1], qsymm[s][2], qsymm[s][3] }; //might be superficial

		//QuatOnVector3D( qq, ndrotp, hp ); //no renormalization necessary because pure rotation not stretching
		multiplyQuaternions(qtest, qq, qtestqsym);

		//ctranspose for unit quaternions, ie. q0, -q1, -q2, -q3
		qtestqsym[1] *= -1.0;
		qtestqsym[2] *= -1.0;
		qtestqsym[3] *= -1.0;

		double hbar[3];
		QuatOnVector3D(qtestqsym, nd, hbar);


		//##MK::antipodal consistency with MTex3.3.2
		if (hbar[2] < 0.0) {
			hbar[0] *= -1.0;
			hbar[1] *= -1.0;
			hbar[2] *= -1.0;
		}


		//now get the azimuth angle theta and assure no throw out of range exceptions by the acos function, but before range limiting
		if (hbar[2] > (1.0 - DOUBLE_ACCURACY)) hbar[2] = 1.0;
		if (hbar[2] < (-1.0 + DOUBLE_ACCURACY)) hbar[2] = -1.0;

		double theta = acos(hbar[2]); //result is [0.0 <= theta <= pi]

		//the fact that theta goes to pi can cause tan singularities if theta --> pi/2 periodicities
		double tt = tan(theta / 2);
		if (tt > (1.0 - DOUBLE_ACCURACY)) { tt = 1.0; }

		double psi = 0.0; //capture discontinuity of the atan2 at 0,0
		if (fabs(hbar[0]) > DOUBLE_ACCURACY&& fabs(hbar[1]) > DOUBLE_ACCURACY) {
			psi = atan2(hbar[1], hbar[0]);
		}

		//psi is limited against pi
		double epps = 1e-10;
		if (psi > (_PI_ - epps)) { psi = _PI_; }
		if (psi < ((-1.0 * _PI_) + epps)) { psi = -1.0 * _PI_; }


		triposp[s][0] = tt * cos(psi);
		triposp[s][1] = tt * sin(psi);
	}

	//check which are in the standard triangle and closest to the center in the standard triangle x >= y > 0
	double minnorm = 10.0;
	double norm = 10.0;
	double xmin = FAIL_MYMATH_NUMERIC;
	double ymin = FAIL_MYMATH_NUMERIC;

	//double sstr = pow( 2, 0.5 ) - 1.0;
	//sstr = pow ( sstr + EPS_ENVIRONMENT, 2); //radius of the sst circle

	for (int s = 0; s < SYMMETRIES_IN_FCC; s++) {
		double x2y2 = SQR(triposp[s][0]) + SQR(triposp[s][1]);
		norm = pow(x2y2, 0.5);

		if (norm <= minnorm && triposp[s][0] >= 0.0 && triposp[s][1] >= 0.0 && triposp[s][0] >= triposp[s][1]) { //&& x2y2 <= sstr ) {
			xmin = triposp[s][0];
			ymin = triposp[s][1];
			minnorm = norm;
		}
	}

	xy[0] = xmin;
	xy[1] = ymin;
}


void mathMethods::bunge2ipfz(double phi1, double PHI, double phi2, unsigned char* rgb, double* pos)
{
	//project to fundamental region
	double bunge[3];
	bunge[0] = phi1;
	bunge[1] = PHI;
	bunge[2] = phi2;
	double qbunge[4];
	double position[2];

	euler2quaternion(bunge, qbunge);
	//cout << "Euler2Quat=" << qbunge[0] << ";" << qbunge[1] << ";" << qbunge[2] << ";" << qbunge[3] << endl;

	project2fundamentalregion_ipfz(qbunge, position);
	//cout << "Position=" << position[0] << ";" << position[1] << endl;

	if (position[0] == FAIL_MYMATH_NUMERIC || position[1] == FAIL_MYMATH_NUMERIC) {
		//color in black, a color otherwise not utilized any mismatch and failure in the function!
		rgb[0] = 0;
		rgb[1] = 0;
		rgb[2] = 0;

		pos[0] = position[0];
		pos[1] = position[1];
		return;
	}

	//heuristic approach to catch numeric cases too close to the IPFZ coloring triangle vertices 
	//"red"
	double xr = 0.0;							double yr = 0.0; //red
	double xg = pow(2.0, 0.5) - 1.0;			double yg = 0.0; //green
	double xb = (0.5 * pow(3.0, 0.5)) - 0.5;	double yb = xb; //blue

	double xo = position[0];					double yo = position[1];

	//heuristic catch to avoid running in on-the-vertex-location singularities 
	/*guideline implementation of the scaling approach Molodov
	scaling suggestion from K. Molodov HKL and TSL implement different white
	points and different strategies with respect to how to rescale the RGB
	mixing, source unknown so heuristic decision, no problem: because
	an IPF colorcoding for the m-3m, 1 symmetry is possible purpose is to show
	qualitative location, coding necessarily is unpractical for sharp ipfz
	fiber textures in particular if close to the vertices, then a power-law
	rescale might exaggerate contrast the stronger the closer to the vertices
	but the less distinct for orientations which are located near the white point*/

	if (fabs(xo - xr) < EPS_PROXIMITY_IPFZ && fabs(yo - yr) < EPS_PROXIMITY_IPFZ) {
		rgb[0] = 255; //assign pure RED
		rgb[1] = 0;
		rgb[2] = 0;

		pos[0] = position[0];
		pos[1] = position[1];
		return;
	}
	if (fabs(xo - xg) < EPS_PROXIMITY_IPFZ && fabs(yo - yg) < EPS_PROXIMITY_IPFZ) {
		rgb[0] = 0;
		rgb[1] = 255; //assign pure GREEN
		rgb[2] = 0;

		pos[0] = position[0];
		pos[1] = position[1];
		return;
	}
	if (fabs(xo - xb) < EPS_PROXIMITY_IPFZ && fabs(yo - yb) < EPS_PROXIMITY_IPFZ) {
		rgb[0] = 0;
		rgb[1] = 0;
		rgb[2] = 255; //assign pure BLUE

		pos[0] = position[0];
		pos[1] = position[1];
		return;
	}

	//identify RGB color in stereographic triangle
	double k1g = (xg - xo) * (yr - yb) - (yg - yo) * (xr - xb);
	//k1g can only be come 0 if: then the lines are parallel or coincident  xg-xo and yg-yo are either both zero or both terms the same because or naturally yr-yb = xr-xb = const  != 0catch k1g <= DOUBLE_ACCURACY)
	//which based on the shape of the triangle can be safely assumed

	double xgbar = ((xg * yo - yg * xo) * (xr - xb) - (xg - xo) * (xr * yb - yr * xb)) / k1g;
	double ygbar = ((xg * yo - yg * xo) * (yr - yb) - (yg - yo) * (xr * yb - yr * xb)) / k1g;
	double absggbar = fabs(pow((SQR(xg - xgbar) + SQR(yg - ygbar)), 0.5));
	double absogbar = fabs(pow((SQR(xo - xgbar) + SQR(yo - ygbar)), 0.5));

	double k1b = (xb - xo) * (yr - yg) - (yb - yo) * (xr - xg);
	double xbbar = ((xb * yo - yb * xo) * (xr - xg) - (xb - xo) * (xr * yg - yr * xg)) / k1b;
	double ybbar = ((xb * yo - yb * xo) * (yr - yg) - (yb - yo) * (xr * yg - yr * xg)) / k1b;
	double absbbbar = fabs(pow((SQR(xb - xbbar) + SQR(yb - ybbar)), 0.5));
	double absobbar = fabs(pow((SQR(xo - xbbar) + SQR(yo - ybbar)), 0.5));

	double k1r = pow((yo / xo), 2.0);
	//x0 == zero already excluded by R vertex proximity check and the coordiante system choice x>=y
	double xrbar = pow((2 + k1r), 0.5);
	xrbar = xrbar - 1;
	xrbar /= (k1r + 1.0);
	double yrbar = yo / xo * xrbar;
	double absrrbar = fabs(pow((SQR(xr - xrbar) + SQR(yr - yrbar)), 0.5));
	double absorbar = fabs(pow((SQR(xo - xrbar) + SQR(yo - yrbar)), 0.5));

	//IPF color stretch approach
	double rrggbb[3];
	rrggbb[0] = pow((absorbar / absrrbar), IPF_COLOR_STRETCH_R);
	rrggbb[1] = pow((absogbar / absggbar), IPF_COLOR_STRETCH_G);
	rrggbb[2] = pow((absobbar / absbbbar), IPF_COLOR_STRETCH_B);

	double maxx = rrggbb[0];
	if (rrggbb[1] > maxx) maxx = rrggbb[1];
	if (rrggbb[2] > maxx) maxx = rrggbb[2];

	//K. Molodov got a better agrreement with the HKL colo though by anisotropically IPF_COLOR_STRETCHING via reverse engineering
	int rr = rrggbb[0] * (1.0 / maxx) * 255;
	int gg = rrggbb[1] * (1.0 / maxx) * 255;
	int bb = rrggbb[2] * (1.0 / maxx) * 255;

	//"RRGGBB=" << rr << ";" << gg << ";" << bb << endl;

	rgb[0] = rr;
	rgb[1] = gg;
	rgb[2] = bb;

	pos[0] = position[0];
	pos[1] = position[1];
}
