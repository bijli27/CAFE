#pragma once
#ifndef _mymath_h_
#define _mymath_h_

#include <stdlib.h>
#include <math.h>"
#include "random.h"
#include <stdint.h>

#define DOUBLE_ACCURACY 		5e-16				//2^-51 as determined with MATLAB eps(2)
#define MYMATH_STANDARD_SEED	-4256
#define MINIMUM_ANGULAR_SPREAD	(0.017453292)		//1.0/180.0_PI_
#define MAXIMUM_MISORI_FCC		(1.09606677)		//FCC 62.8/180.0*_PI_
#define SYMMETRIES_IN_FCC		24

//IPFZ Coloring
#define FAIL_MYMATH_NUMERIC		(-1.0)
#define EPS_PROXIMITY_IPFZ		(0.01)
#define IPF_COLOR_STRETCH_R		(0.5)
#define IPF_COLOR_STRETCH_G		(0.5)
#define IPF_COLOR_STRETCH_B		(0.5)
#define EPS_ENVIRONMENT			(1e-7)
#define RGBRANGE				255

class randomClass;
typedef double Real;

struct Point{
	double x;
	double y;
	double z;
};

class mathMethods //: public randomClass
{
public:
	mathMethods( void ); 		//constructor
	~mathMethods( void );		//destructor

	Real areaPolygon( int *poly, int counter );
	int pnpoly(int nvert, Real *vertx, Real *verty, Real testx, Real testy);
	void bubbleSort ( Real arr [ ], int size );
	void bubbleSortIndex( Real arr [], long idx[], int size );
	void sortInt( long arr [], long size );
	void swap ( Real& x, Real& y );
	void swapInt ( long& x, long& y );
	Real convertConc(Real conc,int j);
	
	//void newOrientationFromReference( Real const *oriOri, Real deviation, Real *newOri);

	void QuatOnVector3D(double* q, double* v, double* r);

	void factorizeIn3( long n, long * factors );
	char isPrime( long n );
	void counting_sort(int *A, int size, int range);
	void K_S_Test( double * data1, unsigned long n1, double * data2, unsigned long n2, double * d, double * prob );
	void sort(int n, double *ra);
	double probks( double alam );
	void newOrientationFromReference( double *oriOri, double deviation, double *newOri );
	void quaternion2Euler( double * quat, double * euler );
	void euler2quaternion( double * euler, double * q );
	void rotateOrientation( double *oriOri, double angle, double u, double v, double w, double *newOri );
	void newOrientationFromReferenceFixedAngularCone(double * oriOri, double maxDev,double angle, double u, double v, double w,double * newOri);
	double misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 );
	void misorientationQuaternionCubic( double* p, double* q, double* quat  );
	void multiplyQuaternions( double *q, double* p, double* r );
	void randomMisorientation( double theta, double* qr  );
	double angleBetweenQuaternions( double * q, double * p );
	void randomMisorientationAxisConsidered(  double * qref, double * qr, double maxDev  );
	double misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 );
	void randomOrientation( double * result );

	//identify orientation by RGB scheme
	void devtorefEuler2RGB(double* bunge, double* ideal, double maxDev, unsigned char* rgb); //blue channel stretch from 0.0 to maxDev in radian, all other orientations white
	void project2fundamentalregion_ipfz(double* qtest, double* xy);
	void bunge2ipfz(double phi1, double PHI, double phi2, unsigned char* rgb, double* pos);
	void linearRegression(Real* x, Real* y, long size, Real* slope, Real* intercept, Real* dev);

	randomClass r;
};

#endif