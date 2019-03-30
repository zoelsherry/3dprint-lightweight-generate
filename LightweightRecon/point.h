#ifndef _POINT_H_
#define _POINT_H_

#include <limits>
#include <iostream>

using namespace std;

typedef double DP;	// Coincide with the definiton in "vecmat.h"
const DP TINY=numeric_limits<DP>::epsilon();

class POINT3D {
    public:
	DP x,y,z;
    	POINT3D(const DP xx=0., const DP yy=0., const DP zz=0.);
	POINT3D(const POINT3D &point);
	
	POINT3D & operator=(const POINT3D q);
	POINT3D operator/(const DP c);
	POINT3D operator*(const DP c);
	POINT3D & operator/=(const DP c);

	DP norm()	{	return sqrt(pow(x,2)+pow(y,2)+pow(z,2)); }	
};

typedef POINT3D VECTOR;
typedef POINT3D POINT;

#endif //for _POINT_H_


