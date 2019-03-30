#include <math.h>
#include <iomanip>
#include "point.h"

//Defnitions for POINT3D
POINT3D::POINT3D(const DP xx, const DP yy, const DP zz) 
{
	x=xx;
	y=yy;
	z=zz;
}
	
POINT3D::POINT3D(const POINT3D &point) 
{
	if(this != &point)
	{
	    x=point.x;
	    y=point.y;
	    z=point.z;
	}
}

POINT3D& POINT3D::operator=(const POINT3D temp)
{
	x=temp.x;
	y=temp.y;
	z=temp.z;
	return *this;
}

POINT3D POINT3D::operator/(const DP c)
{
	DP d=c;	

	if (d == 0.0) d=TINY;
	return POINT3D(x/d,y/d,z/d);
}	

POINT3D& POINT3D::operator/=(const DP c)
{
	*this=(*this)/c;
	return *this;
}