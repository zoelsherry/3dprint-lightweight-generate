/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef GEOMETRY_INCLUDED
#define GEOMETRY_INCLUDED

#include <stdio.h>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <unordered_map>
#include <climits>


template<class Real>
Real Random(void);

template< class Real >
struct Point3D
{
	Real coords[3];
	Point3D( void ) { coords[0] = coords[1] = coords[2] = Real(0); }
	Point3D( Real v ) { coords[0] = coords[1] = coords[2] = v; }
	template< class _Real > Point3D( _Real v0 , _Real v1 , _Real v2 ){ coords[0] = Real(v0) , coords[1] = Real(v1) , coords[2] = Real(v2); }
	template< class _Real > Point3D( const Point3D< _Real >& p ){ coords[0] = Real( p[0] ) , coords[1] = Real( p[1] ) , coords[2] = Real( p[2] ); }
	inline       Real& operator[] ( int i )       { return coords[i]; }
	inline const Real& operator[] ( int i ) const { return coords[i]; }
	inline Point3D  operator - ( void ) const { Point3D q ; q.coords[0] = -coords[0] , q.coords[1] = -coords[1] , q.coords[2] = -coords[2] ; return q; }

	template< class _Real > inline Point3D& operator += ( Point3D< _Real > p ){ coords[0] += Real(p.coords[0]) , coords[1] += Real(p.coords[1]) , coords[2] += Real(p.coords[2]) ; return *this; }
	template< class _Real > inline Point3D  operator +  ( Point3D< _Real > p ) const { Point3D q ; q.coords[0] = coords[0] + Real(p.coords[0]) , q.coords[1] = coords[1] + Real(p.coords[1]) , q.coords[2] = coords[2] + Real(p.coords[2]) ; return q; }
	template< class _Real > inline Point3D& operator *= ( _Real r ) { coords[0] *= Real(r) , coords[1] *= Real(r) , coords[2] *= Real(r) ; return *this; }
	template< class _Real > inline Point3D  operator *  ( _Real r ) const { Point3D q ; q.coords[0] = coords[0] * Real(r) , q.coords[1] = coords[1] * Real(r) , q.coords[2] = coords[2] * Real(r) ; return q; }
	// TODO: supply the left multiple!

	template< class _Real > inline Point3D& operator -=	 ( Point3D< _Real > p ){ return ( (*this)+=(-p) ); }
	template< class _Real > inline Point3D  operator -  ( Point3D< _Real > p ) const { return (*this)+(-p); }
	template< class _Real > inline Point3D& operator /= ( _Real r ){ return ( (*this)*=Real(1./r) ); }
	template< class _Real > inline Point3D  operator /  ( _Real r ) const 
	{
		_Real d = r;
		if (d == 0.0) d = std::numeric_limits<_Real>::epsilon();
		return (*this) * ( Real(1.)/d ); 
	}

	template< class _Real > inline bool operator == (Point3D< _Real > p) 
	{ 
		if (fabs(coords[0] - p.coords[0]) < std::numeric_limits<_Real>::epsilon() &&
			fabs(coords[1] - p.coords[1]) < std::numeric_limits<_Real>::epsilon() && 
			fabs(coords[2] - p.coords[2]) < std::numeric_limits<_Real>::epsilon() )
			return true; 
		else return false; 
	}

	template< class _Real > inline bool operator == (Point3D< _Real > p) const
	{
		if (fabs(coords[0] - p.coords[0]) < std::numeric_limits<_Real>::epsilon() &&
			fabs(coords[1] - p.coords[1]) < std::numeric_limits<_Real>::epsilon() &&
			fabs(coords[2] - p.coords[2]) < std::numeric_limits<_Real>::epsilon())
			return true;
		else return false;
	}

	template< class _Real > inline Point3D& operator = (const Point3D< _Real > p)
	{
		coords[0] = p.coords[0];
		coords[1] = p.coords[1];
		coords[2] = p.coords[2];
		return *this;
	}

	template< class _Real > inline bool operator < (Point3D< _Real > p)
	{
		if (coords[0] < p.coords[0] && fabs(coords[0] - p.coords[0]) > std::numeric_limits<_Real>::epsilon()) return true;
		else if (fabs(coords[0] - p.coords[0]) < std::numeric_limits<_Real>::epsilon() && coords[1] < p.coords[1] && fabs(coords[1] - p.coords[1]) > std::numeric_limits<_Real>::epsilon()) return true;
		else if (fabs(coords[0] - p.coords[0]) < std::numeric_limits<_Real>::epsilon() &&
			fabs(coords[1] - p.coords[1]) < std::numeric_limits<_Real>::epsilon() && coords[2] < p.coords[2] && fabs(coords[2] - p.coords[2]) > std::numeric_limits<_Real>::epsilon()) return true;
		else return false;
	}

	template< class _Real > inline bool operator < (Point3D< _Real > p) const
	{
		if (coords[0] < p.coords[0] && fabs(coords[0] - p.coords[0]) > std::numeric_limits<_Real>::epsilon()) return true;
		else if (fabs(coords[0] - p.coords[0]) < std::numeric_limits<_Real>::epsilon() && coords[1] < p.coords[1] && fabs(coords[1] - p.coords[1]) > std::numeric_limits<_Real>::epsilon()) return true;
		else if (fabs(coords[0] - p.coords[0]) < std::numeric_limits<_Real>::epsilon() &&
			fabs(coords[1] - p.coords[1]) < std::numeric_limits<_Real>::epsilon() && coords[2] < p.coords[2] && fabs(coords[2] - p.coords[2]) > std::numeric_limits<_Real>::epsilon()) return true;
		else return false;
	}

	template< class _Real > inline bool operator > (Point3D< _Real > p)
	{
		if (coords[0] > p.coords[0] && fabs(coords[0] - p.coords[0]) > std::numeric_limits<_Real>::epsilon()) return true;
		else if (fabs(coords[0] - p.coords[0]) < std::numeric_limits<_Real>::epsilon() && coords[1] > p.coords[1] && fabs(coords[1] - p.coords[1]) > std::numeric_limits<_Real>::epsilon()) return true;
		else if (fabs(coords[0] - p.coords[0]) < std::numeric_limits<_Real>::epsilon() &&
			fabs(coords[1] - p.coords[1]) < std::numeric_limits<_Real>::epsilon() && coords[2] > p.coords[2] && fabs(coords[2] - p.coords[2]) > std::numeric_limits<_Real>::epsilon()) return true;
		else return false;
	}


	static Real Dot( const Point3D< Real >& p1 , const Point3D< Real >& p2 ){ return p1.coords[0]*p2.coords[0] + p1.coords[1]*p2.coords[1] + p1.coords[2]*p2.coords[2]; }
	template< class Real1 , class Real2 >
	static Real Dot( const Point3D< Real1 >& p1 , const Point3D< Real2 >& p2 ){ return Real( p1.coords[0]*p2.coords[0] + p1.coords[1]*p2.coords[1] + p1.coords[2]*p2.coords[2] ); }


	static Point3D Normalize(Point3D< Real >& p)
	{
		Real r = sqrt(pow(p.coords[0], 2) + pow(p.coords[1], 2) + pow(p.coords[2], 2));

		return p /= r;
	}


};

struct gridsize
{
	int nx;
	int ny;
	int nz;
};



template< class Real >
struct OrientedPoint3D
{
	Point3D< Real > p , n;
	OrientedPoint3D( Point3D< Real > pp=Point3D< Real >() , Point3D< Real > nn=Point3D< Real >() ) : p(pp) , n(nn) { ; }
	template< class _Real > OrientedPoint3D( const OrientedPoint3D< _Real >& p ) : OrientedPoint3D( Point3D< Real >( p.p ) , Point3D< Real >( p.n ) ){ ; }

	template< class _Real > inline OrientedPoint3D& operator += ( OrientedPoint3D< _Real > _p ){ p += _p.p , n += _p.n ; return *this; }
	template< class _Real > inline OrientedPoint3D  operator +  ( OrientedPoint3D< _Real > _p ) const { return OrientedPoint3D< Real >( p+_p.p , n+_p.n ); }
	template< class _Real > inline OrientedPoint3D& operator *= ( _Real r ) { p *= r , n *= r ; return *this; }
	template< class _Real > inline OrientedPoint3D  operator *  ( _Real r ) const { return OrientedPoint3D< Real >( p*r , n*r ); }

	template< class _Real > inline OrientedPoint3D& operator -= ( OrientedPoint3D< _Real > p ){ return ( (*this)+=(-p) ); }
	template< class _Real > inline OrientedPoint3D  operator -  ( OrientedPoint3D< _Real > p ) const { return (*this)+(-p); }
	template< class _Real > inline OrientedPoint3D& operator /= ( _Real r ){ return ( (*this)*=Real(1./r) ); }
	template< class _Real > inline OrientedPoint3D  operator /  ( _Real r ) const { return (*this) * ( Real(1.)/r ); }

	template< class _Real > inline bool operator == ( OrientedPoint3D< _Real > _p ) { if (p[0] == _p.p[0] && p[1] == _p.p[1] && p[2] == _p.p[2]) return true; else return false; }
	template< class _Real > inline bool operator < ( OrientedPoint3D< _Real > _p )
	{
		if (p[0] < _p.p[0]) return true;
		else if (p[0] == _p.p[0] && p[1] < _p.p[1]) return true;
		else if (p[0] == _p.p[0] && p[1] == _p.p[1] && p[2] < _p.p[2]) return true;
		else return false;
	}
	template< class _Real > inline bool operator > ( OrientedPoint3D< _Real > _p )
	{
		if (p[0] > _p.p[0]) return true;
		else if (p[0] == _p.p[0] && p[1] > _p.p[1]) return true;
		else if (p[0] == _p.p[0] && p[1] == _p.p[1] && p[2] > _p.p[2]) return true;
		else return false;
	}

};


template<class Real>
Point3D<Real> RandomBallPoint(void);

template<class Real>
Point3D<Real> RandomSpherePoint(void);

template<class Real>
double Length(const Point3D<Real>& p);

template<class Real>
double SquareLength(const Point3D<Real>& p);

template<class Real>
double Distance(const Point3D<Real>& p1,const Point3D<Real>& p2);

template<class Real>
double SquareDistance(const Point3D<Real>& p1,const Point3D<Real>& p2);

template <class Real>
void CrossProduct(const Point3D<Real>& p1,const Point3D<Real>& p2,Point3D<Real>& p);




/////////////////////////////////////////////////////////////////////


template<class Real>
Real Random(void) { return Real(rand()) / RAND_MAX; }

template<class Real>
Point3D<Real> RandomBallPoint(void) {
	Point3D<Real> p;
	while (1) {
		p.coords[0] = Real(1.0 - 2.0*Random<Real>());
		p.coords[1] = Real(1.0 - 2.0*Random<Real>());
		p.coords[2] = Real(1.0 - 2.0*Random<Real>());
		double l = SquareLength(p);
		if (l <= 1) { return p; }
	}
}
template<class Real>
Point3D<Real> RandomSpherePoint(void) {
	Point3D<Real> p = RandomBallPoint<Real>();
	Real l = Real(Length(p));
	p.coords[0] /= l;
	p.coords[1] /= l;
	p.coords[2] /= l;
	return p;
}

template<class Real>
double SquareLength(const Point3D<Real>& p) { return p.coords[0] * p.coords[0] + p.coords[1] * p.coords[1] + p.coords[2] * p.coords[2]; }

template<class Real>
double Length(const Point3D<Real>& p) { return sqrt(SquareLength(p)); }

template<class Real>
double SquareDistance(const Point3D<Real>& p1, const Point3D<Real>& p2) {
	return (p1.coords[0] - p2.coords[0])*(p1.coords[0] - p2.coords[0]) + (p1.coords[1] - p2.coords[1])*(p1.coords[1] - p2.coords[1]) + (p1.coords[2] - p2.coords[2])*(p1.coords[2] - p2.coords[2]);
}

template<class Real>
double Distance(const Point3D<Real>& p1, const Point3D<Real>& p2) { return sqrt(SquareDistance(p1, p2)); }

template <class Real>
void CrossProduct(const Point3D<Real>& p1, const Point3D<Real>& p2, Point3D<Real>& p) {
	p.coords[0] = p1.coords[1] * p2.coords[2] - p1.coords[2] * p2.coords[1];
	p.coords[1] = -p1.coords[0] * p2.coords[2] + p1.coords[2] * p2.coords[0];
	p.coords[2] = p1.coords[0] * p2.coords[1] - p1.coords[1] * p2.coords[0];
}


#endif // GEOMETRY_INCLUDED
