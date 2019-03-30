#ifndef DISTANCE_3D_H
#define DISTANCE_3D_H

#include <vector>
#include "Geometry.h"

template<typename T>
class SignedDistSolver
{

public:
	SignedDistSolver() {};
	SignedDistSolver(gridsize ngrid, std::vector<Point3D<T>> nvertices, std::vector<Point3D<T>> nnormal, int nPIXEL)
		:XPIXEL(ngrid.nx), YPIXEL(ngrid.ny), ZPIXEL(ngrid.nz), vertices(nvertices), normal(nnormal), PIXEL(nPIXEL) {};
	~SignedDistSolver();
	void initial_dist(T*** u, int*** cp, bool *** lock);
	void updatelock(bool*** lock, int i, int j, int k);
	void updatelock(bool*** lock, int i, int j, int k, T*** u);

	void distance_sweeping(T*** u0, T*** ud, int*** cp, bool*** lock);
	void update(T*** u, int*** cp, int i, int j, int k, bool*** lock);
	void update3(T*** u, int i, int j, int k, int m);
	bool update4(T*** u, int i, int j, int k, int m);

	void calc_signed_dist(T*** dis0, int*** idx, T*** phi);
	static int sign_function(T pr);

	void unify(T*** u, int N);

private:
	int XPIXEL;
	int YPIXEL;
	int ZPIXEL;
	int PIXEL;
	std::vector<Point3D<T>> vertices;
	std::vector<Point3D<T>> normal;


};

#endif // !DISTANCE_3D_H

