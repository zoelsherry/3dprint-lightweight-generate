#ifndef _PARAMETER3D_H
#define _PARAMETER3D_H

#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <ctime>
#include "Geometry.h"



/***********constants*********/
#define CONSTANT 3.141592653589793/180.0
#define sigma  0.008
#define error pow(10.,-10)

typedef struct heap
{
	int i, j, k;
}heapnode;


template<typename T>
class DataUnifier
{
public:
	DataUnifier() {};
	DataUnifier(std::vector<Point3D<T>> _vertices, std::vector<Point3D<T>> _normal, int _PIXEL);
	~DataUnifier();
	gridsize unify_data(T scale, int box);
	int simplify_data();
	std::vector<Point3D<T>> get_simple_vertices();
	std::vector<Point3D<T>> get_simple_normal();

private:
	int PIXEL;
	gridsize grid;
	std::vector<Point3D<T>> vertices;
	std::vector<Point3D<T>> normal;

};


#endif


