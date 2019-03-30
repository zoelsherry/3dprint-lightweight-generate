#ifndef ISOSURFACE_H
#define ISOSURFACE_H

#include "Mesh.h"
#include "MeshBuilder.h"
#include "MarchingCubes.h"
#include <cmath>



template<typename T>
class Isosurface
{
public:
	Isosurface(){}
	Isosurface(unsigned int nCellsX, unsigned int nCellsY, unsigned int nCellsZ):m_nCellsX(nCellsX),m_nCellsY(nCellsY),m_nCellsZ(nCellsZ){}
	~Isosurface();

	Mesh<T> GenerateIsosurface(T*** nphi);
	Mesh<T> GenerateSTLIsosurface(T*** nphi);
	void ExtractTriangles(GridCell<T> cell);
	void ExtractSTLTriangles(GridCell<T> cell);

	Point3D<T> VertexInterp(const Point3D<T> &p1, const Point3D<T> &p2, T valp1, T valp2);
protected:
	unsigned int m_nCellsX;
	unsigned int m_nCellsY;
	unsigned int m_nCellsZ;
	MeshBuilder<T> builder;


	
};

#endif // !ISOSURFACE_H
