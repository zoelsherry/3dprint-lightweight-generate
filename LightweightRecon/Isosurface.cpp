#include "Isosurface.h"
#include <iostream>
using namespace std;

template<typename T>
Isosurface<T>::~Isosurface()
{
	
}

template<typename T>
Mesh<T> Isosurface<T>::GenerateIsosurface(T*** nphi)
{
	for (int i = 0; i < m_nCellsX - 1; i++)
	{
		for (int j = 0; j < m_nCellsY - 1; j++)
		{
			for (int k = 0; k < m_nCellsZ - 1; k++)
			{
			
				GridCell<T> cell;
				cell.val[7] = nphi[i][j][k];
				cell.val[6] = nphi[i + 1][j][k];
				cell.val[5] = nphi[i + 1][j][k + 1];
				cell.val[4] = nphi[i][j][k + 1];
				cell.val[3] = nphi[i][j + 1][k];
				cell.val[2] = nphi[i + 1][j + 1][k];
				cell.val[1] = nphi[i + 1][j + 1][k + 1];
				cell.val[0] = nphi[i][j + 1][k + 1];
				
				/*cell.p[7] = Point3D<double>(x[i], y[j], z[k]);
				cell.p[6] = Point3D<double>(x[i + 1], y[j], z[k]);
				cell.p[5] = Point3D<double>(x[i + 1], y[j], z[k + 1]);
				cell.p[4] = Point3D<double>(x[i], y[j], z[k + 1]);
				cell.p[3] = Point3D<double>(x[i], y[j + 1], z[k]);
				cell.p[2] = Point3D<double>(x[i + 1], y[j + 1], z[k]);
				cell.p[1] = Point3D<double>(x[i + 1], y[j + 1], z[k + 1]);
				cell.p[0] = Point3D<double>(x[i], y[j + 1], z[k + 1]);*/
									
				cell.p[7] = Point3D<T>(i+1    , j+1    , k+1    );
				cell.p[6] = Point3D<T>(i+1 + 1, j+1    , k+1    );
				cell.p[5] = Point3D<T>(i+1 + 1, j+1    , k+1 + 1);
				cell.p[4] = Point3D<T>(i+1    , j+1    , k+1 + 1);
				cell.p[3] = Point3D<T>(i+1    , j+1 + 1, k+1    );
				cell.p[2] = Point3D<T>(i+1 + 1, j+1 + 1, k+1    );
				cell.p[1] = Point3D<T>(i+1 + 1, j+1 + 1, k+1 + 1);
				cell.p[0] = Point3D<T>(i+1    , j+1 + 1, k+1 + 1);

				ExtractTriangles(cell);
			}
		}
	}

	return builder.getMesh();
}

template<typename T>
Mesh<T> Isosurface<T>::GenerateSTLIsosurface(T *** nphi)
{
	for (int i = 0; i < m_nCellsX - 1; i++)
	{
		for (int j = 0; j < m_nCellsY - 1; j++)
		{
			for (int k = 0; k < m_nCellsZ - 1; k++)
			{

				GridCell<T> cell;
				cell.val[7] = nphi[i][j][k];
				cell.val[6] = nphi[i + 1][j][k];
				cell.val[5] = nphi[i + 1][j][k + 1];
				cell.val[4] = nphi[i][j][k + 1];
				cell.val[3] = nphi[i][j + 1][k];
				cell.val[2] = nphi[i + 1][j + 1][k];
				cell.val[1] = nphi[i + 1][j + 1][k + 1];
				cell.val[0] = nphi[i][j + 1][k + 1];

				/*cell.p[7] = Point3D<double>(x[i], y[j], z[k]);
				cell.p[6] = Point3D<double>(x[i + 1], y[j], z[k]);
				cell.p[5] = Point3D<double>(x[i + 1], y[j], z[k + 1]);
				cell.p[4] = Point3D<double>(x[i], y[j], z[k + 1]);
				cell.p[3] = Point3D<double>(x[i], y[j + 1], z[k]);
				cell.p[2] = Point3D<double>(x[i + 1], y[j + 1], z[k]);
				cell.p[1] = Point3D<double>(x[i + 1], y[j + 1], z[k + 1]);
				cell.p[0] = Point3D<double>(x[i], y[j + 1], z[k + 1]);*/

				cell.p[7] = Point3D<T>(i + 1, j + 1, k + 1);
				cell.p[6] = Point3D<T>(i + 1 + 1, j + 1, k + 1);
				cell.p[5] = Point3D<T>(i + 1 + 1, j + 1, k + 1 + 1);
				cell.p[4] = Point3D<T>(i + 1, j + 1, k + 1 + 1);
				cell.p[3] = Point3D<T>(i + 1, j + 1 + 1, k + 1);
				cell.p[2] = Point3D<T>(i + 1 + 1, j + 1 + 1, k + 1);
				cell.p[1] = Point3D<T>(i + 1 + 1, j + 1 + 1, k + 1 + 1);
				cell.p[0] = Point3D<T>(i + 1, j + 1 + 1, k + 1 + 1);

				ExtractSTLTriangles(cell);
			}
		}
	}

	return builder.getMesh();
}


template<typename T>
void Isosurface<T>::ExtractTriangles(GridCell<T> cell)
{

	//Determine the index into the edge table which
	//tells us which vertices are inside of the surface
	int CubeIndex;
	CubeIndex = 0;
	//CubeIndex = MarchingCubes::GetIndex(cell.val, 0.0f);
	if (cell.val[0] > 0.0f) CubeIndex |= 1;
	if (cell.val[1] > 0.0f) CubeIndex |= 2;
	if (cell.val[2] > 0.0f) CubeIndex |= 4;
	if (cell.val[3] > 0.0f) CubeIndex |= 8;
	if (cell.val[4] > 0.0f) CubeIndex |= 16;
	if (cell.val[5] > 0.0f) CubeIndex |= 32;
	if (cell.val[6] > 0.0f) CubeIndex |= 64;
	if (cell.val[7] > 0.0f) CubeIndex |= 128;

	//Cube is entirely in/out of the surface
	if (MarchingCubes::edgeMask[CubeIndex] == 0)
		return;

	//Find the vertices where the surface intersects the cube
	Point3D<double> VertexList[12];
	
	if (MarchingCubes::edgeMask[CubeIndex] & 1)
		VertexList[0] =
		VertexInterp(cell.p[0], cell.p[1], cell.val[0], cell.val[1]);
	if (MarchingCubes::edgeMask[CubeIndex] & 2)
		VertexList[1] =
		VertexInterp(cell.p[1], cell.p[2], cell.val[1], cell.val[2]);
	if (MarchingCubes::edgeMask[CubeIndex] & 4)
		VertexList[2] =
		VertexInterp(cell.p[2], cell.p[3], cell.val[2], cell.val[3]);
	if (MarchingCubes::edgeMask[CubeIndex] & 8)
		VertexList[3] =
		VertexInterp(cell.p[3], cell.p[0], cell.val[3], cell.val[0]);
	if (MarchingCubes::edgeMask[CubeIndex] & 16)
		VertexList[4] =
		VertexInterp(cell.p[4], cell.p[5], cell.val[4], cell.val[5]);
	if (MarchingCubes::edgeMask[CubeIndex] & 32)
		VertexList[5] =
		VertexInterp(cell.p[5], cell.p[6], cell.val[5], cell.val[6]);
	if (MarchingCubes::edgeMask[CubeIndex] & 64)
		VertexList[6] =
		VertexInterp(cell.p[6], cell.p[7], cell.val[6], cell.val[7]);
	if (MarchingCubes::edgeMask[CubeIndex] & 128)
		VertexList[7] =
		VertexInterp(cell.p[7], cell.p[4], cell.val[7], cell.val[4]);
	if (MarchingCubes::edgeMask[CubeIndex] & 256)
		VertexList[8] =
		VertexInterp(cell.p[0], cell.p[4], cell.val[0], cell.val[4]);
	if (MarchingCubes::edgeMask[CubeIndex] & 512)
		VertexList[9] =
		VertexInterp(cell.p[1], cell.p[5], cell.val[1], cell.val[5]);
	if (MarchingCubes::edgeMask[CubeIndex] & 1024)
		VertexList[10] =
		VertexInterp(cell.p[2], cell.p[6], cell.val[2], cell.val[6]);
	if (MarchingCubes::edgeMask[CubeIndex] & 2048)
		VertexList[11] =
		VertexInterp(cell.p[3], cell.p[7], cell.val[3], cell.val[7]);
	
	Point3D<T> p0;
	Point3D<T> p1;
	Point3D<T> p2;
	for (int i = 0; MarchingCubes::triangles[CubeIndex][i] != -1; i += 3)
	{
		p0 = VertexList[MarchingCubes::triangles[CubeIndex][i + 0]];
		p1 = VertexList[MarchingCubes::triangles[CubeIndex][i + 1]];
		p2 = VertexList[MarchingCubes::triangles[CubeIndex][i + 2]];

		builder.AddTriangle(p0, p1, p2);

	}
	//

	/*cout << p0[0] << "," << p0[1] << "," << p0[2] << endl;
	cout << p1[0] << "," << p1[1] << "," << p1[2] << endl;
	cout << p2[0] << "," << p2[1] << "," << p2[2] << endl;*/



}

template<typename T>
void Isosurface<T>::ExtractSTLTriangles(GridCell<T> cell)
{
	//Determine the index into the edge table which
	//tells us which vertices are inside of the surface
	int CubeIndex;
	CubeIndex = 0;
	//CubeIndex = MarchingCubes::GetIndex(cell.val, 0.0f);
	if (cell.val[0] > 0.0f) CubeIndex |= 1;
	if (cell.val[1] > 0.0f) CubeIndex |= 2;
	if (cell.val[2] > 0.0f) CubeIndex |= 4;
	if (cell.val[3] > 0.0f) CubeIndex |= 8;
	if (cell.val[4] > 0.0f) CubeIndex |= 16;
	if (cell.val[5] > 0.0f) CubeIndex |= 32;
	if (cell.val[6] > 0.0f) CubeIndex |= 64;
	if (cell.val[7] > 0.0f) CubeIndex |= 128;

	//Cube is entirely in/out of the surface
	if (MarchingCubes::edgeMask[CubeIndex] == 0)
		return;

	//Find the vertices where the surface intersects the cube
	Point3D<double> VertexList[12];

	if (MarchingCubes::edgeMask[CubeIndex] & 1)
		VertexList[0] =
		VertexInterp(cell.p[0], cell.p[1], cell.val[0], cell.val[1]);
	if (MarchingCubes::edgeMask[CubeIndex] & 2)
		VertexList[1] =
		VertexInterp(cell.p[1], cell.p[2], cell.val[1], cell.val[2]);
	if (MarchingCubes::edgeMask[CubeIndex] & 4)
		VertexList[2] =
		VertexInterp(cell.p[2], cell.p[3], cell.val[2], cell.val[3]);
	if (MarchingCubes::edgeMask[CubeIndex] & 8)
		VertexList[3] =
		VertexInterp(cell.p[3], cell.p[0], cell.val[3], cell.val[0]);
	if (MarchingCubes::edgeMask[CubeIndex] & 16)
		VertexList[4] =
		VertexInterp(cell.p[4], cell.p[5], cell.val[4], cell.val[5]);
	if (MarchingCubes::edgeMask[CubeIndex] & 32)
		VertexList[5] =
		VertexInterp(cell.p[5], cell.p[6], cell.val[5], cell.val[6]);
	if (MarchingCubes::edgeMask[CubeIndex] & 64)
		VertexList[6] =
		VertexInterp(cell.p[6], cell.p[7], cell.val[6], cell.val[7]);
	if (MarchingCubes::edgeMask[CubeIndex] & 128)
		VertexList[7] =
		VertexInterp(cell.p[7], cell.p[4], cell.val[7], cell.val[4]);
	if (MarchingCubes::edgeMask[CubeIndex] & 256)
		VertexList[8] =
		VertexInterp(cell.p[0], cell.p[4], cell.val[0], cell.val[4]);
	if (MarchingCubes::edgeMask[CubeIndex] & 512)
		VertexList[9] =
		VertexInterp(cell.p[1], cell.p[5], cell.val[1], cell.val[5]);
	if (MarchingCubes::edgeMask[CubeIndex] & 1024)
		VertexList[10] =
		VertexInterp(cell.p[2], cell.p[6], cell.val[2], cell.val[6]);
	if (MarchingCubes::edgeMask[CubeIndex] & 2048)
		VertexList[11] =
		VertexInterp(cell.p[3], cell.p[7], cell.val[3], cell.val[7]);

	Point3D<T> p0;
	Point3D<T> p1;
	Point3D<T> p2;
	for (int i = 0; MarchingCubes::triangles[CubeIndex][i] != -1; i += 3)
	{
		p0 = VertexList[MarchingCubes::triangles[CubeIndex][i + 0]];
		p1 = VertexList[MarchingCubes::triangles[CubeIndex][i + 1]];
		p2 = VertexList[MarchingCubes::triangles[CubeIndex][i + 2]];

		builder.m_mesh.m_trimesh.push_back(triMeshFace<T>(p0, p1, p2));

	}
}

/*
Linearly interpolate the position where an isosurface cuts
an edge between two vertices, each with their own scalar value
*/
template<typename T>
Point3D<T> Isosurface<T>::VertexInterp(const Point3D<T> &p1, const Point3D<T> &p2, T valp1, T valp2)
{
	return (p1 + (p2 - p1) * (-valp1 / (valp2 - valp1)));
}

template Isosurface<float>::Isosurface();
template Isosurface<float>::~Isosurface();
template Isosurface<float>::Isosurface(unsigned int, unsigned int, unsigned int);
template Mesh<float> Isosurface<float>::GenerateIsosurface(float ***);
template Mesh<float> Isosurface<float>::GenerateSTLIsosurface(float ***);



template Isosurface<double>::Isosurface();
template Isosurface<double>::~Isosurface();
template Isosurface<double>::Isosurface(unsigned int, unsigned int, unsigned int);
template Mesh<double> Isosurface<double>::GenerateIsosurface(double ***);
template Mesh<double> Isosurface<double>::GenerateSTLIsosurface(double ***);


