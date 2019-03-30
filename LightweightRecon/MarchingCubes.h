#ifndef MARCHINGCUBES_H
#define MARCHINGCUBES_H

#include <vector>
#include "Geometry.h"


#define NEW_ORDERING 1



class Cube {
public:
	const static unsigned int CORNERS = 8, EDGES = 12, FACES = 6;

	static int  CornerIndex(int x, int y, int z);
	static void FactorCornerIndex(int idx, int& x, int& y, int& z);
	static int  EdgeIndex(int orientation, int i, int j);
	static void FactorEdgeIndex(int idx, int& orientation, int& i, int &j);
	static int  FaceIndex(int dir, int offSet);
	static int  FaceIndex(int x, int y, int z);
	static void FactorFaceIndex(int idx, int& x, int &y, int& z);
	static void FactorFaceIndex(int idx, int& dir, int& offSet);

	static int  AntipodalCornerIndex(int idx);
	static int  FaceReflectCornerIndex(int idx, int faceIndex);
	static int  FaceReflectEdgeIndex(int idx, int faceIndex);
	static int	FaceReflectFaceIndex(int idx, int faceIndex);
	static int	EdgeReflectCornerIndex(int idx, int edgeIndex);
	static int	EdgeReflectEdgeIndex(int edgeIndex);

	static int  FaceAdjacentToEdges(int eIndex1, int eIndex2);
	static void FacesAdjacentToEdge(int eIndex, int& f1Index, int& f2Index);

	static void EdgeCorners(int idx, int& c1, int &c2);
	static void FaceCorners(int idx, int& c1, int &c2, int& c3, int& c4);

	static bool IsEdgeCorner(int cIndex, int e);
	static bool IsFaceCorner(int cIndex, int f);
};




class MarchingCubes
{
	static void SetVertex(int e, const double values[Cube::CORNERS], double iso);
	static unsigned char GetFaceIndex(const double values[Cube::CORNERS], double iso, int faceIndex);

	static void SetVertex(int e, const float values[Cube::CORNERS], float iso);
	static unsigned char GetFaceIndex(const float values[Cube::CORNERS], float iso, int faceIndex);

public:
	static unsigned char GetFaceIndex(unsigned char mcIndex, int faceIndex);
	static double Interpolate(double v1, double v2);
	static float Interpolate(float v1, float v2);
	const static unsigned int MAX_TRIANGLES = 5;
	static const int edgeMask[1 << Cube::CORNERS];
	static const int triangles[1 << Cube::CORNERS][3 * MAX_TRIANGLES + 1];
	static const int cornerMap[Cube::CORNERS];
	static double vertexList[Cube::EDGES][3];

	static int AddTriangleIndices(int mcIndex, int* triangles);

	static unsigned char GetIndex(const double values[Cube::CORNERS], double iso);
	static int AddTriangleIndices(const double v[Cube::CORNERS], double isoValue, int* triangles);

	static unsigned char GetIndex(const float values[Cube::CORNERS], float iso);
	static int AddTriangleIndices(const float v[Cube::CORNERS], float isoValue, int* triangles);

};

#endif // !MARCHINGCUBES_H
