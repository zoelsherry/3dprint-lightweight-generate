#ifndef MESH_H
#define MESH_H

#include "Geometry.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>


//mesh cell
template<typename T>
struct GridCell
{

	Point3D<T> p[8];    //position of each corner of the grid in world space
	T val[8];    //value of the function at this grid 

};


struct TriFaceIndex
{
	int I[3];
	TriFaceIndex() {}
	TriFaceIndex(int I0, int I1, int i2)
	{
		I[0] = I0;
		I[1] = I1;
		I[2] = i2;
	}

};


template<typename T>
struct triMeshFace
{
	Point3D<T> vec1;
	Point3D<T> vec2;
	Point3D<T> vec3;
	triMeshFace()
	{
		vec1 = vec2 = vec3 = Point3D<T>(0.0, 0.0, 0.0);
	}
	triMeshFace(Point3D<T> _vec1, Point3D<T> _vec2, Point3D<T> _vec3)
		:vec1(_vec1), vec2(_vec2), vec3(_vec3) {}

};


template<typename T>
struct STLtriMeshFace
{
	Point3D<T> normal;
	Point3D<T> vec1;
	Point3D<T> vec2;
	Point3D<T> vec3;
	STLtriMeshFace()
	{
		normal = vec1 = vec2 = vec3 = Point3D<T>(0.0, 0.0, 0.0);
	}
	STLtriMeshFace(Point3D<T> _normal, Point3D<T> _vec1, Point3D<T> _vec2, Point3D<T> _vec3)
		:normal(_normal),vec1(_vec1),vec2(_vec2),vec3(_vec3){}

};

template<typename T>
class Mesh
{
public:
	Mesh();
	~Mesh();

	int AddVertex(Point3D<T> po);
	int AddFace(TriFaceIndex tri);
	void ClearDate();
	

	std::vector<Point3D<T>> m_normal;
	std::vector<Point3D<T>> m_vertices;
	std::vector<TriFaceIndex> m_faces_index;
	std::vector<triMeshFace<T>> m_trimesh;
	std::vector<STLtriMeshFace<T>> m_stltrimesh;
private:

};




template<typename T>
Mesh<T>::Mesh()
{
}

template<typename T>
Mesh<T>::~Mesh()
{
}

template<typename T>
int Mesh<T>::AddVertex(Point3D<T> po)
{
	int index = m_vertices.size();
	m_vertices.push_back(po);
	return index;
}
template<typename T>
int Mesh<T>::AddFace(TriFaceIndex tri)
{
	int index = m_faces_index.size();
	m_faces_index.push_back(tri);
	return index;
}

template<typename T>
void Mesh<T>::ClearDate()
{
	m_vertices.clear();
	m_faces_index.clear();
}



#endif // !MESH_H
