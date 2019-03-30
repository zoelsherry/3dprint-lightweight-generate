#ifndef MESHBUILDER_H
#define MESHBUILDER_H

#include <map>

#include "Mesh.h"
#include "Geometry.h"
#include <cassert>

#include <iostream>

//namespace std
//{
//	template<> struct hash<Point3D<double>>
//	{
//		unsigned int operator()(const Point3D<double> & s) const noexcept
//		{
//
//			//
//			unsigned int hash = 534381;
//			int exponent;
//			double mantissa;
//			for (int i = 0; i < 3; i++)
//				mantissa = frexpf(s[i], &exponent);
//			hash += (unsigned int)((2 * mantissa - 1) * (~0U));;
//			return hash;
//		}
//	};
//}
template<typename T>
class Isosurface;

template<typename T>
class MeshBuilder
{
	friend Isosurface<T>;
public:
	MeshBuilder();
	//MeshBuilder(const std::string& filename);
	~MeshBuilder();

	void AddTriangle(Point3D<T> p0, Point3D<T> p1, Point3D<T> p2);
	Mesh<T> getMesh();
	void Clear();

	/*void load_stl(const std::string& filename);
	void ReadASCII(const std::string& filename);
	void ReadBinary(const std::string& filename);
	void export_stl(const std::string& filename, std::vector<STLtriMeshFace<float>> stltrimesh);
	void export_stl(const std::string& filename, std::vector<triMeshFace<float>> trimesh);
	void export_stl(const std::string& filename, std::vector<Point3D<float>> _vertices, std::vector<TriFaceIndex> _faces_index);
	Point3D<T> parse_point(std::ifstream& s);
	float parse_float(std::ifstream& s);*/

private:
	Mesh<T> m_mesh;
	//std::unordered_map<Point3D<double>, int> m_hashmap;
	std::map<Point3D<double>, int> m_hashmap;

	//
	int vertex_count = 0;
	//
};


template<typename T>
MeshBuilder<T>::MeshBuilder()
{
}


template<typename T>
MeshBuilder<T>::~MeshBuilder()
{

}

template<typename T>
void MeshBuilder<T>::AddTriangle(Point3D<T> p0, Point3D<T> p1, Point3D<T> p2)
{
	int p0idx;
	int p1idx;
	int p2idx;
	int index;
	//std::unordered_map<Point3D<double>, int>::iterator got = m_hashmap.find(p0);
	std::map<Point3D<double>, int>::iterator got = m_hashmap.find(p0);
	if (got == m_hashmap.end())
	{
		p0idx = m_mesh.AddVertex(p0);
		m_hashmap.emplace(p0, p0idx);
	}
	else
	{
		index = m_hashmap.at(p0);
		assert(index >= 0);
		p0idx = index;
	}

	got = m_hashmap.find(p1);
	if (got == m_hashmap.end())
	{
		p1idx = m_mesh.AddVertex(p1);
		m_hashmap.emplace(p1, p1idx);
	}
	else
	{
		index = m_hashmap.at(p1);
		assert(index >= 0);
		p1idx = index;
	}

	got = m_hashmap.find(p2);
	if (got == m_hashmap.end())
	{
		p2idx = m_mesh.AddVertex(p2);
		m_hashmap.emplace(p2, p2idx);
	}
	else
	{
		index = m_hashmap.at(p2);
		assert(index >= 0);
		p2idx = index;
	}

	//cout << m_hashmap.size() << endl;
	TriFaceIndex tri(p0idx, p1idx, p2idx);
	m_mesh.AddFace(tri);

}

template<typename T>
Mesh<T> MeshBuilder<T>::getMesh()
{
	return m_mesh;
}

template<typename T>
void MeshBuilder<T>::Clear()
{
	m_hashmap.clear();
}

// read stl file data into algorithm
//Binary Stl format
//UINT8[80] ?Header
//UINT32 ?Number of triangles
//foreach triangle
//REAL32[3] ?Normal vector
//REAL32[3] ?Vertex 1
//REAL32[3] ?Vertex 2
//REAL32[3] ?Vertex 3
//UINT16 ?Attribute byte count
//end
//template<typename T>
//void MeshBuilder<T>::load_stl(const std::string & filename)
//{
//	std::ifstream fin(filename, std::ios::in);
//	if (!fin.good())
//	{
//		std::cout << "Error: cannot open the STL file " << filename << std::endl;
//		return;
//	}
//
//	std::string headStr;
//
//	getline(fin, headStr, ' ');
//	fin.close();
//
//	if (headStr.empty())
//	{
//		std::cout << "Error: the STL file is wrong " << filename << std::endl;
//		return;
//	}
//
//	if (headStr[0] == 's')
//	{
//		ReadASCII(filename);
//	}
//	else
//	{
//		ReadBinary(filename);
//	}
//
//}
//
//template<typename T>
//void MeshBuilder<T>::ReadASCII(const std::string & filename)
//{
//}
//
//template<typename T>
//void MeshBuilder<T>::ReadBinary(const std::string & filename)
//{
//	// read binary file
//	// file open for reading
//	// operations are performed in binary mode rather than text
//	std::ifstream fin(filename, std::ios::in | std::ios::binary);
//	if (!fin.good())
//	{
//		std::cout << "Error: cannot open the STL file " << filename << std::endl;
//		return;
//	}
//	// header
//	char header_info[80] = "";
//	char n_triangles[4];
//	fin.read(header_info, 80);
//	fin.read(n_triangles, 4);
//	unsigned int* r = (unsigned int*)n_triangles;
//	unsigned int num_triangles = *r;
//	// number of triangles
//	std::cout << "Num Triangles: " << num_triangles << "\n";
//	for (unsigned int i = 0; i < num_triangles; i++)
//	{
//		// normal vector of triangle
//		auto normal = parse_point(fin);
//		// vertice of triangle
//		auto v1 = parse_point(fin);
//		auto v2 = parse_point(fin);
//		auto v3 = parse_point(fin);
//		//push into triangle list
//		m_mesh.m_stltrimesh.push_back(STLtriMeshFace<T>(normal, v1, v2, v3));
//		AddTriangle(v1, v2, v3);
//		//attributes
//		char dummy[2];
//
//		fin.read(dummy, 2);
//	}
//	std::cout << "Successfully import stl file" << "\n";
//}
//
//template<typename T>
//void MeshBuilder<T>::export_stl(const std::string & filename, std::vector<STLtriMeshFace<float>> stltrimesh)
//{
//	//binary file
//	std::string header_info = "solid " + filename + "-output";
//	char head[80];
//	std::strncpy(head, header_info.c_str(), sizeof(head));
//	//attributes
//	char attribute[2] = "0";
//	//number of triangles
//	int nTrilong = stltrimesh.size();
//	std::ofstream ofs(filename, std::ofstream::binary | std::ofstream::trunc);
//	if (!ofs.good())
//	{
//		std::cout << "cannot open output file" << filename << std::endl;
//	}
//	//write header into stl file
//	ofs.write(head, sizeof(head));
//	//write number of triangles into stl file
//	ofs.write((char*)&nTrilong, sizeof(int));
//	for (unsigned int i = 0; i < stltrimesh.size(); ++i)
//	{
//		
//		// normal vector coordinates
//
//		ofs.write((char*)&stltrimesh[i].normal[0], sizeof(float));
//		ofs.write((char*)&stltrimesh[i].normal[1], sizeof(float));
//		ofs.write((char*)&stltrimesh[i].normal[2], sizeof(float));
//
//
//		//vertex 1 coordinates
//		ofs.write((char*)&stltrimesh[i].vec1[0], sizeof(float));
//		ofs.write((char*)&stltrimesh[i].vec1[1], sizeof(float));
//		ofs.write((char*)&stltrimesh[i].vec1[2], sizeof(float));
//
//		//vertex 2 coordinates
//		ofs.write((char*)&stltrimesh[i].vec2[0], sizeof(float));
//		ofs.write((char*)&stltrimesh[i].vec2[1], sizeof(float));
//		ofs.write((char*)&stltrimesh[i].vec2[2], sizeof(float));
//
//		//vertex 3 coordinates
//		ofs.write((char*)&stltrimesh[i].vec3[0], sizeof(float));
//		ofs.write((char*)&stltrimesh[i].vec3[1], sizeof(float));
//		ofs.write((char*)&stltrimesh[i].vec3[2], sizeof(float));
//
//		//Attibutes
//		ofs.write(attribute, sizeof(short));
//	}
//
//	ofs.close();
//	std::cout << "Successfully export stl file." << std::endl;
//}
//
//template<typename T>
//void MeshBuilder<T>::export_stl(const std::string & filename, std::vector<triMeshFace<float>> trimesh)
//{
//	//binary file
//	std::string header_info = "solid " + filename + "-output";
//	char head[80];
//	std::strncpy(head, header_info.c_str(), sizeof(head));
//	//attributes
//	char attribute[2] = "0";
//	//number of triangles
//	int nTrilong = trimesh.size();
//	std::ofstream ofs(filename, std::ofstream::binary | std::ofstream::trunc);
//	if (!ofs.good())
//	{
//		std::cout << "cannot open output file" << filename << std::endl;
//	}
//	//write header into stl file
//	ofs.write(head, sizeof(head));
//	//write number of triangles into stl file
//	ofs.write((char*)&nTrilong, sizeof(int));
//	for (unsigned int i = 0; i < trimesh.size(); ++i)
//	{
//		Point3D<float> vec1 = trimesh[i].vec1;
//		Point3D<float> vec2 = trimesh[i].vec2;
//		Point3D<float> vec3 = trimesh[i].vec3;
//
//		Point3D<float> edge1 = vec3 - vec1;
//		Point3D<float> edge2 = vec2 - vec1;
//		Point3D<float> normal = Point3D<float>
//			(edge1.coords[1] * edge2.coords[2] - edge1.coords[2] * edge2.coords[1],
//				edge1.coords[2] * edge2.coords[0] - edge1.coords[0] * edge2.coords[2],
//				edge1.coords[0] * edge2.coords[1] - edge1.coords[1] * edge2.coords[0]) * (-1.0f);
//		Point3D<float>::Normalize(normal);
//		// normal vector coordinates
//
//		ofs.write((char*)&normal[0], sizeof(float));
//		ofs.write((char*)&normal[1], sizeof(float));
//		ofs.write((char*)&normal[2], sizeof(float));
//
//
//		//vertex 1 coordinates
//		ofs.write((char*)&vec1[0], sizeof(float));
//		ofs.write((char*)&vec1[1], sizeof(float));
//		ofs.write((char*)&vec1[2], sizeof(float));
//
//		//vertex 2 coordinates
//		ofs.write((char*)&vec2[0], sizeof(float));
//		ofs.write((char*)&vec2[1], sizeof(float));
//		ofs.write((char*)&vec2[2], sizeof(float));
//
//		//vertex 3 coordinates
//		ofs.write((char*)&vec3[0], sizeof(float));
//		ofs.write((char*)&vec3[1], sizeof(float));
//		ofs.write((char*)&vec3[2], sizeof(float));
//
//		//Attibutes
//		ofs.write(attribute, sizeof(short));
//	}
//
//	ofs.close();
//	std::cout << "Successfully export stl file." << std::endl;
//}
//
//template<typename T>
//void MeshBuilder<T>::export_stl(const std::string & filename, std::vector<Point3D<float>> _vertices, std::vector<TriFaceIndex> _faces_index)
//{
//	//binary file
//	std::string header_info = "solid " + filename + "-output";
//	char head[80];
//	std::strncpy(head, header_info.c_str(), sizeof(head));
//	//attributes
//	char attribute[2] = "0";
//	//number of triangles
//	int nTrilong = _faces_index.size();
//	std::ofstream ofs(filename, std::ofstream::binary | std::ofstream::trunc);
//	if (!ofs.good())
//	{
//		std::cout << "cannot open output file" << filename << std::endl;
//	}
//	//write header into stl file
//	ofs.write(head, sizeof(head));
//	//write number of triangles into stl file
//	ofs.write((char*)&nTrilong, sizeof(int));
//	for (unsigned int i = 0; i < _faces_index.size(); ++i)
//	{
//		Point3D<float> vec1 = _vertices[_faces_index[i].I[0]];
//		Point3D<float> vec2 = _vertices[_faces_index[i].I[1]];
//		Point3D<float> vec3 = _vertices[_faces_index[i].I[2]];
//		
//		Point3D<float> edge1 = vec3 - vec1;
//		Point3D<float> edge2 = vec2 - vec1;
//		Point3D<float> normal = Point3D<float>
//			(edge1.coords[1] * edge2.coords[2] - edge1.coords[2] * edge2.coords[1],
//				edge1.coords[2] * edge2.coords[0] - edge1.coords[0] * edge2.coords[2],
//				edge1.coords[0] * edge2.coords[1] - edge1.coords[1] * edge2.coords[0]) * (-1.0f);
//		Point3D<float>::Normalize(normal);
//		// normal vector coordinates
//
//		ofs.write((char*)&normal[0], sizeof(float));
//		ofs.write((char*)&normal[1], sizeof(float));
//		ofs.write((char*)&normal[2], sizeof(float));
//
//
//		//vertex 1 coordinates
//		ofs.write((char*)&vec1[0], sizeof(float));
//		ofs.write((char*)&vec1[1], sizeof(float));
//		ofs.write((char*)&vec1[2], sizeof(float));
//
//		//vertex 2 coordinates
//		ofs.write((char*)&vec2[0], sizeof(float));
//		ofs.write((char*)&vec2[1], sizeof(float));
//		ofs.write((char*)&vec2[2], sizeof(float));
//
//		//vertex 3 coordinates
//		ofs.write((char*)&vec3[0], sizeof(float));
//		ofs.write((char*)&vec3[1], sizeof(float));
//		ofs.write((char*)&vec3[2], sizeof(float));
//
//		//Attibutes
//		ofs.write(attribute, sizeof(short));
//	}
//
//	ofs.close();
//	std::cout << "Successfully export stl file." << std::endl;
//}
//
////parse the vector
//template<typename T>
//Point3D<T> MeshBuilder<T>::parse_point(std::ifstream& s)
//{
//	Point3D<T> vec;
//	vec[0] = parse_float(s);
//	vec[1] = parse_float(s);
//	vec[2] = parse_float(s);
//
//	return  vec;
//};
//
////parse the binary data and convert to float 
//template<typename T>
//float MeshBuilder<T>::parse_float(std::ifstream& s)
//{
//	char f_buf[sizeof(float)];
//	s.read(f_buf, 4);
//	float* fptr = (float*)f_buf;
//
//	return *fptr;
//};


#endif // !MESHBUILDER_H
