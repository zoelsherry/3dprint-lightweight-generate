#ifndef STLSTREAM_H
#define STLSTREAM_H

#include "Geometry.h"
#include "Mesh.h"
#include <map>
#include <cassert>
#include <iostream>

template<typename T>
class STLMeshInStream
{
public:
	STLMeshInStream();
	STLMeshInStream(const std::string& _filename);
	~STLMeshInStream();

	void AddTriangle(Point3D<T> p0, Point3D<T> p1, Point3D<T> p2);
	Mesh<T> getMesh();
	void Clear();

	void load_stl();
	void ReadASCII();
	void ReadBinary();
	Point3D<T> parse_point(std::ifstream& s);
	float parse_float(std::ifstream& s);

	void unity_normal();

private:
	const std::string filename;
	Mesh<T> m_mesh;
	//std::unordered_map<Point3D<double>, int> m_hashmap;
	std::map<Point3D<T>, int> m_map;

	//
	int vertex_count = 0;
	//
};


template<typename T>
class STLMeshOutStream
{
public:
	STLMeshOutStream();
	STLMeshOutStream(const std::string& _filename);
	~STLMeshOutStream();


	void export_stl(std::vector<STLtriMeshFace<T>> stltrimesh);
	void export_stl(std::vector<triMeshFace<T>> trimesh);
	void export_stl(std::vector<Point3D<T>> _vertices, std::vector<TriFaceIndex> _faces_index);
	

private:
	const std::string filename;

};

///////////////////////////////////////////////////////////////////////
//STLMeshInStream                                                    //
///////////////////////////////////////////////////////////////////////
template<typename T>
STLMeshInStream<T>::STLMeshInStream()
{
}

template<typename T>
STLMeshInStream<T>::STLMeshInStream(const std::string& _filename) :filename("./" + _filename + ".stl")
{
	load_stl();
	unity_normal();
}

template<typename T>
STLMeshInStream<T>::~STLMeshInStream()
{
}

template<typename T>
void STLMeshInStream<T>::AddTriangle(Point3D<T> p0, Point3D<T> p1, Point3D<T> p2)
{
	int p0idx;
	int p1idx;
	int p2idx;
	int index;
	//std::unordered_map<Point3D<double>, int>::iterator got = m_hashmap.find(p0);
	typename std::map<Point3D<T>, int>::iterator got = m_map.find(p0);
	if (got == m_map.end())
	{
		p0idx = m_mesh.AddVertex(p0);
		m_map.emplace(p0, p0idx);
	}
	else
	{
		index = m_map.at(p0);
		assert(index >= 0);
		p0idx = index;
	}

	got = m_map.find(p1);
	if (got == m_map.end())
	{
		p1idx = m_mesh.AddVertex(p1);
		m_map.emplace(p1, p1idx);
	}
	else
	{
		index = m_map.at(p1);
		assert(index >= 0);
		p1idx = index;
	}

	got = m_map.find(p2);
	if (got == m_map.end())
	{
		p2idx = m_mesh.AddVertex(p2);
		m_map.emplace(p2, p2idx);
	}
	else
	{
		index = m_map.at(p2);
		assert(index >= 0);
		p2idx = index;
	}

	//cout << m_hashmap.size() << endl;
	TriFaceIndex tri(p0idx, p1idx, p2idx);
	m_mesh.AddFace(tri);
}

template<typename T>
Mesh<T> STLMeshInStream<T>::getMesh()
{
	return m_mesh;
}

template<typename T>
void STLMeshInStream<T>::Clear()
{
	m_map.clear();
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
template<typename T>
void STLMeshInStream<T>::load_stl()
{
	std::ifstream fin(filename, std::ios::in);
	if (!fin.good())
	{
		std::cout << "Error: cannot open the STL file " << filename << std::endl;
		return;
	}

	std::string headStr;

	getline(fin, headStr, ' ');
	fin.close();

	if (headStr.empty())
	{
		std::cout << "Error: the STL file is wrong " << filename << std::endl;
		return;
	}

	if (headStr[0] == 's')
	{
		ReadASCII();
	}
	else
	{
		ReadBinary();
	}
}

template<typename T>
void STLMeshInStream<T>::ReadASCII()
{

}

template<typename T>
void STLMeshInStream<T>::ReadBinary()
{
	// read binary file
	// file open for reading
	// operations are performed in binary mode rather than text
	std::ifstream fin(filename, std::ios::in | std::ios::binary);
	if (!fin.good())
	{
		std::cout << "Error: cannot open the STL file " << filename << std::endl;
		return;
	}
	// header
	char header_info[80] = "";
	char n_triangles[4];
	fin.read(header_info, 80);
	fin.read(n_triangles, 4);
	unsigned int* r = (unsigned int*)n_triangles;
	unsigned int num_triangles = *r;
	// number of triangles
	std::cout << "Num Triangles: " << num_triangles << "\n";
	for (unsigned int i = 0; i < num_triangles; i++)
	{
		// normal vector of triangle
		auto normal = parse_point(fin);
		// vertice of triangle
		auto v1 = parse_point(fin);
		auto v2 = parse_point(fin);
		auto v3 = parse_point(fin);
		// push into triangle list
		m_mesh.m_stltrimesh.push_back(STLtriMeshFace<T>(normal, v1, v2, v3));
		// push into noraml list
		AddTriangle(v1, v2, v3);
		// attributes
		char dummy[2];

		fin.read(dummy, 2);
	}
	std::cout << "Successfully import stl file" << "\n";
}


//parse the vector
template<typename T>
Point3D<T> STLMeshInStream<T>::parse_point(std::ifstream& s)
{
	Point3D<T> vec;
	vec[0] = (T)parse_float(s);
	vec[1] = (T)parse_float(s);
	vec[2] = (T)parse_float(s);

	return  vec;
};

//parse the binary data and convert to 
template<typename T>
float STLMeshInStream<T>::parse_float(std::ifstream& s)
{
	char f_buf[sizeof(float)];
	s.read(f_buf, 4);
	float* fptr = (float*)f_buf;

	return *fptr;
};


template<typename T>
void STLMeshInStream<T>::unity_normal()
{
	m_mesh.m_normal.resize(m_mesh.m_vertices.size());
	int i = 0;
	for (auto p : m_mesh.m_faces_index)
	{
		m_mesh.m_normal[p.I[0]] += m_mesh.m_stltrimesh[i].normal;
		m_mesh.m_normal[p.I[1]] += m_mesh.m_stltrimesh[i].normal;
		m_mesh.m_normal[p.I[2]] += m_mesh.m_stltrimesh[i].normal;

		i++;
	}

	for (auto &p : m_mesh.m_normal)
	{
		Point3D<T>::Normalize(p);
	}
}



///////////////////////////////////////////////////////////////////////
//STLMeshInStream                                                    //
///////////////////////////////////////////////////////////////////////

template<typename T>
STLMeshOutStream<T>::STLMeshOutStream()
{

}

template<typename T>
STLMeshOutStream<T>::STLMeshOutStream(const std::string& _filename) :filename("./" + _filename + ".stl")
{

}

template<typename T>
STLMeshOutStream<T>::~STLMeshOutStream()
{
}

template<typename T>
void STLMeshOutStream<T>::export_stl(std::vector<STLtriMeshFace<T>> stltrimesh)
{
	//binary file
	std::string header_info = "solid " + filename + "-output";
	char head[80];
	std::strncpy(head, header_info.c_str(), sizeof(head));
	//attributes
	char attribute[2] = "0";
	//number of triangles
	int nTrilong = stltrimesh.size();
	std::ofstream ofs(filename, std::ofstream::binary | std::ofstream::trunc);
	if (!ofs.good())
	{
		std::cout << "cannot open output file" << filename << std::endl;
	}
	//write header into stl file
	ofs.write(head, sizeof(head));
	//write number of triangles into stl file
	ofs.write((char*)&nTrilong, sizeof(int));
	for (unsigned int i = 0; i < stltrimesh.size(); ++i)
	{

		// normal vector coordinates

		ofs.write((char*)&(float)stltrimesh[i].normal[0], sizeof(float));
		ofs.write((char*)&(float)stltrimesh[i].normal[1], sizeof(float));
		ofs.write((char*)&(float)stltrimesh[i].normal[2], sizeof(float));


		//vertex 1 coordinates
		ofs.write((char*)&(float)stltrimesh[i].vec1[0], sizeof(float));
		ofs.write((char*)&(float)stltrimesh[i].vec1[1], sizeof(float));
		ofs.write((char*)&(float)stltrimesh[i].vec1[2], sizeof(float));

		//vertex 2 coordinates
		ofs.write((char*)&(float)stltrimesh[i].vec2[0], sizeof(float));
		ofs.write((char*)&(float)stltrimesh[i].vec2[1], sizeof(float));
		ofs.write((char*)&(float)stltrimesh[i].vec2[2], sizeof(float));

		//vertex 3 coordinates
		ofs.write((char*)&(float)stltrimesh[i].vec3[0], sizeof(float));
		ofs.write((char*)&(float)stltrimesh[i].vec3[1], sizeof(float));
		ofs.write((char*)&(float)stltrimesh[i].vec3[2], sizeof(float));

		//Attibutes
		ofs.write(attribute, sizeof(short));
	}

	ofs.close();
	std::cout << "Successfully export stl file." << std::endl;
}

template<typename T>
void STLMeshOutStream<T>::export_stl(std::vector<triMeshFace<T>> trimesh)
{
	//binary file
	std::string header_info = "solid " + filename + "-output";
	char head[80];
	std::strncpy(head, header_info.c_str(), sizeof(head));
	//attributes
	char attribute[2] = "0";
	//number of triangles
	int nTrilong = trimesh.size();
	std::ofstream ofs(filename, std::ofstream::binary | std::ofstream::trunc);
	if (!ofs.good())
	{
		std::cout << "cannot open output file" << filename << std::endl;
	}
	//write header into stl file
	ofs.write(head, sizeof(head));
	//write number of triangles into stl file
	ofs.write((char*)&nTrilong, sizeof(int));
	for (unsigned int i = 0; i < trimesh.size(); ++i)
	{
		Point3D<float> vec1 = (Point3D<float>)trimesh[i].vec1;
		Point3D<float> vec2 = (Point3D<float>)trimesh[i].vec2;
		Point3D<float> vec3 = (Point3D<float>)trimesh[i].vec3;

		Point3D<float> edge1 = vec3 - vec1;
		Point3D<float> edge2 = vec2 - vec1;
		Point3D<float> normal = Point3D<float>
			(edge1.coords[1] * edge2.coords[2] - edge1.coords[2] * edge2.coords[1],
				edge1.coords[2] * edge2.coords[0] - edge1.coords[0] * edge2.coords[2],
				edge1.coords[0] * edge2.coords[1] - edge1.coords[1] * edge2.coords[0]) * (-1.0f);
		Point3D<float>::Normalize(normal);
		// normal vector coordinates

		ofs.write((char*)&normal[0], sizeof(float));
		ofs.write((char*)&normal[1], sizeof(float));
		ofs.write((char*)&normal[2], sizeof(float));


		//vertex 1 coordinates
		ofs.write((char*)&vec1[0], sizeof(float));
		ofs.write((char*)&vec1[1], sizeof(float));
		ofs.write((char*)&vec1[2], sizeof(float));

		//vertex 2 coordinates
		ofs.write((char*)&vec2[0], sizeof(float));
		ofs.write((char*)&vec2[1], sizeof(float));
		ofs.write((char*)&vec2[2], sizeof(float));

		//vertex 3 coordinates
		ofs.write((char*)&vec3[0], sizeof(float));
		ofs.write((char*)&vec3[1], sizeof(float));
		ofs.write((char*)&vec3[2], sizeof(float));

		//Attibutes
		ofs.write(attribute, sizeof(short));
	}

	ofs.close();
	std::cout << "Successfully export stl file." << std::endl;
}

template<typename T>
void STLMeshOutStream<T>::export_stl(std::vector<Point3D<T>> _vertices, std::vector<TriFaceIndex> _faces_index)
{
	//binary file
	std::string header_info = "solid " + filename + "-output";
	char head[80];
	std::strncpy(head, header_info.c_str(), sizeof(head));
	//attributes
	char attribute[2] = "0";
	//number of triangles
	int nTrilong = _faces_index.size();
	std::ofstream ofs(filename, std::ofstream::binary | std::ofstream::trunc);
	if (!ofs.good())
	{
		std::cout << "cannot open output file" << filename << std::endl;
	}
	//write header into stl file
	ofs.write(head, sizeof(head));
	//write number of triangles into stl file
	ofs.write((char*)&nTrilong, sizeof(int));
	for (unsigned int i = 0; i < _faces_index.size(); ++i)
	{
		Point3D<float> vec1 = (Point3D<float>)_vertices[_faces_index[i].I[0]];
		Point3D<float> vec2 = (Point3D<float>)_vertices[_faces_index[i].I[1]];
		Point3D<float> vec3 = (Point3D<float>)_vertices[_faces_index[i].I[2]];

		Point3D<float> edge1 = vec3 - vec1;
		Point3D<float> edge2 = vec2 - vec1;
		Point3D<float> normal = Point3D<float>
			(edge1.coords[1] * edge2.coords[2] - edge1.coords[2] * edge2.coords[1],
				edge1.coords[2] * edge2.coords[0] - edge1.coords[0] * edge2.coords[2],
				edge1.coords[0] * edge2.coords[1] - edge1.coords[1] * edge2.coords[0]) * (-1.0f);
		Point3D<float>::Normalize(normal);
		// normal vector coordinates

		ofs.write((char*)&normal[0], sizeof(float));
		ofs.write((char*)&normal[1], sizeof(float));
		ofs.write((char*)&normal[2], sizeof(float));


		//vertex 1 coordinates
		ofs.write((char*)&vec1[0], sizeof(float));
		ofs.write((char*)&vec1[1], sizeof(float));
		ofs.write((char*)&vec1[2], sizeof(float));

		//vertex 2 coordinates
		ofs.write((char*)&vec2[0], sizeof(float));
		ofs.write((char*)&vec2[1], sizeof(float));
		ofs.write((char*)&vec2[2], sizeof(float));

		//vertex 3 coordinates
		ofs.write((char*)&vec3[0], sizeof(float));
		ofs.write((char*)&vec3[1], sizeof(float));
		ofs.write((char*)&vec3[2], sizeof(float));

		//Attibutes
		ofs.write(attribute, sizeof(short));
	}

	ofs.close();
	std::cout << "Successfully export stl file." << std::endl;
}


#endif // !STLSTREAM_H

