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
prior writften permission. 

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

#ifndef POINT_STREAM_INCLUDED
#define POINT_STREAM_INCLUDED
#include <vector>
#include <iostream>
#include "Ply.h"
#include "Geometry.h"

#undef DEBUG
//#define DEBUG
//#undef WRITE
#define WRITE

template<typename Real>
struct boundingbox
{
	Real max_x, min_x;
	Real max_y, min_y;
	Real max_z, min_z;
};

template< class Real >
class OrientedPointInStream
{

public:
	virtual ~OrientedPointInStream( void ){}
	virtual void reset( void ) = 0;
	virtual bool nextPoint(OrientedPoint3D<Real> &pointDate) = 0;
	virtual boundingbox<Real> getboundingBox() = 0;
	virtual std::vector<OrientedPoint3D<Real>> getPointData() = 0;

};


template< class Real >
class PLYOrientedPointInStream : public OrientedPointInStream< Real >
{
	char* _fileName;
	PlyFile* _ply;
	int _nr_elems;
	char **_elist;

	int _pCount, _pIdx;
	void _free(void);

	
	std::vector<OrientedPoint3D<Real>> point_data;


public:
	PLYOrientedPointInStream(const char* fileName);
	~PLYOrientedPointInStream(void);
	void reset(void);
	bool nextPoint(OrientedPoint3D<Real> &pointDate);
	std::vector<OrientedPoint3D<Real>> getPointData();
	boundingbox<Real> getboundingBox();

};


//////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef WRITE
template< class Real >
class PointOutStream
{
	char* _fileName;
	PlyFile* _ply;

	int _pCount, _pIdx;

	void _free(void);

public:
	PointOutStream(const char* fileName);
	~PointOutStream(void);
	void write_to_ply(const std::vector<Point3D<Real>>& verts, const std::vector<PlyFace>& faces);

};
#endif //WRITE
//////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////
// PLYOrientedPointInStream //
////////////////////////////
template< class Real >
PLYOrientedPointInStream< Real >::PLYOrientedPointInStream(const char* fileName)
{
	_fileName = new char[strlen(fileName) + 1];
	strcpy(_fileName, fileName);
	_ply = NULL;
	reset();
}
template< class Real >
void PLYOrientedPointInStream< Real >::reset(void)
{
	int fileType;
	float version;
	PlyProperty** plist;
	if (_ply) _free();
	_ply = ply_open_for_reading(_fileName, &_nr_elems, &_elist, &fileType, &version);
	if (!_ply)
	{
		fprintf(stderr, "[ERROR] Failed to open ply file for reading: %s\n", _fileName);
		exit(0);
	}
	bool foundVertices = false;
	for (int i = 0; i<_nr_elems; i++)
	{
		int num_elems;
		int nr_props;
		char* elem_name = _elist[i];
		plist = ply_get_element_description(_ply, elem_name, &num_elems, &nr_props);
		if (!plist)
		{
			fprintf(stderr, "[ERROR] Failed to get element description: %s\n", elem_name);
			exit(0);
		}

		if (equal_strings("vertex", elem_name))
		{
			foundVertices = true;
			_pCount = num_elems, _pIdx = 0;
			for (int i = 0; i<PlyOrientedVertex< Real >::ReadComponents; i++)
				if (!ply_get_property(_ply, elem_name, &(PlyOrientedVertex< Real >::ReadProperties[i])))
				{
					fprintf(stderr, "[ERROR] Failed to find property in ply file: %s\n", PlyOrientedVertex< Real >::ReadProperties[i].name);
					exit(0);
				}


		}
		for (int j = 0; j<nr_props; j++)
		{
			free(plist[j]->name);
			free(plist[j]);
		}
		free(plist);
		if (foundVertices) break;
	}
	if (!foundVertices)
	{
		fprintf(stderr, "[ERROR] Could not find vertices in ply file\n");
		exit(0);
	}
}
template< class Real >
void PLYOrientedPointInStream< Real >::_free(void)
{
	if (_ply) ply_close(_ply), _ply = NULL;
	if (_elist)
	{
		for (int i = 0; i<_nr_elems; i++) free(_elist[i]);
		free(_elist);
	}
}
template< class Real >
PLYOrientedPointInStream< Real >::~PLYOrientedPointInStream(void)
{
	_free();
	if (_fileName) delete[] _fileName, _fileName = NULL;
}
template< class Real >
bool PLYOrientedPointInStream< Real >::nextPoint(OrientedPoint3D<Real> &pointDate)
{
	if (_pIdx<_pCount)
	{
		PlyOrientedVertex< Real > op;
		ply_get_element(_ply, (void *)&op);
		pointDate.p = op.point;
		pointDate.n = op.normal;
		_pIdx++;
#ifdef DEBUG
		/*std::cout << "vertex: " << std::endl;
		std::cout << pointCoord[i][0] << " " << pointCoord[i][1] << " " << pointCoord[i][2] << std::endl;
		std::cout << "normal: " << std::endl;
		std::cout << op.normal[0] << " " << op.normal[1] << " " << op.normal[2] << std::endl;*/
#endif	// DEBUG
		return true;
	}
	else return false;
}

template<class Real>
std::vector<OrientedPoint3D<Real>> PLYOrientedPointInStream<Real>::getPointData()
{
	OrientedPoint3D<Real> p;
	for (int i = 0; i < _pCount; i++)
	{
		nextPoint(p);
		point_data.push_back(p);
#ifdef DEBUG
		//std::cout << "vertex: " << std::endl;
		//std::cout << pointCoord[i][0] << " " << pointCoord[i][1] << " " << pointCoord[i][2] << std::endl;
#endif // DEBUG

	}
	return point_data;
}

template<class Real>
boundingbox<Real> PLYOrientedPointInStream<Real>::getboundingBox()
{
	std::vector<Real> nx;
	std::vector<Real> ny;
	std::vector<Real> nz;
	for (size_t i = 0; i < point_data.size(); i++)
	{
		nx.push_back(point_data[i].p[0]);
		ny.push_back(point_data[i].p[1]);
		nz.push_back(point_data[i].p[2]);
	}
	auto minmax_x = std::minmax_element(nx.begin(), nx.end());
	auto minmax_y = std::minmax_element(ny.begin(), ny.end());
	auto minmax_z = std::minmax_element(nz.begin(), nz.end());

	boundingbox<Real> bound_box;
	bound_box.min_x = *minmax_x.first;
	bound_box.max_x = *minmax_x.second;
	bound_box.min_y = *minmax_y.first;
	bound_box.max_y = *minmax_y.second;
	bound_box.min_z = *minmax_z.first;
	bound_box.max_z = *minmax_z.second;

	return bound_box;

}

//////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////
// PLYOrientedPointOutStream //
///////////////////////////////

#ifdef WRITE
template<class Real>
inline void PointOutStream<Real>::_free(void)
{
	if (_ply) ply_close(_ply), _ply = NULL;

}
template<class Real>
PointOutStream<Real>::PointOutStream(const char * fileName)
{
	_fileName = new char[strlen(fileName) + 1];
	strcpy(_fileName, fileName);
	_ply = NULL;
}

template< class Real >
PointOutStream<Real>::~PointOutStream(void)
{
	_free();
	if (_fileName) delete[] _fileName, _fileName = NULL;
}

template<class Real>
void PointOutStream<Real>::write_to_ply(const std::vector<Point3D<Real>>& verts, const std::vector<PlyFace>& faces)
{
	int fileType;
	if (_ply) _free();
	float version;
	const char *elem_names[] = { "vertex" , "face" };
	if (_ply) _free();
#if 1
	_ply = ply_open_for_writing(_fileName, 2, elem_names, PLY_ASCII, &version);
#else
	_ply = ply_open_for_writing(_fileName, 2, elem_names, PLY_BINARY_BE, &version);
#endif
	if (!_ply)
	{
		fprintf(stderr, "[ERROR] Failed to open ply file for writing: %s\n", _fileName);
		exit(0);
	}

	int nverts = verts.size();
	int nfaces = faces.size();

	int num_elems;
	int nr_props;
	/* describe what properties go into the vertex and face elements */
	ply_element_count(_ply, "vertex", nverts);
	for (int i = 0; i < PlyVertex< Real >::WriteComponents; i++)
		if (!ply_describe_property(_ply, "vertex", &(PlyVertex< Real >::ReadProperties[i])))
		{
			fprintf(stderr, "[ERROR] Failed to write property in ply file: %s\n", PlyOrientedVertex< Real >::ReadProperties[i].name);
			exit(0);
		}

	ply_element_count(_ply, "face", nfaces);
	ply_describe_property(_ply, "face", &face_props[0]);


	/* write a comment and an object information field */
	ply_put_comment(_ply, _strdup("author: NOPLAYER"));
	ply_put_obj_info(_ply, _strdup("random information"));

	/* we have described exactly what we will put in the file, so */
	/* we are now done with the header info */
	ply_header_complete(_ply);

	/* set up and write the vertex elements */
	ply_put_element_setup(_ply, "vertex");
	for (int i = 0; i < nverts; i++)
		ply_put_element(_ply, (void *)&verts[i]);

	/* set up and write the face elements */
	ply_put_element_setup(_ply, "face");
	for (int i = 0; i < nfaces; i++)
		ply_put_element(_ply, (void *)&faces[i]);

	/* close the PLY file */
	ply_close(_ply);
	
}

#endif // WRITE





#endif // POINT_STREAM_INCLUDED


