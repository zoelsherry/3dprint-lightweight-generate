/*
Copyright (c) 2018, Joy Smith (Zhengyuan Shi)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Xi'an Jiaotong University nor the names of its contributors
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

#include "DataInit.h"
#include "PointStream.h"
#include "SurfaceProcess.h"
#include "Matrix.h"
#include "splinterp.h"
#include "RangeImpl.h"
#include "Isosurface.h"
#include "distance_3D.h"
#include "StructureFiller.h"
#include "StlStream.h"
#include <vector>

#undef DEBUG

using namespace std;


typedef OrientedPointInStream< double > PointInStream;


void main(int argc, char** argv)
{
	clock_t start,stop;	
	double time;

	//int PIXEL=128;
	int PIXEL = 200;

	float scale=1.1;
	
	//cin>>PIXEL;

//////////////////////////////////////////////////////////////////////////////////

	
	cout << "\nReading and simplifying data ... ";	
	start = clock();	  
	
	// Read .ply data
	//PointInStream* point_stream;
	//point_stream = new PLYOrientedPointInStream<double>("dragon");
	//std::vector<OrientedPoint3D<double>> point_data;
	//point_data = point_stream->getPointData();

	//vector<Point3D<double>> vertices;
	//vector<Point3D<double>> normal;
	//for (OrientedPoint3D<double> p_ : point_data)
	//{
	//	Point3D<double> p(p_.p[0], p_.p[1], p_.p[2]);
	//	Point3D<double> n(p_.n[0], p_.n[1], p_.n[2]);
	//	vertices.push_back(p);
	//	normal.push_back(n);
	//}

	// read .stl data
	STLMeshInStream<double>* stl_istream = new STLMeshInStream<double>("bunny");
	Mesh<double> mesh_in = stl_istream->getMesh();

	vector<Point3D<double>> vertices = mesh_in.m_vertices;
	vector<Point3D<double>> normal = mesh_in.m_normal;

	//
	cout << vertices.size() << endl;
	cout << normal.size() << endl;
	//

	int num_original=vertices.size();
	int v_num = vertices.size();

	DataUnifier<double> initial_data(vertices, normal, PIXEL);
	gridsize grid = initial_data.unify_data(scale, 1); // scale the data and return the gridsize
	v_num = initial_data.simplify_data(); // simplify the data
	vertices = initial_data.get_simple_vertices();
	normal = initial_data.get_simple_normal();


	stop = clock();	
	time = (double(stop)-double(start))/CLOCKS_PER_SEC;
	cout << "done " << time << " sec(s)." << endl;
	cout<<"there are "<<num_original<<" points"<<endl;
	cout<<v_num<<" points after simplification\n";

#ifdef DEBUG
	std::ofstream out_vertices;
	out_vertices.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try
	{
		out_vertices.open("./outdata/vertices_simply.txt");
		for (int i = 0; i < vertices.size(); i++)
			out_vertices << vertices[i][0] << "," << vertices[i][1] << "," << vertices[i][2] << std::endl;

		out_vertices.close();
	}
	catch (std::ofstream::failure e)
	{
		std::cout << "ERROR::SIGNEDDIST::FILE FAIL TO WRITE!" << std::endl;
	}

	std::ofstream out_normal;
	out_normal.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try
	{
		out_normal.open("./outdata/normal_simply.txt");
		for (int i = 0; i < normal.size(); i++)
			out_normal << normal[i][0] << "," << normal[i][1] << "," << normal[i][2] << std::endl;

		out_normal.close();
	}
	catch (std::ofstream::failure e)
	{
		std::cout << "ERROR::SIGNEDDIST::FILE FAIL TO WRITE!" << std::endl;
	}
#endif // DEBUG

//////////////////////////////////////////////////////////////////////////////////
	cout << "\nGrid Size: " << grid.nx << "x" << grid.ny << "x" << grid.nz << endl;

	cout << "\nComputing Signed Distance Function ... ";	
	start = clock();
	int*** idx = matrix3D<int>(grid.nx + 2, grid.ny + 2, grid.nz + 2);
	double*** dis0 = matrix3D<double>(grid.nx, grid.ny, grid.nz);
	double*** disd = matrix3D<double>(grid.nx + 4, grid.ny + 4, grid.nz + 4);
	bool*** lock = matrix3D<bool>(grid.nx + 2, grid.ny + 2, grid.nz + 2);
	
	SignedDistSolver<double> sdist(grid, vertices, normal, PIXEL);
	sdist.initial_dist(dis0, idx, lock);
	sdist.distance_sweeping(dis0, disd, idx, lock);
	sdist.unify(disd, 2);

	double*** phi = matrix3D<double>(grid.nx, grid.ny, grid.nz);
	sdist.calc_signed_dist(dis0, idx, phi);


	delete_matrix3D<bool>(lock);
	delete_matrix3D<double>(disd);
	delete_matrix3D<double>(dis0);
	delete_matrix3D<int>(idx);

	
	stop = clock();
	time = (double(stop)-double(start))/CLOCKS_PER_SEC;
	cout << "done " << time << " sec(s)." << endl;

#ifdef DEBUG
	std::ofstream out_phi;
	out_phi.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try
	{
		out_phi.open("./outdata/phi_first.txt");
		for (int i = 0; i < grid.nx; i++)
			for (int j = 0; j < grid.ny; j++)
				for (int k = 0; k < grid.nz; k++)
					out_phi << phi[i][j][k] << std::endl;

		out_phi.close();
	}
	catch (std::ofstream::failure e)
	{
		std::cout << "ERROR::SIGNEDDIST::FILE FAIL TO WRITE!" << std::endl;
	}
#endif // DEBUG

//////////////////////////////////////////////////////////////////////////////////


	cout << "Surface computing is running ......" << endl;
	start = clock();

	// TODO: the maxiter is to be changed
	SmoothSurface<double>(grid, phi, 50);

	stop = clock();
	time = (double(stop) - double(start)) / CLOCKS_PER_SEC;
	cout << "done! " << time << " sec(s)." << endl;

#ifdef DEBUG
	std::ofstream out_phi_smooth;
	out_phi_smooth.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try
	{
		out_phi_smooth.open("./outdata/phi_smooth.txt");
		for (int i = 0; i < grid.nx; i++)
			for (int j = 0; j < grid.ny; j++)
				for (int k = 0; k < grid.nz; k++)
					out_phi_smooth << phi[i][j][k] << std::endl;

		out_phi_smooth.close();
	}
	catch (std::ofstream::failure e)
	{
		std::cout << "ERROR::SIGNEDDIST::FILE FAIL TO WRITE!" << std::endl;
	}
#endif // DEBUG


////////////////////////////////////////////////////////////////////////////////////
// fill lightweight structure													  //
////////////////////////////////////////////////////////////////////////////////////


	int flag0 = 1;
	int flag = 0;
	int nphi_num = 0; // cout the size of lightweight structure
	vector<double> nphi_v;
	while (flag0)
	{
		cout << "please input the function (e.g. P,D,G,F,I): ";
		string funtype;
		cin >> funtype;
		cout << "please input the volumn fraction (e.g. 0.05:0.05:0.95): ";
		double f_volfra;
		cin >> f_volfra;
		int i_volfra = 100 * f_volfra;

		if (i_volfra > 50)
		{
			i_volfra = 100 - i_volfra;
			flag = 1;
		}
		string s_volfra = to_string(i_volfra);
		if (i_volfra < 10)
		{
			s_volfra = "0" + s_volfra;
		}

		string tag = funtype + s_volfra;
		//
		cout << tag << endl;
		//
		ifstream fp;
		fp.exceptions(std::ifstream::badbit); // delete the `failbit` otherwise error
		try
		{
			fp.open("various_surface/" + tag + ".m", ifstream::in);
			double a;
			while (fp >> a)
			{
				if (flag == 1)
				{
					a = 1 - a;
				}
				nphi_v.push_back(a);
			}
		}
		catch (ifstream::failure e)
		{
			cout << "Caught an exception: " << e.what() << endl;
			cout << "ERROR::SURFACE::The data source does not exist! Please input again! " << endl;
		}
		flag0 = 0;
	}

	for (auto &p : nphi_v)
	{
		p = p * 2 - 1;
	}
	//
	//cout << nphi_v.back() << endl;
	//cout << nphi_v.size() << endl;
	//

	double*** phi2 = matrix3D<double>(grid.nx, grid.ny, grid.nz);
	for (int i = 0; i < grid.nx; i++)
		for (int j = 0; j < grid.ny; j++)
			for (int k = 0; k < grid.nz; k++)
			{
				phi2[i][j][k] = phi[i][j][k];

			}
	start = clock();
	StructureFiller<double> struct_fill(nphi_v, 8);
	struct_fill.generate_whole(grid, phi);
	struct_fill.generate_inner(grid, phi2);

	cout << "Computing is finished!" << endl;
	stop = clock();
	time = (double(stop) - double(start)) / CLOCKS_PER_SEC;
	cout << "done! " << time << " sec(s)." << endl;

#ifdef DEBUG
	std::ofstream out_phi_structure;
	out_phi_structure.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try
	{
		out_phi_structure.open("./outdata/phi_structure.txt");
		for (int i = 0; i < grid.nx; i++)
			for (int j = 0; j < grid.ny; j++)
				for (int k = 0; k < grid.nz; k++)
					out_phi_structure << phi[i][j][k] << std::endl;

		out_phi_structure.close();
	}
	catch (std::ofstream::failure e)
	{
		std::cout << "ERROR::SIGNEDDIST::FILE FAIL TO WRITE!" << std::endl;
	}
#endif // DEBUG
	
//////////////////////////////////////////////////////////////////////////////////////
//marchingcubes to extract the isosurface										    //
//////////////////////////////////////////////////////////////////////////////////////

	//export .ply

	//start = clock();
	//cout << "Extract whole isosurface ......" << endl;
	//Isosurface<double> surface_w(grid.nx, grid.ny, grid.nz);
	//Mesh<double> mesh = surface_w.GenerateIsosurface(phi);
	//std::vector<PlyFace> ply_face;
	//PlyFace _face;
	//for (int i = 0; i < mesh.m_faces_index.size(); i++)
	//{
	//	_face.vertices = new int[3];
	//	_face.nr_vertices = int(3);
	//	for (int j = 0; j<_face.nr_vertices; j++)
	//		_face.vertices[j] = mesh.m_faces_index[i].I[j];
	//	ply_face.push_back(_face);
	//}
	//PointOutStream<double>* point_ostream;
	//point_ostream = new PointOutStream<double>("dragon_whole_out");
	//point_ostream->write_to_ply(mesh.m_vertices, ply_face);

	//stop = clock();
	//time = (double(stop) - double(start)) / CLOCKS_PER_SEC;
	//cout << "done! " << time << " sec(s)." << endl;



	//start = clock();
	//cout << "Extract inner isosurface ......" << endl;
	//Isosurface<double> surface_in(grid.nx, grid.ny, grid.nz);
	//Mesh<double> mesh2 = surface_in.GenerateIsosurface(phi2);
	//std::vector<PlyFace> ply_face2;
	//PlyFace _face2;
	//for (int i = 0; i < mesh2.m_faces_index.size(); i++)
	//{
	//	_face2.vertices = new int[3];
	//	_face2.nr_vertices = int(3);
	//	for (int j = 0; j < _face2.nr_vertices; j++)
	//		_face2.vertices[j] = mesh2.m_faces_index[i].I[j];
	//	ply_face2.push_back(_face2);
	//}
	//PointOutStream<double>* point_ostream2;
	//point_ostream2 = new PointOutStream<double>("dragon_inner_out");
	//point_ostream2->write_to_ply(mesh2.m_vertices, ply_face2);

	//stop = clock();
	//time = (double(stop) - double(start)) / CLOCKS_PER_SEC;
	//cout << "done! " << time << " sec(s)." << endl;

	// export .stl
	start = clock();
	cout << "Extract whole isosurface ......" << endl;
	Isosurface<double> surface_w(grid.nx, grid.ny, grid.nz);
	Mesh<double> mesh = surface_w.GenerateSTLIsosurface(phi);
	
	STLMeshOutStream<double>* stl_ostream = new STLMeshOutStream<double>("bunny_whole_out");
	stl_ostream->export_stl(mesh.m_trimesh);

	stop = clock();
	time = (double(stop) - double(start)) / CLOCKS_PER_SEC;
	cout << "done! " << time << " sec(s)." << endl;



	start = clock();
	cout << "Extract inner isosurface ......" << endl;
	Isosurface<double> surface_in(grid.nx, grid.ny, grid.nz);
	Mesh<double> mesh2 = surface_in.GenerateSTLIsosurface(phi2);
	
	STLMeshOutStream<double>* stl_ostream2 = new STLMeshOutStream<double>("bunny_inner_out");
	stl_ostream2->export_stl(mesh2.m_trimesh);

	stop = clock();
	time = (double(stop) - double(start)) / CLOCKS_PER_SEC;
	cout << "done! " << time << " sec(s)." << endl;
	system("pause");
}