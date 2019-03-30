#include "DataInit.h"
#include <iostream>
#include <algorithm>

#include "Matrix.h"


template<typename T>
DataUnifier<T>::DataUnifier(std::vector<Point3D<T>> _vertices, std::vector<Point3D<T>> _normal, int _PIXEL)
{
	vertices = _vertices;
	normal = _normal;
	grid.nx = 0;
	grid.ny = 0;
	grid.nz = 0;
	PIXEL = _PIXEL;
}

template<typename T>
DataUnifier<T>::~DataUnifier()
{
}

template<typename T>
gridsize DataUnifier<T>::unify_data(T scale, int box)
{
	int i;
	T minx,maxx,miny,maxy,minz,maxz;
	
	minx=maxx=vertices[0][0];
	miny=maxy=vertices[0][1];
	minz=maxz=vertices[0][2];
	
	Point3D<T> tmp;
	int v_num = vertices.size();
	for(i=1; i<v_num; i++){ // find the max and min coord
		tmp=vertices[i];
		if(minx>tmp[0]) minx=tmp[0];
		if(maxx<tmp[0]) maxx=tmp[0];
		if(miny>tmp[1]) miny=tmp[1];
		if(maxy<tmp[1]) maxy=tmp[1];
		if(minz>tmp[2]) minz=tmp[2];
		if(maxz<tmp[2]) maxz=tmp[2];
	}
	
	T midx=0.5*(minx+maxx), midy=0.5*(miny+maxy), midz=0.5*(minz+maxz);
	for(i=0; i<v_num; i++){
		vertices[i][0]-=midx;	
		vertices[i][1]-=midy;	
		vertices[i][2]-=midz;
	}

	T r;
	T widx=maxx-minx, widy=maxy-miny, widz=maxz-minz;
	r = std::max(widx, std::max(widy, widz));
	int XPIXEL = 2*int(0.5*(PIXEL*widx/r +1))+1;
	int YPIXEL = 2*int(0.5*(PIXEL*widy/r +1))+1;
	int ZPIXEL = 2*int(0.5*(PIXEL*widz/r +1))+1;

	r = scale*r;
	for(i=0;i<v_num;i++)	
		vertices[i]/=r;

	maxx -= midx;	minx -= midx;
	maxy -= midy;  miny -= midy; 
	maxz -= midz;	minz -= midz;

	T xmin0 = 0.5*(XPIXEL-1)/PIXEL;
	T ymin0 = 0.5*(YPIXEL-1)/PIXEL;
	T zmin0 = 0.5*(ZPIXEL-1)/PIXEL;

	i = box; // calculate the boundingbox
	int vxmin = std::max(floor(PIXEL*(minx/r+xmin0))-i-1, (T)0.0);
	int vxmax = std::min(ceil(PIXEL*(maxx/r+xmin0))+i-1, (T)XPIXEL);
	int vymin = std::max(floor(PIXEL*(miny/r+ymin0))-i-1, (T)0.0);
	int vymax = std::min(ceil(PIXEL*(maxy/r+ymin0))+i-1, (T)YPIXEL);
	int vzmin = std::max(floor(PIXEL*(minz/r+zmin0))-i-1, (T)0.0);
	int vzmax = std::min(ceil(PIXEL*(maxz/r+zmin0))+i-1, (T)ZPIXEL);

	grid.nx = XPIXEL;
	grid.ny = YPIXEL;
	grid.nz = ZPIXEL;

	return grid;

//	cout<<"Grid Size: "<<XPIXEL<<"x"<<YPIXEL<<"x"<<ZPIXEL;


}

template<typename T>
int DataUnifier<T>::simplify_data()
{
	if (0 == grid.nx || 0 == grid.ny || 0 == grid.nz)
	{
		std::cout << "please unify the data first !" << std::endl;
		exit(-1);
	}

	int*** num = matrix3D<int>(grid.nx, grid.ny, grid.nz);
	int*** index = matrix3D<int>(grid.nx, grid.ny, grid.nz);

	int i(0), j(0), k(0), t(0), local(0), pt_count(0);

	T xmin = 0.5*(grid.nx-1)/PIXEL;
	T ymin = 0.5*(grid.ny-1)/PIXEL;
	T zmin = 0.5*(grid.nz-1)/PIXEL;
	T r;

	Point3D<T> tmp, tmp1;
	std::vector<Point3D<T>> ver, nor;
	ver.clear();	nor.clear();

	int v_num = vertices.size();

	for(t=0 ;t<v_num; t++){
		tmp = vertices[t];
		tmp1 = normal[t];

		i = (int)(floor((tmp[0]+xmin)*PIXEL));	
		j = (int)(floor((tmp[1]+ymin)*PIXEL));	
		k = (int)(floor((tmp[2]+zmin)*PIXEL));

		if(num[i][j][k]){
			num[i][j][k]++;
			r = 1./num[i][j][k];

			local = index[i][j][k];
			ver[local][0] = (1-r) * ver[local][0] + r * tmp[0];	
			ver[local][1] = (1-r) * ver[local][1] + r * tmp[1];
			ver[local][2] = (1-r) * ver[local][2] + r * tmp[2];

			nor[local][0] = (1 - r) * nor[local][0] + r * tmp1[0];
			nor[local][1] = (1 - r) * nor[local][1] + r * tmp1[1];
			nor[local][2] = (1 - r) * nor[local][2] + r * tmp1[2];
		}
		else{
			num[i][j][k]++;
			index[i][j][k] = pt_count;
			pt_count++;

			ver.push_back(tmp);
			nor.push_back(tmp1);
		}
	}

	v_num = pt_count;
	
	vertices.clear();
	normal.clear();
	
	vertices = ver;
	ver.clear();
	normal = nor;		
	nor.clear();

//	cout<<endl<<v_num<<" points after simplification\n";

	// TODO: it seems to something wrong
	for(t=0; t<v_num; t++){ // normalize the normal
		//ormal[t] = Point3D<T>::Normalize(normal[t]);
		Point3D<T>::Normalize(normal[t]);

	}

	delete_matrix3D<int>(num);	
	delete_matrix3D<int>(index);

	return v_num;
}

template<typename T>
std::vector<Point3D<T>> DataUnifier<T>::get_simple_vertices()
{
	return vertices;
}

template<typename T>
std::vector<Point3D<T>> DataUnifier<T>::get_simple_normal()
{
	return normal;
}


template DataUnifier<float>::DataUnifier();
template DataUnifier<float>::DataUnifier(std::vector<Point3D<float>> , std::vector<Point3D<float>> , int );
template DataUnifier<float>::~DataUnifier();
template int DataUnifier<float>::simplify_data();
template gridsize DataUnifier<float>::unify_data(float, int);
template std::vector<Point3D<float>> DataUnifier<float>::get_simple_vertices();
template std::vector<Point3D<float>> DataUnifier<float>::get_simple_normal();




template DataUnifier<double>::DataUnifier();
template DataUnifier<double>::DataUnifier(std::vector<Point3D<double>>, std::vector<Point3D<double>>, int);
template DataUnifier<double>::~DataUnifier();
template int DataUnifier<double>::simplify_data();
template gridsize DataUnifier<double>::unify_data(double, int);
template std::vector<Point3D<double>> DataUnifier<double>::get_simple_vertices();
template std::vector<Point3D<double>> DataUnifier<double>::get_simple_normal();

