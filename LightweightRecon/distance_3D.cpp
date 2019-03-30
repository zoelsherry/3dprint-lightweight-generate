#include "distance_3D.h"
#include "Matrix.h"
#include <algorithm>
#include <limits>

//#define max
//#define min

template<typename T>
SignedDistSolver<T>::~SignedDistSolver()
{
	vertices.clear();
	normal.clear();
}

template<typename T>
void SignedDistSolver<T>::calc_signed_dist(T *** dis0, int *** idx, T *** phi)
{
	T*** ip = matrix3D<T>(XPIXEL + 4, YPIXEL + 4, ZPIXEL + 4);
	T x, y, z, r, xmin, ymin, zmin;
	xmin = 0.5*(XPIXEL - 1) / PIXEL;
	ymin = 0.5*(YPIXEL - 1) / PIXEL;
	zmin = 0.5*(ZPIXEL - 1) / PIXEL;

	T max_dist = 0;
	for (int i = 0; i < XPIXEL; i++) {
		for (int j = 0; j < YPIXEL; j++) {
			for (int k = 0; k < ZPIXEL; k++) {
				max_dist = std::max(max_dist, fabs(dis0[i][j][k]));
				x = i * 1. / PIXEL - xmin - vertices[idx[i + 1][j + 1][k + 1]][0];
				y = j * 1. / PIXEL - ymin - vertices[idx[i + 1][j + 1][k + 1]][1];
				z = k * 1. / PIXEL - zmin - vertices[idx[i + 1][j + 1][k + 1]][2];
				ip[i + 2][j + 2][k + 2] = (x*normal[idx[i + 1][j + 1][k + 1]][0] + y * normal[idx[i + 1][j + 1][k + 1]][1] + z * normal[idx[i + 1][j + 1][k + 1]][2]);
			}
		}
	}


	//max_dist *= 2;
	//for(i=2; i<XPIXEL+2; i++)
	//	for(j=2; j<YPIXEL+2; j++)
	//		for(k=2; k<ZPIXEL+2; k++)
	//			ip[i][j][k] = ip[i][j][k]/max_dist + 0.5;


	//////////////////////////////////////////////////////////////////////////////

	int*** si = matrix3D<int>(XPIXEL, YPIXEL, ZPIXEL);

	for (int i = 0; i < XPIXEL; i++)
		for (int j = 0; j < YPIXEL; j++)
			for (int k = 0; k < ZPIXEL; k++)
			{
				si[i][j][k] = sign_function(ip[i + 2][j + 2][k + 2]);
				// TODO: 3 should be adjusted
				phi[i][j][k] = si[i][j][k] * dis0[i][j][k] * 200.0f;
			}

	delete_matrix3D<T>(ip);	
	delete_matrix3D<int>(si);




}

template<typename T>
int SignedDistSolver<T>::sign_function(T pr)
{
	if (pr > 0)
	{
		return 1;
	}
	else if (fabs(pr) < std::numeric_limits<T>::epsilon())
	{
		return 0;
	}
	else
	{
		return -1;

	}
}




//initialize $u0$ where 0 at data points and 100 otherwise
template<typename T>
void SignedDistSolver<T>::initial_dist(T*** u, int*** cp, bool *** lock)
{
	int i, j, k, m;
	T r;
	Point3D<T> tmp;

	for(i=0; i<XPIXEL; i++)
		for(j=0; j<YPIXEL; j++)
			for(k=0; k<ZPIXEL; k++)
				u[i][j][k] = 100;

	for(i=0; i<XPIXEL+2; i++)
		for(j=0; j<YPIXEL+2; j++)
			for(k=0; k<ZPIXEL+2; k++)
				cp[i][j][k] = -1;

	T xmin, ymin, zmin, x, y, z;
	xmin = 0.5*(XPIXEL-1)/PIXEL;
	ymin = 0.5*(YPIXEL-1)/PIXEL;
	zmin = 0.5*(ZPIXEL-1)/PIXEL;
	
	for(m=0; m<vertices.size(); m++)	{
		tmp = vertices[m];
		i = (int)((tmp[0]+xmin)*PIXEL);
		j = (int)((tmp[1]+ymin)*PIXEL);
		k = (int)((tmp[2]+zmin)*PIXEL);

		x=i*1./PIXEL- xmin - tmp[0]; // grid size is 1/PIXEL ????
		y=j*1./PIXEL- ymin - tmp[1];
		z=k*1./PIXEL- zmin - tmp[2];
		r = sqrt( pow(x,2)+pow(y,2)+pow(z,2) );
		if(r<u[i][j][k]){
			u[i][j][k] = r;
			cp[i+1][j+1][k+1] = m; // the number of vertices
		}
		lock[i+1][j+1][k+1] = 1;
		updatelock(lock, i+1, j+1, k+1);

		x=(i+1)*1./PIXEL- xmin - tmp[0];
		r = sqrt( pow(x,2)+pow(y,2)+pow(z,2) );
		if(r<u[i+1][j][k]){
			u[i+1][j][k] = r;
			cp[i+2][j+1][k+1] = m;
		}
		lock[i+2][j+1][k+1] = 1;
		updatelock(lock, i+2, j+1, k+1);

		y=(j+1)*1./PIXEL- ymin - tmp[1];
		r = sqrt( pow(x,2)+pow(y,2)+pow(z,2) );
		if(r<u[i+1][j+1][k]){
			u[i+1][j+1][k] = r;
			cp[i+2][j+2][k+1] = m;
		}
		lock[i+2][j+2][k+1] = 1;
		updatelock(lock, i+2, j+2, k+1);

		z=(k+1)*1./PIXEL- zmin - tmp[2];
		r = sqrt( pow(x,2)+pow(y,2)+pow(z,2) );
		if(r<u[i+1][j+1][k+1]){
			u[i+1][j+1][k+1] = r;
			cp[i+2][j+2][k+2] = m;
		}
		lock[i+2][j+2][k+2] = 1;
		updatelock(lock, i+2, j+2, k+2);

		x=i*1./PIXEL- xmin - tmp[0];
		r = sqrt( pow(x,2)+pow(y,2)+pow(z,2) );
		if(r<u[i][j+1][k+1]){
			u[i][j+1][k+1] = r;
			cp[i+1][j+2][k+2] = m;
		}
		lock[i+1][j+2][k+2] = 1;
		updatelock(lock, i+1, j+2, k+2);

		y=j*1./PIXEL- ymin - tmp[1];
		r = sqrt( pow(x,2)+pow(y,2)+pow(z,2) );
		if(r<u[i][j][k+1]){
			u[i][j][k+1] = r;
			cp[i+1][j+1][k+2] = m;
		}
		lock[i+1][j+1][k+2] = 1;
		updatelock(lock, i+1, j+1, k+2);

		x=(i+1)*1./PIXEL- xmin - tmp[0];
		r = sqrt( pow(x,2)+pow(y,2)+pow(z,2) );
		if(r<u[i+1][j][k+1]){
			u[i+1][j][k+1] = r;
			cp[i+2][j+1][k+2] = m;
		}
		lock[i+2][j+1][k+2] = 1;
		updatelock(lock, i+2, j+1, k+2);

		x=i*1./PIXEL- xmin - tmp[0];
		y=(j+1)*1./PIXEL- ymin - tmp[1];
		z=k*1./PIXEL- zmin - tmp[2];
		r = sqrt( pow(x,2)+pow(y,2)+pow(z,2) );
		if(r<u[i][j+1][k]){
			u[i][j+1][k] = r;
			cp[i+1][j+2][k+1] = m;
		}
		lock[i+1][j+2][k+1] = 1;
		updatelock(lock, i+1, j+2, k+1);
	}
}

//solution of Godunov upwind difference scheme
template<typename T>
void SignedDistSolver<T>::update(T*** u, int*** cp, int i, int j, int k, bool*** lock)
{
	if(cp[i-1][j][k]!=-1){
		if(cp[i][j][k] == -1) {
			cp[i][j][k] = cp[i-1][j][k];
			update3(u, i, j, k, cp[i][j][k]); // update the distance u
			updatelock(lock, i, j, k, u); // extend the lock = 1
		}
		else
			if(cp[i][j][k] != cp[i-1][j][k]){
				if( update4(u, i, j, k, cp[i-1][j][k])==1 ){
					cp[i][j][k] = cp[i-1][j][k];
					updatelock(lock, i, j, k, u);
				}
			}
	}

	if(cp[i+1][j][k]!=-1){
		if(cp[i][j][k] == -1) {
			cp[i][j][k] = cp[i+1][j][k];
			update3(u, i, j, k, cp[i][j][k]);
			updatelock(lock, i, j, k, u);
		}
		else
			if(cp[i][j][k] != cp[i+1][j][k]){
				if( update4(u, i, j, k, cp[i+1][j][k])==1 ){
					cp[i][j][k] = cp[i+1][j][k];
					updatelock(lock, i, j, k, u);
				}
			}
	}
	
	if(cp[i][j-1][k]!=-1){
		if(cp[i][j][k] == -1) {
			cp[i][j][k] = cp[i][j-1][k];
			update3(u, i, j, k, cp[i][j][k]);
			updatelock(lock, i, j, k, u);
		}
		else
			if(cp[i][j][k] != cp[i][j-1][k]){
				if( update4(u, i, j, k, cp[i][j-1][k])==1 ){
					cp[i][j][k] = cp[i][j-1][k];
					updatelock(lock, i, j, k, u);
				}
			}
	}

	if(cp[i][j+1][k]!=-1){
		if(cp[i][j][k] == -1) {
			cp[i][j][k] = cp[i][j+1][k];
			update3(u, i, j, k, cp[i][j][k]);
			updatelock(lock, i, j, k, u);
		}
		else
			if(cp[i][j][k] != cp[i][j+1][k]){
				if( update4(u, i, j, k, cp[i][j+1][k])==1 ){
					cp[i][j][k] = cp[i][j+1][k];
					updatelock(lock, i, j, k, u);
				}
			}
	}

	if(cp[i][j][k-1]!=-1){
		if(cp[i][j][k] == -1) {
			cp[i][j][k] = cp[i][j][k-1];
			update3(u, i, j, k, cp[i][j][k]);
			updatelock(lock, i, j, k, u);
		}
		else
			if(cp[i][j][k] != cp[i][j][k-1]){
				if( update4(u, i, j, k, cp[i][j][k-1])==1 ){
					cp[i][j][k] = cp[i][j][k-1];
					updatelock(lock, i, j, k, u);
				}
			}
	}

	if(cp[i][j][k+1]!=-1){
		if(cp[i][j][k] == -1) {
			cp[i][j][k] = cp[i][j][k+1];
			update3(u, i, j, k, cp[i][j][k]);
			updatelock(lock, i, j, k, u);
		}
		else
			if(cp[i][j][k] != cp[i][j][k+1]){
				if( update4(u, i, j, k, cp[i][j][k+1])==1 ){
					cp[i][j][k] = cp[i][j][k+1];
					updatelock(lock, i, j, k, u);
				}
			}
	}
}

template<typename T>
void SignedDistSolver<T>::updatelock(bool*** lock, int i, int j, int k)
{
	lock[i-1][j][k] = 1;
	lock[i+1][j][k] = 1;
	lock[i][j-1][k] = 1;
	lock[i][j+1][k] = 1;
	lock[i][j][k-1] = 1;
	lock[i][j][k+1] = 1;
}

template<typename T>
void SignedDistSolver<T>::updatelock(bool*** lock, int i, int j, int k, T*** u)
{
	if ( lock[i-1][j][k] == 0 && u[i-1][j][k]>u[i][j][k] ) // extend lock = 1 
		lock[i-1][j][k] = 1;
	if ( lock[i+1][j][k] == 0 && u[i+1][j][k]>u[i][j][k] )
		lock[i+1][j][k] = 1;
	if ( lock[i][j-1][k] == 0 && u[i][j-1][k]>u[i][j][k] )
		lock[i][j-1][k] = 1;
	if ( lock[i][j+1][k] == 0 && u[i][j+1][k]>u[i][j][k] )
		lock[i][j+1][k] = 1;
	if ( lock[i][j][k-1] == 0 && u[i][j][k-1]>u[i][j][k] )
		lock[i][j][k-1] = 1;
	if ( lock[i][j][k+1] == 0 && u[i][j][k+1]>u[i][j][k] )
		lock[i][j][k+1] = 1;
}

//Fast Sweeping Algorithm (using 8 sweeps to guarantee convergence)
template<typename T>
void SignedDistSolver<T>::distance_sweeping(T*** u0, T*** ud, int*** cp, bool*** lock)
{
	int i, j, k;
	T*** u = matrix3D<T>(XPIXEL + 2, YPIXEL + 2, ZPIXEL + 2); // 0 and xpixel+1 is boundary

//set boundary condition
	for(i=0; i<XPIXEL+2; i++){
		for(j=0; j<YPIXEL+2; j++){
			u[i][j][0] = 100;	u[i][j][ZPIXEL+1] = 100;
		}
	}

	for(i=0; i<XPIXEL+2; i++){
		for(k=0; k<ZPIXEL+2; k++){
			u[i][0][k] = 100;	u[i][YPIXEL+1][k] = 100;
		}
	}

	for(j=0; j<YPIXEL+2; j++){
		for(k=0; k<ZPIXEL+2; k++){
			u[0][j][k] = 100;	u[XPIXEL+1][j][k] = 100;
		}
	}

	for(i=0; i<XPIXEL; i++)
		for(j=0; j<YPIXEL; j++)
			for(k=0; k<ZPIXEL; k++)
				u[i+1][j+1][k+1] = u0[i][j][k];

//(1) i=1:n, j=1:n
	for(i=1; i<XPIXEL+1; i++){
		for(j=1; j<YPIXEL+1; j++){
			for(k=1; k<ZPIXEL+1; k++){
				if( lock[i][j][k] ){
					update(u, cp, i, j, k, lock);		
					lock[i][j][k] = 0;
				}
			}
		}
	}

//(2) i=1:n, j=n:1
	for(i=1; i<XPIXEL+1; i++){
		for(j=1; j<YPIXEL+1; j++){
			for(k=ZPIXEL; k>0; k--){
				if( lock[i][j][k] ){
					update(u, cp, i, j, k, lock);		
					lock[i][j][k] = 0;
				}
			}
		}
	}

//(3) i=n:1, j=1:n
	for(i=1; i<XPIXEL+1; i++){
		for(j=YPIXEL; j>0; j--){
			for(k=1; k<ZPIXEL+1; k++){
				if( lock[i][j][k] ){
					update(u, cp, i, j, k, lock);		
					lock[i][j][k] = 0;
				}
			}
		}
	}

//(4) i=n:1, j=n:1
	for(i=1; i<XPIXEL+1; i++){
		for(j=YPIXEL; j>0; j--){
			for(k=ZPIXEL; k>0; k--){
				if( lock[i][j][k] ){
					update(u, cp, i, j, k, lock);		
					lock[i][j][k] = 0;
				}
			}
		}
	}

//(5)
	for(i=XPIXEL; i>0; i--){
		for(j=1; j<YPIXEL+1; j++){
			for(k=1; k<ZPIXEL+1; k++){
				if( lock[i][j][k] ){
					update(u, cp, i, j, k, lock);		
					lock[i][j][k] = 0;
				}
			}
		}
	}

//(6) i=1:n, j=1:n
	for(i=XPIXEL; i>0; i--){
		for(j=1; j<YPIXEL+1; j++){
			for(k=ZPIXEL; k>0; k--){
				if( lock[i][j][k] ){
					update(u, cp, i, j, k, lock);		
					lock[i][j][k] = 0;
				}
			}
		}
	}

//(7) 
	for(i=XPIXEL; i>0; i--){
		for(j=YPIXEL; j>0; j--){
			for(k=1; k<ZPIXEL+1; k++){
				if( lock[i][j][k] ){
					update(u, cp, i, j, k, lock);		
					lock[i][j][k] = 0;
				}
			}
		}
	}

//(8)
	for(i=XPIXEL; i>0; i--){
		for(j=YPIXEL; j>0; j--){
			for(k=ZPIXEL; k>0; k--){
				if( lock[i][j][k] ){
					update(u, cp, i, j, k, lock);		
					lock[i][j][k] = 0;
				}
			}
		}
	}

	for(i=0; i<XPIXEL; i++)	 
		for(j=0; j<YPIXEL; j++)
			for(k=0; k<ZPIXEL; k++){
				ud[i+2][j+2][k+2] = u[i+1][j+1][k+1];
				u0[i][j][k] = u[i+1][j+1][k+1];
			}
	
	delete_matrix3D<T>(u);
}

template<typename T>
void SignedDistSolver<T>::update3(T*** u, int i, int j, int k, int m)
{
	Point3D<T> tmp = vertices[m];

	T xmin, ymin, zmin, x, y, z;

	xmin = 0.5*(XPIXEL-1)/PIXEL;
	ymin = 0.5*(YPIXEL-1)/PIXEL;
	zmin = 0.5*(ZPIXEL-1)/PIXEL;
	
	x=(i-1)*1./PIXEL- xmin - tmp[0];
	y=(j-1)*1./PIXEL- ymin - tmp[1];
	z=(k-1)*1./PIXEL- zmin - tmp[2];
	
	u[i][j][k] = sqrt( pow(x, 2) + pow(y, 2) + pow(z, 2) );
}


template<typename T>
bool SignedDistSolver<T>::update4(T*** u, int i, int j, int k, int m)
{
	Point3D<T> tmp = vertices[m];

	T xmin, ymin, zmin, x, y, z, r;

	xmin = 0.5*(XPIXEL-1)/PIXEL;
	ymin = 0.5*(YPIXEL-1)/PIXEL;
	zmin = 0.5*(ZPIXEL-1)/PIXEL;
	
	x=(i-1)*1./PIXEL- xmin - tmp[0];
	y=(j-1)*1./PIXEL- ymin - tmp[1];
	z=(k-1)*1./PIXEL- zmin - tmp[2];
	
	r = sqrt( pow(x, 2) + pow(y, 2) + pow(z, 2) );

	if(r<u[i][j][k]){
		u[i][j][k] = r;
		return 1;
	}
	else
		return 0;
}

template<typename T>
void SignedDistSolver<T>::unify(T*** u, int n)
{
	int i, j, k;
	T m = 0;
	for (i = n; i<XPIXEL + n; i++)
		for (j = n; j<YPIXEL + n; j++)
			for (k = n; k<ZPIXEL + n; k++)
				m = std::max(m, (T)fabs(u[i][j][k]));
	for (i = n; i<XPIXEL + n; i++)
		for (j = n; j<YPIXEL + n; j++)
			for (k = n; k<ZPIXEL + n; k++)
				u[i][j][k] /= m;
}


template SignedDistSolver<float>::SignedDistSolver();
template SignedDistSolver<float>::~SignedDistSolver();
template void SignedDistSolver<float>::initial_dist(float***, int***, bool ***);
template void SignedDistSolver<float>::distance_sweeping(float***, float***, int***, bool***);
template void SignedDistSolver<float>::calc_signed_dist(float***, int***, float***);
template SignedDistSolver<float>::SignedDistSolver(gridsize, std::vector<Point3D<float>>, std::vector<Point3D<float>>, int);
template void SignedDistSolver<float>::unify(float***, int);

template SignedDistSolver<double>::SignedDistSolver();
template SignedDistSolver<double>::~SignedDistSolver();
template void SignedDistSolver<double>::initial_dist(double***, int***, bool ***);
template void SignedDistSolver<double>::distance_sweeping(double***, double***, int***, bool***);
template void SignedDistSolver<double>::calc_signed_dist(double***, int***, double***);
template SignedDistSolver<double>::SignedDistSolver(gridsize, std::vector<Point3D<double>>, std::vector<Point3D<double>>, int);
template void SignedDistSolver<double>::unify(double***, int);
