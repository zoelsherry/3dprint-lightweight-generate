#include "SurfaceProcess.h"
#include "Matrix.h"
#include <iostream>

using std::cout;
using std::endl;

extern int PIXEL;

template<typename T>
void SmoothSurface(gridsize grid, T*** phi, int maxiter)
{
	
	// TODO: 1 should adjust
	float mp = 1;
	for (int i = 0; i < grid.nx; i++)
	{
		for (int j = 0; j < grid.ny; j++)
		{
			for (int k = 0; k < grid.nz; k++)
			{
				phi[i][j][k] = -tanh(phi[i][j][k] / mp * atanh(0.95));
			}
		}
	}

	// TODO: 2 should adjust
	float mp_distance = 8;
	T*** gf1 = matrix3D<T>(grid.nx, grid.ny, grid.nz);
	for (int i = 0; i < grid.nx; i++)
	{
		for (int j = 0; j < grid.ny; j++)
		{
			for (int k = 0; k < grid.nz; k++)
			{
				gf1[i][j][k] = pow(tanh(phi[i][j][k] / mp_distance * atanh(0.95)), 2);
			}
		}
	}

	//driechlet boudary condition
	for (int i = 0; i < grid.nx; i++)
		for (int j = 0; j < grid.ny; j++)
		{
			phi[i][j][grid.nz - 1] = -1;
			phi[i][j][0] = -1;
		}
	for (int i = 0; i < grid.nx; i++)
		for (int k = 0; k < grid.nz; k++)
		{
			phi[i][grid.ny - 1][k] = -1;
			phi[i][0][k] = -1;
		}
	for (int j = 0; j < grid.ny; j++)
		for (int k = 0; k < grid.nz; k++)
		{
			phi[grid.nx - 1][j][k] = -1;
			phi[0][j][k] = -1;
		}

	T epsilon = pow(5.0, 2);
	T dt = 0.15;
	int it = 1;
	for (int i = 0; i < maxiter; i++)
	{
		cout << it << endl;
		for (int i = 1; i < grid.nx - 1; i++)
		{
			for (int j = 1; j < grid.ny - 1; j++)
			{
				for (int k = 1; k < grid.nz - 1; k++)
				{
					phi[i][j][k] = phi[i][j][k] + 0.5*dt*(
						(gf1[i + 1][j][k] + gf1[i][j][k]) * (phi[i + 1][j][k] - phi[i][j][k]) + (gf1[i - 1][j][k] + gf1[i][j][k]) * (phi[i - 1][j][k] - phi[i][j][k])
						+ (gf1[i][j + 1][k] + gf1[i][j][k]) * (phi[i][j + 1][k] - phi[i][j][k]) + (gf1[i][j - 1][k] + gf1[i][j][k]) * (phi[i][j - 1][k] - phi[i][j][k])
						+ (gf1[i][j][k + 1] + gf1[i][j][k]) * (phi[i][j][k + 1] - phi[i][j][k]) + (gf1[i][j][k - 1] + gf1[i][j][k]) * (phi[i][j][k - 1] - phi[i][j][k])
						);

					phi[i][j][k] = phi[i][j][k] / sqrt(exp(-2.0 * dt / epsilon * gf1[i][j][k]) + pow(phi[i][j][k], 2) * (1.0 - exp(-2.0 * dt / epsilon * gf1[i][j][k])));

				}
			}
		}

		it++;

	}

}

template void SmoothSurface<float>(gridsize, float***, int);
template void SmoothSurface<double>(gridsize, double***, int);