#include "StructureFiller.h"
#include "splinterp.h"
#include "Matrix.h"
#include "RangeImpl.h"

template<typename T>
StructureFiller<T>::StructureFiller()
{
	
}

template<typename T>
StructureFiller<T>::StructureFiller(std::vector<T> _nphi, int _interp_num)
{
	interp_num = _interp_num;
	nun = matrix3D<T>(interp_num, interp_num, interp_num);
	cell_interp(_nphi);


}

template<typename T>
StructureFiller<T>::~StructureFiller()
{
	delete_matrix3D<T>(nun);
}

template<typename T>
void StructureFiller<T>::generate_whole(gridsize grid, T *** phi)
{
	int nxp = grid.nx / interp_num - 2;
	int nyp = grid.ny / interp_num - 2;
	int nzp = grid.nz / interp_num - 2;

	T*** phi_t = matrix3D<T>(grid.nx, grid.ny, grid.nz);

	for (int i = 0; i < grid.nx; i++)
		for (int j = 0; j < grid.ny; j++)
			for (int k = 0; k < grid.nz; k++)
			{

				phi_t[i][j][k] = phi[i][j][k];

			}

	for (int ix = 0; ix < nxp; ix++)
		for (int jy = 0; jy < nyp; jy++)
			for (int kz = 0; kz < nzp; kz++)
			{
				int num = 0;
				int num_tar = interp_num * interp_num * interp_num * 9;
				for (int i = 0; i < interp_num; i++)
					for (int j = 0; j < interp_num; j++)
						for (int k = 0; k < interp_num; k++)
						{
							if (phi_t[ix * interp_num + i + 1][jy * interp_num + j + 1][kz * interp_num + k + 1] > 0)
								num++;

							if (phi_t[ix * interp_num + i + 2][jy * interp_num + j][kz * interp_num + k + 2] > 0)
								num++;
							if (phi_t[ix * interp_num + i + 2][jy * interp_num + j + 2][kz * interp_num + k + 2] > 0)
								num++;
							if (phi_t[ix * interp_num + i][jy * interp_num + j][kz * interp_num + k + 2] > 0)
								num++;
							if (phi_t[ix * interp_num + i][jy * interp_num + j + 2][kz * interp_num + k + 2] > 0)
								num++;

							if (phi_t[ix * interp_num + i + 2][jy * interp_num + j][kz * interp_num + k] > 0)
								num++;
							if (phi_t[ix * interp_num + i + 2][jy * interp_num + j + 2][kz * interp_num + k] > 0)
								num++;
							if (phi_t[ix * interp_num + i][jy * interp_num + j][kz * interp_num + k] > 0)
								num++;
							if (phi_t[ix * interp_num + i][jy * interp_num + j + 2][kz * interp_num + k] > 0)
								num++;

						}

				if (num_tar == num)
				{
					for (int i = 0; i < interp_num; i++)
						for (int j = 0; j < interp_num; j++)
							for (int k = 0; k < interp_num; k++)
								phi[ix * interp_num + i + 1][jy * interp_num + j + 1][kz * interp_num + k + 1] = nun[i][j][k];
				}
			}

	delete_matrix3D<T>(phi_t);

}

template<typename T>
void StructureFiller<T>::generate_inner(gridsize grid, T *** phi)
{
	int nxp = grid.nx / interp_num - 2;
	int nyp = grid.ny / interp_num - 2;
	int nzp = grid.nz / interp_num - 2;

	T*** phi_t = matrix3D<T>(grid.nx, grid.ny, grid.nz);

	for (int i = 0; i < grid.nx; i++)
		for (int j = 0; j < grid.ny; j++)
			for (int k = 0; k < grid.nz; k++)
			{

				phi_t[i][j][k] = phi[i][j][k];

			}


	for (int ix = 0; ix < nxp; ix++)
		for (int jy = 0; jy < nyp; jy++)
			for (int kz = 0; kz < nzp; kz++)
			{
				int num = 0;
				for (int i = 0; i < interp_num; i++)
					for (int j = 0; j < interp_num; j++)
						for (int k = 0; k < interp_num; k++)
						{
							if (phi_t[ix * interp_num + i][jy * interp_num + j][kz * interp_num + k] > 0)
								num++;
						}
				if (num > 1)
				{
					for (int i = 0; i < interp_num; i++)
						for (int j = 0; j < interp_num; j++)
							for (int k = 0; k < interp_num; k++)
								phi[ix * interp_num + i][jy * interp_num + j][kz * interp_num + k] =
								-(1 - phi_t[ix * interp_num + i][jy * interp_num + j][kz * interp_num + k]) + phi_t[ix * interp_num + i][jy * interp_num + j][kz * interp_num + k] * nun[i][j][k];
				}
			}

	delete_matrix3D<T>(phi_t);

}


template<typename T>
void StructureFiller<T>::cell_interp(std::vector<T> _nphi)
{
	double ms_ori = exp(log(_nphi.size()) / 3); // can not use int because lose accuracy
	//
	//cout << ceil(ms_ori) << endl;
	//
	ms_ori = ceil(ms_ori); // transform the accurate int

	int len = ms_ori * ms_ori * ms_ori;
	T *data = new T[len];
	for (size_t i = 0; i < _nphi.size(); i++)
	{
		data[i] = _nphi[i];
	}

	int len_tar = interp_num * interp_num * interp_num;
	T *linex = new T[interp_num];
	T *liney = new T[interp_num];
	T *linez = new T[interp_num];
	linex = crange::template linespace<T>(0, ms_ori - 1, interp_num);
	liney = crange::template linespace<T>(0, ms_ori - 1, interp_num);
	linez = crange::template linespace<T>(0, ms_ori - 1, interp_num);
	T** xxyyzz = matrix2D<T>(3, len_tar);
	xxyyzz = crange::template meshgrid<T>(linex, liney, linez, interp_num, interp_num, interp_num);

	T *nun_t = new T[len_tar];
	splinterp::parallel_interp3(splinterp::interp3_F<T>, data, ms_ori, ms_ori, ms_ori, xxyyzz[0], xxyyzz[1], xxyyzz[2], len_tar, nun_t, 0);

	int nun_c = 0;
	for (int i = 0; i < interp_num; i++)
	{
		for (int j = 0; j < interp_num; j++)
		{
			for (int k = 0; k < interp_num; k++)
			{
				nun[i][j][k] = nun_t[nun_c++];
			}
		}
	}

	delete[] data;
	delete[] linex;
	delete[] liney;
	delete[] linez;

	delete_matrix2D<T>(xxyyzz);
	delete[] nun_t;

#ifdef DEBUG
	std::ofstream out_structure_unit;
	out_structure_unit.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try
	{
		out_structure_unit.open("./outdata/structure_unit.txt");
		for (int i = 0; i < ms_tar; i++)
			for (int j = 0; j < ms_tar; j++)
				for (int k = 0; k < ms_tar; k++)
					out_structure_unit << nun[i][j][k] << std::endl;

		out_structure_unit.close();
	}
	catch (std::ofstream::failure e)
	{
		std::cout << "ERROR::SIGNEDDIST::FILE FAIL TO WRITE!" << std::endl;
	}
#endif // DEBUG



}

template StructureFiller<float>::StructureFiller();
template StructureFiller<float>::~StructureFiller();
template StructureFiller<float>::StructureFiller(std::vector<float> , int );
template void StructureFiller<float>::generate_whole(gridsize, float***);
template void StructureFiller<float>::generate_inner(gridsize, float***);


template StructureFiller<double>::StructureFiller();
template StructureFiller<double>::~StructureFiller();
template StructureFiller<double>::StructureFiller(std::vector<double>, int);
template void StructureFiller<double>::generate_whole(gridsize, double***);
template void StructureFiller<double>::generate_inner(gridsize, double***);
