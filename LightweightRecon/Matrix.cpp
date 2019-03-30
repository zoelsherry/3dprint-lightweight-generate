#include "Matrix.h"


template<typename T>
T** matrix2D(int row, int col)
{
	T** m = new T*[(size_t)row];
	m[0] = new T[(size_t)(row * col)]();
	
	for (int i = 1; i < row; i++)
	{
		m[i] = m[i - 1] + col;
	}
	memset(m[0], 0, row * col * sizeof(T));
	return m;
}

template<typename T>
void delete_matrix2D(T** m)
{
	if (m != nullptr)
	{
		delete[] m[0];
		delete[] m;
		m = nullptr;
	}
}

template<typename T>
T*** matrix3D(int row, int col, int dep)
{
	T*** m = new T**[(size_t)row];
	m[0] = new T*[(size_t)(row * col)];
	m[0][0] = new T[(size_t)(row * col * dep)]();
	
	for (int i = 1; i < row; i++)
	{
		m[i] = m[i - 1] + col;
		m[i][0] = m[i - 1][0] + col * dep;
		for (int j = 1; j < col; j++)
		{
			m[i][j] = m[i][j - 1] + dep;
		}
	}
	for (int i = 1; i < col; i++)
	{
		m[0][i] = m[0][i - 1] + dep;
	}
	memset(m[0][0], 0, row * col * dep * sizeof(T));
	return m;
}

template<typename T>
void delete_matrix3D(T *** m)
{
	if (m != nullptr)
	{
		delete[] m[0][0];
		delete[] m[0];
		delete[] m;
		m = nullptr;
	}
}


template int** matrix2D<int>(int, int);
template void delete_matrix2D<int>(int**);
template float** matrix2D<float>(int, int);
template void delete_matrix2D<float>(float**);
template double** matrix2D<double>(int, int);
template void delete_matrix2D<double>(double**);
template bool** matrix2D<bool>(int, int);
template void delete_matrix2D<bool>(bool**);

template int*** matrix3D<int>(int, int, int);
template void delete_matrix3D<int>(int***);
template float*** matrix3D<float>(int, int, int);
template void delete_matrix3D<float>(float***);
template double*** matrix3D<double>(int, int, int);
template void delete_matrix3D<double>(double***);
template bool*** matrix3D<bool>(int, int, int);
template void delete_matrix3D<bool>(bool***);

