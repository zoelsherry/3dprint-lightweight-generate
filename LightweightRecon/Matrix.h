#ifndef MATRIX_H
#define MATRIX_H

#include <string>

template<typename T>
T** matrix2D(int row, int col);

template<typename T>
void delete_matrix2D(T **m);

template<typename T>
T*** matrix3D(int row, int col, int dep);

template<typename T>
void delete_matrix3D(T ***m);

// TODO(1): I will use it later!!!

//template<typename T>
//class Matrix2D
//{
//public:
//	Matrix2D() {};
//	Matrix2D(int row, int col);
//	~Matrix2D();
//	Matrix2D(const Matrix2D<T> &ma);
//	virtual void reset();
//
//	int nrow;
//	int ncol;
//	int size;
//private:
//	T** m;
//	T* operator[](int k);
//	void operator=(const Matrix2D<T> &mat);
//
//};
//
//template<typename T>
//class Matrix3D
//{
//public	:
//	Matrix3D() {};
//	Matrix3D(int row, int col, int depth);
//	~Matrix3D();
//	Matrix3D(const Matrix3D<T> &mat);
//	virtual void reset();
//	
//	int nrow;
//	int ncol;
//	int ndepth;
//	int size;
//private:
//	T*** m;
//	T operator()(int i, int j, int k);
//	void operator=(const Matrix3D<T> &mat);
//};


#endif // !MATRIX_H

