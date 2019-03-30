#ifndef STRUCTUREFILLER_H
#define STRUCTUREFILLER_H

#include <vector>
#include "Geometry.h"

template<typename T>
class StructureFiller
{
public:
	StructureFiller();
	StructureFiller(std::vector<T> _nphi, int _interp_num);
	~StructureFiller();

	void generate_whole(gridsize grid, T*** phi);

	void generate_inner(gridsize grid, T*** phi);

private:
	void cell_interp(std::vector<T> _nphi);


	int interp_num;
	T*** nun;



};

#endif // !STRUCTUREFILLER_H

