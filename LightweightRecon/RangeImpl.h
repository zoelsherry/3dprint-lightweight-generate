#ifndef RANGEIMPL
#define RANGEIMPL

#include <iostream>
#include "Matrix.h"

namespace crange
{

	template<typename value_t>
	value_t* linespace(value_t begin, value_t end, int inter_num)
	{
		value_t* line = new value_t[inter_num];
		for (int i = 0; i < inter_num; i++)
		{
			line[i] = begin + i * (end - begin) / (inter_num - 1);
		}
		
		return line;
	}

	template<typename value_t>
	value_t** meshgrid(value_t* xx, value_t* yy, value_t* zz, int length_x, int length_y, int length_z)
	{
		value_t** mgrid = matrix2D<value_t>(3, length_x * length_y * length_z);
		for (int i = 0; i < length_z; i++)
		{
			for (int j = 0; j < length_x; j++)
			{
				for (int k = 0; k < length_y; k++)
				{
					mgrid[0][i*length_x*length_y + j*length_y +k] = xx[j];
				}
			}
		}

		for (int i = 0; i < length_z; i++)
		{
			for (int j = 0; j < length_x; j++)
			{
				for (int k = 0; k < length_y; k++)
				{
					mgrid[1][i*length_x*length_y + j * length_y + k] = yy[k];
				}
			}
		}

		for (int i = 0; i < length_z; i++)
		{
			for (int j = 0; j < length_x; j++)
			{
				for (int k = 0; k < length_y; k++)
				{
					mgrid[2][i*length_x*length_y + j * length_y + k] = zz[i];
				}
			}
		}

		return mgrid;
	}

	template<typename value_t>
	class RangeImpl
	{
		class Iterator;
	public:
		RangeImpl(value_t begin, value_t end, value_t step = 1) :m_begin(begin), m_end(end), m_step(step)
		{
			if (step>0 && m_begin >= m_end)

				throw std::logic_error("end must greater than begin.");

			else if (step<0 && m_begin <= m_end)

				throw std::logic_error("end must less than begin.");

			m_step_end = (m_end - m_begin) / m_step + 1  ; // plus one more!

		/*	if (m_begin + m_step_end * m_step != m_end)
			{
				m_step_end++;
			}*/

		}

		Iterator begin()
		{
			return Iterator(0, *this);
		}

		Iterator end()
		{
			return Iterator(m_step_end, *this);    // return next to the last one!!!
		}

		value_t operator[](int s)
		{
			return m_begin + s * m_step;
		}

		int size()
		{

			return m_step_end;
		}

	private:
		value_t m_begin;
		value_t m_end;
		value_t m_step;
		int m_step_end;

		class Iterator
		{
		public:
			Iterator(int start, RangeImpl& range) : m_current_step(start), m_range(range)
			{
				m_current_value = m_range.m_begin + m_current_step * m_range.m_step;
			}

			value_t operator*() { return m_current_value; }

			const Iterator* operator++()
			{
				m_current_value += m_range.m_step;
				m_current_step++;
				return this;
			}

			bool operator==(const Iterator& other)
			{
				return m_current_step == other.m_current_step;
			}

			bool operator!=(const Iterator& other)
			{
				return m_current_step != other.m_current_step;
			}

			const Iterator* operator--()
			{
				m_current_value -= m_range.m_step;
				m_current_step--;
				return this;
			}

		private:
			value_t m_current_value;
			int m_current_step;
			RangeImpl& m_range;

		};

	};

	template<typename T, typename V>
	auto Range(T begin, T end, V stepsize)->RangeImpl<T>
	{
		return RangeImpl<T>(begin, end, stepsize);
	}

	template<typename T>
	RangeImpl<T> Range(T begin, T end)
	{
		return RangeImpl<T>(begin, end, 1);
	}

	template<typename T>
	RangeImpl<T> Range(T end)
	{
		return RangeImpl<T>(T(), end, 1);
	}

}

#endif // !RANGEIMPL

