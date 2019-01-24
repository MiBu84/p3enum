/*
 * MBMat.hpp
 *
 *  Created on: 29.05.2018
 *      Author: Michael
 */

#ifndef SRC_MBMAT_HPP_
#define SRC_MBMAT_HPP_

#include "MBVec.hpp"
#include <vector>
#include <iostream>
#include <boost/align/aligned_allocator.hpp>

using namespace std;

namespace MB {
	template <typename T> class Mat {
	public:
		Mat() {
			_x=0;
			_y=0;
		}

		Mat(int x, int y) {
			_data.reserve(y);
			_data.resize(y);
			for(int i=0; i<y; i++) {
				Vec<T> tvec(x);
				_data[i] = tvec;
			}
			_x=x;
			_y=y;
		}

		void SetDims(int x, int y) {
			_data.clear();
			for(int i=0; i<y; i++) {
				Vec<T>row_tvec(x);
				_data.push_back(row_tvec);
			}
			_x=x;
			_y=y;
		}

		long NumCols() const {
			return _x;
		}

		long NumRows() const {
			return _y;
		}

		// Operators
		template <typename TT> friend std::ostream& operator<<( std::ostream&, const Mat<TT>& );
	    const Vec<T>& operator [](int y) const    {return _data[y];}
	    Vec<T> & operator [](int y) {return _data[y];}

		const T* data(int y) const {
			return _data[y].data();
		}

	private:
		//std::vector<MB::Vec<T>, boost::alignment::aligned_allocator<T, 64>> _data;
		std::vector<MB::Vec<T>> _data;
		size_t _x;
		size_t _y;
	};

	template <typename T>
	std::ostream& operator<<( std::ostream& o, const Mat<T>& v) {
		o << "[";
		for(size_t i=0; i<v._y; i++) {
			o << "[ ";
			for(size_t j=0; j < v._x; j++) {
				o << v._data[i][j] << " ";
			}
			o << "] ";
			if(i<v._y-1)
				o << endl;
		}
		o << "]" << endl;
		return o;
	}

}
#endif /* SRC_MBMAT_HPP_ */
