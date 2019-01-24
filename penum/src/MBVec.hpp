/*
 * MBVec.hpp
 *
 *  Created on: 29.05.2018
 *      Author: Michael
 */

#ifndef SRC_MBVEC_HPP_
#define SRC_MBVEC_HPP_

#include <vector>
#include <ostream>
#include <boost/align/aligned_allocator.hpp>

using namespace std;

namespace MB {

template <typename T> class Vec {
public:
	Vec() {

	}

	Vec(int si) {
		this->_data.resize(si);
		this->_data.reserve(si);
		for(auto it=_data.begin();it!=_data.end();++it) {
			*it = 0;
		}
	}

	Vec(const Vec& v_in) {
		size_t s = v_in.length();
		_data.reserve(s);
		_data.resize(s);

		for(size_t i=0; i<s; i++) {
			_data[i] = v_in._data[i];
		}
	}

	Vec& operator=( const Vec& other) {
		size_t s = other.length();
		_data.reserve(s);
		_data.resize(s);

		for(int i=0; i<s; i++) {
			_data[i] = other._data[i];
		}

		return *this;
	}



	~Vec() {
		_data.clear();
	}

	void SetLength(int si) {
		_data.clear();
		this->_data.resize(si);
		this->_data.reserve(si);
		for(auto it=_data.begin();it!=_data.end();++it) {
			*it = 0;
		}
	}

	size_t length() const {
		return _data.size();
	}

	const T* data() const {
		return _data.data();
	}

	T* dat() {
		return _data.data();
	}

	// Operators
	template <typename TT> friend std::ostream& operator<<( std::ostream&, const Vec<TT>& );
	const T operator [](int i) const    {return _data[i];}
    T & operator [](int i) {return _data[i];}

private:
	//vector<T> _data;
	std::vector<T, boost::alignment::aligned_allocator<T, 64>> _data;
};

template <typename T>
std::ostream& operator<<( std::ostream& o, const Vec<T>& v) {
	for(auto it=v._data.begin(); it!=v._data.end(); ++it) {
		o << *it << " ";
	}
	o << endl;
   return o;
}

}

#endif /* SRC_MBVEC_HPP_ */
