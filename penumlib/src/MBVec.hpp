/*
 * MBVec.hpp
 *
 *  Created on: 29.05.2018
 *      Author: Michael
 */

#ifndef SRC_MBVEC2_HPP_
#define SRC_MBVEC2_HPP_

#include <vector>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <boost/align/aligned_allocator.hpp>

using namespace std;

namespace MB {

template <typename T> class MBVec {
public:
	MBVec() {
		len = 0;
	}

	MBVec(int si) {
		this->_data.resize(si);
		this->_data.reserve(si);
		for(auto it=_data.begin();it!=_data.end();++it) {
			*it = 0;
		}
		len = 0;
	}

	MBVec(int si, double* vals) {
		this->_data.resize(si);
		this->_data.reserve(si);
		int i=0;
		for(auto it=_data.begin();it!=_data.end();++it) {
			*it = vals[i];
			i++;
		}
		len = 0;
	}

	MBVec(int si, double* vals, int s, int e) {
		this->_data.resize(si);
		this->_data.reserve(si);
		int cnt = 0;
		for(int i=s; i<=e; i++) {
			this->_data[cnt] = vals[i];
			cnt++;
		}
		len = 0;
	}

	MBVec(const MBVec& v_in) {
		size_t s = v_in.length();
		_data.reserve(s);
		_data.resize(s);

		for(size_t i=0; i<s; i++) {
			_data[i] = v_in._data[i];
		}
		len = v_in.len;
	}

	MBVec& operator=( const MBVec& other) {
		size_t s = other.length();
		_data.reserve(s);
		_data.resize(s);

		for(size_t i=0; i<s; i++) {
			_data[i] = other._data[i];
		}
		len = other.len;
		return *this;
	}



	~MBVec() {
		_data.clear();
		len = 0;
	}


	void SetLength(int si) {
		_data.clear();
		this->_data.resize(si);
		this->_data.reserve(si);
		for(auto it=_data.begin();it!=_data.end();++it) {
			*it = 0;
		}
	}

	std::ostream& toStream( std::ostream& o, const MBVec<T>& v, int s, int e ) {
		for(int i=s; i<=e; i++) {
			int datval = (int)_data[i];
			o << setw(4) <<  setw(4) << setfill(' ') << datval << " ";
			std::cout << setw(4) <<  setw(4) << setfill(' ') << datval << " ";
		}
		return o;
	}

	size_t length() const {
		return _data.size();
	}

	const T* data() const {
		return _data.data();
	}

	T* data() {
		return _data.data();
	}

	// Operators
	template <typename TT> friend std::ostream& operator<<( std::ostream&, const MBVec<TT>& );


	template <typename TT>
	bool operator==(const MBVec<TT>& v2) const {
		bool equ = true;
		bool sign_reversed = false;
		bool found_first_non_zero = false;

		if(v2.length() != length())
			return false;

		else {
			for(size_t i=0; i<v2.length(); i++) {
				if(!found_first_non_zero && v2[i] == -this->_data[i] && abs(v2[i]) > 1e-10) {
					sign_reversed = true;
					found_first_non_zero = true;
				}

				if(v2[i] != this->_data[i] && !sign_reversed)
					return false;

				if(v2[i] != -this->_data[i] && sign_reversed)
					return false;

			}
		}

		return equ;
	}

	template <typename TT>
	bool operator <(const MBVec<TT>& v2) const {
		if(length() < v2.length())
			return true;

		if(length() > v2.length())
			return false;

		for(int i=0; i < v2.length(); i++) {
			if(_data[i] < v2[i])
				return true;
		}
		return false;
	}


	const T operator [](int i) const    {return _data[i];}
    T & operator [](int i) {return _data[i];}

private:
	vector<T> _data;
	//std::vector<T, boost::alignment::aligned_allocator<T, 64>> _data;
	double len;
};

template <typename T>
std::ostream& operator<<( std::ostream& o, const MBVec<T>& v) {
	int pos = 0;
	for(auto it=v._data.begin(); it!=v._data.end(); ++it) {

		T datval = (T)*it;
		o << "["<<pos<<"]" << setw(4) <<  setw(4) << setfill(' ') << datval << " ";
		pos++;
	}
	//o << endl;
   return o;
}



}; // End Namespace MB


#endif /* SRC_MBVEC2_HPP_ */
