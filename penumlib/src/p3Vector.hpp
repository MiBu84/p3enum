/*
 * p3Vector.hpp
 *
 *  Created on: 17.06.2020
 *      Author: mburger
 */

#ifndef P3VECTOR_H
#define P3VECTOR_H

#include <cassert>
#include <iostream>
#include <vector>

using namespace std;

namespace p3 {

template<class T>
class p3Vector {
public:
	p3Vector() {
		_data.resize(0);
		_dim=0;
	}

	p3Vector(unsigned int dim) {
		_data.resize(dim);

		for(unsigned int i=0; i<dim; i++) {
			_data[i] = T(0);
		}

		_dim=dim;
	}

	p3Vector(const p3Vector& v) {
		_data = v._data;
		_dim=v._dim;
	}

	p3Vector(const std::vector<T>& stdv) {
		_data = stdv;
		_dim= stdv.size();
	}

	// getters and setters
	unsigned int size() const {
		return _dim;
	}

	unsigned int length() const {
		return size();
	}

	void resize(unsigned int s)  {
		_data.resize(s);
		_dim = s;
	}

	void SetLength(unsigned int s)  {
		_data.resize(s);
		_dim = s;
	}

	// operators
	p3Vector operator=(const p3Vector & v) {
		_data = v._data;
		_dim=v._dim;
		return *this;
	}

	T& operator[] (int index) {
		return _data[index];
	}

	const T& operator[] (int index) const {
		return _data[index];
	}

	T operator* (const p3Vector& x) {
		T sum = T(0);
		const T* vv1 = _data.data();
		const T* vv2 = x._data.data();

		for(int i=0; i<x._dim; i++) {
			sum += vv1[i] * vv2[i];
		}

		return sum;
	}

	T operator* (const vector<T>& x) {
		T sum = T(0);

		for(unsigned int i=0; i<x.size(); i++) {
			sum += _data[i] * x[i];
		}

		return sum;
	}

	p3Vector operator* (const T& x) {
		p3Vector resvec = p3Vector(_dim);

		for(int i=0; i<_dim; i++) {
			resvec[i] = x * _data[i];
		}

		return resvec;
	}

	friend p3Vector operator* (const T alpha, const p3Vector<T>& v) {
		p3Vector resvec = p3Vector(v._dim);

		for(int i=0; i<v._dim; i++) {
			resvec[i] = alpha * v[i];
		}

		return resvec;
	}

	friend p3Vector operator* (const T alpha, const vector<T>& v) {
		p3Vector resvec = p3Vector(v.size());

		for(int i=0; i<v.size(); i++) {
			resvec[i] = alpha * v[i];
		}

		return resvec;
	}

	friend p3Vector operator* (const p3Vector<T>& v1, const p3Vector<T>& v2) {
		T sum = T(0);

		for(unsigned int i=0; i<v1.size(); i++) {
			sum += v1[i] * v2[i];
		}

		return sum;
	}

	p3Vector operator- (const p3Vector& x) {
		p3Vector resvec = p3Vector(x._dim);

		for(int i=0; i<x._dim; i++) {
			resvec[i] = _data[i] - x[i];
		}
		return resvec;
	}

	p3Vector operator+ (const p3Vector& x) {
		p3Vector resvec = p3Vector(x._dim);

		for(int i=0; i<x._dim; i++) {
			resvec[i] = _data[i] + x[i];
		}

		return resvec;
	}

	friend std::ostream& operator<<(ostream& os, const p3Vector<T>& dt) {
		os << "[ ";

		for (unsigned int i=0; i<dt.size(); i++) {
			os << dt[i] << " ";
		}

		os << "]" << endl;

		return os;
	}

private:
	std::vector<T> _data;
	int _dim;
};


} // namespace p3

#endif


