/*
 * ShortVectorStorage.cpp
 *
 *  Created on: 05.06.2020
 *      Author: mburger
 */

#include "ShortVectorStorage.hpp"
#include <iostream>
#include <cassert>

void ShortVectorStorage::reset() {
	this->_dim = -1;
	this->_vecs.clear();
}

int ShortVectorStorage::insertVector(double* invec, const double inlen) {
	if(_dim < 1)
		cerr << "ShortVectorStorage::insertVector not initialized with dimension." << endl;

	assert (_dim > 0);

	vector<int> insertvector;

	insertvector.reserve(_dim);
	insertvector.resize(_dim);
	for(int i=0; i<_dim; i++) {
		insertvector[i] = (int)invec[i];
	}

	// Check whether the vector is already contained
	for (auto it=_vecs.begin(); it!=_vecs.end(); ++it) {
		if(vectorEquality(insertvector, it->second))
			return 1;
	}

	// If new, then insert it
	this->_vecs.insert(std::pair<double, vector<int>> (inlen, insertvector));

	return 0;

}

void ShortVectorStorage::printLengthInMap() {
	cout << "Have " << _vecs.size() << " vectors:" << endl;
	for (auto it=_vecs.begin(); it != _vecs.end(); ++it) {
		cout << it->first << " ; ";
	}
	cout << endl;
}

void ShortVectorStorage::printMap() {
	cout << "Have " << _vecs.size() << " vectors:" << endl;
	for (auto it=_vecs.begin(); it != _vecs.end(); ++it) {
		cout << it->first << " : [ ";
		for (auto it_in=it->second.begin(); it_in != it->second.end(); ++it_in) {
			cout << *it_in << " ";
		}
		cout << "]" << endl;

	}
	cout << endl;
}

short ShortVectorStorage::vectorEquality(vector<int>& v1, vector<int>& v2) {
	if(v1.size() != v2.size())
		cerr << "ShortVectorStorage::vectorEquality vec1.size="
		<< v1.size() << " but v2.size=" << v2.size() << endl;
	assert (v1.size() == v2.size());

	auto it1=v1.begin();
	auto it2=v2.begin();
	//int pos = 1;

	for(; it1!=v1.end(); ++it1, ++it2) {
		if (*it1 > *it2 || *it1 < *it2)
			return 0;
		//pos++;
	}

	return -1;
}

