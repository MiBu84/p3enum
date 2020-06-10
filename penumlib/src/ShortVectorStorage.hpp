/*
 * ShortVectorStorage.hpp
 *
 *  Created on: 05.06.2020
 *      Author: mburger
 */

#ifndef SRC_SHORTVECTORSTORAGE_HPP_
#define SRC_SHORTVECTORSTORAGE_HPP_

#include <map>
#include <vector>

using namespace std;

class ShortVectorStorage {
public:
	static ShortVectorStorage& getInstance()
	{
		static ShortVectorStorage instance;
		return instance;
	}

	void reset();
	int insertVector(double* invec, const double inlen);
	void printLengthInMap();
	void printMap();
	void setDimension(int dimension) {_dim = dimension;}

private:
	static ShortVectorStorage* _instance;
	map<double, vector<int>> _vecs;
	int _dim;

	short vectorEquality(vector<int>& v1, vector<int>& v2);
	ShortVectorStorage() {
		this->_dim = -1;
	}
};


#endif /* SRC_SHORTVECTORSTORAGE_HPP_ */
