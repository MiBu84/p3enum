/*
 * VectorStorage.hpp
 *
 *  Created on: 18.07.2018
 *      Author: mburger
 */

#ifndef SRC_VECTORSTORAGE_HPP_
#define SRC_VECTORSTORAGE_HPP_

#include <cstring>
#include <string>
#include <map>
#include <fstream>
#include "MBVec.hpp"

using namespace MB;
using namespace std;

typedef std::map<MBVec<double>, double> vectable;
typedef std::pair<MBVec<double>, double> vectabentry;

class VectorStorage {
public:
	static VectorStorage& getInstance()
	{
		static VectorStorage instance;
		return instance;
	}

	static VectorStorage* _instance;

	vectable tab1;
	vectable tab2;
	vectable tab3;
	int cand_count;
	unsigned long long node_cnt;
	unsigned long long mult_cnt;

	void printVecTable(const vectable& vectab, string filename) {
		std::ofstream outfile;
		outfile.open(filename, std::ios_base::app);

		for(auto it=vectab.begin(); it!=vectab.end(); it++) {
			outfile << it->first;
			outfile	<< "\t" << it->second << endl;
		}

		outfile.close();
	}

private:
	VectorStorage() {
		cand_count = 0;
		node_cnt = 0;
		mult_cnt = 0;
	}

	VectorStorage(const VectorStorage&);
	~VectorStorage() {}



};


#endif /* SRC_VECTORSTORAGE_HPP_ */
