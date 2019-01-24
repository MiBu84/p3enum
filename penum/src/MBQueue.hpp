/*
 * MBQueue.hpp
 *
 *  Created on: 04.06.2018
 *      Author: MichaelB
 */

#ifndef SRC_MBQUEUE_HPP_
#define SRC_MBQUEUE_HPP_

#include <queue>
#include "MBVec.hpp"

#ifdef USE_MB_TYPES
#define VECTYPE MB::Vec
#define MATTYPE MB::Mat
#include "MBMat.hpp"
using namespace MB;
#else
#define VECTYPE NTL::Vec
#define MATTYPE NTL::Mat
using namespace NTL;
#endif

using namespace std;

namespace MB {
template <typename T> class MBVecQueue {
public:
	MBVecQueue(int size) {
		_data.resize(size);
		_data.reserve(size);
		_write_pointer=0;
		_read_pointer=0;
		_size=size;
		_elems=0;
		_thresval=0.25 * size;

		//cout << "Creating MBRingQueue with " << _size << " space." << endl;
	}

	/**
	 * Returns the number of elements after insertion
	 * If vector could not be inserted then -1 indicates that queue is full
	 */
	int push(const VECTYPE<T>& v) {
		{
			if(_elems < _size) {
				_data[_write_pointer] = v;
				_elems++;
				// Go round if we reached the end of the vector
				_write_pointer=(_write_pointer+1) % _size;
				return _elems;

			}

			return -1;
		}
	}

	VECTYPE<T>& next () {
		{
			if(_elems > 0 && _read_pointer >= 0) {
				int old_pointer = _read_pointer;
				_elems--;
				// Go round if we reached the end of the vector
				_read_pointer=(_read_pointer+1) % _size;
				//cout << printf("Requesting entry %d.\n", _elems+1);
				return _data[old_pointer];
			}
			else {
				cerr << printf("Requesting illegal entry.\n");
				return _data[0];
			}
		}
	}

	/**
	 * Indicate whether queue is full
	 * To be sure there is a slide security buffer between size and number of elements
	 */
	bool isFull() const {
		if(_elems < _size)
			return false;
		else
			return true;
	}

	bool isEmpty() const {
		if(_elems > 0)
			return false;
		else
			return true;
	}

	int elemsSize() const {
		return _elems;
	}

	bool isBelowThres() const {
		if(_elems > _thresval)
			return false;
		else
			return true;
	}

	bool checkDuplicates() {
		bool dubs = false;
		int ridx_out = _read_pointer;
		// Run over all registered vectors
		for(int i =_elems; i>0; i--) {

			// Check against other vectors
			int ridx_in = _read_pointer;
			for(int j =_elems; j>0; j--) {
				// Don't compare with itself
				if(ridx_out == ridx_in) {
					; // Do nothing
				}

				else {
					if(_data[ridx_out] ==_data[ridx_in]  ) {
						cout << "Fire" << endl;
					}
				}
				ridx_in=(ridx_in+1) % _size;
			}

			ridx_out=(ridx_out+1) % _size;
		}
		return dubs;
	}

private:
	std::vector<VECTYPE<T> > _data;
	int _write_pointer;
	int _read_pointer;
	int _size;
	int _elems;
	int _thresval;
};
}

#endif /* SRC_MBQUEUE_HPP_ */
