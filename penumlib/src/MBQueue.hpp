/*
 * MBQueue.hpp
 *
 *  Created on: 04.06.2018
 *      Author: MichaelB
 */

#ifndef SRC_MBQUEUE2_HPP_
#define SRC_MBQUEUE2_HPP_

#include <queue>
#include <iostream>
#include "MBVec.hpp"

using namespace std;

namespace MB {

class MBVecQueue3 {
public:
	MBVecQueue3(int size) {
		_write_pointer=0;
		_read_pointer=0;
		_size=size;
		_elems=0;
		_thresval=0.25 * size;

		_data.resize(size+2);
		_data.reserve(size+2);

		//cout << "Creating MBRingQueue with " << _size << " space." << endl;
	}


	MBVecQueue3() {
		_write_pointer=0;
		_read_pointer=0;
		_size=0;
		_elems=0;
		_thresval=0.25 * 0;
	}

	int push(const MBVec<double>& v, const int start=-1, const int end=-1) {
		{
			//cout << "Pushing: " << _elems << endl;
			_data[_write_pointer] = v;
			if(_elems < _size) {
				_elems++;

				// Go round if we reached the end of the vector
				_write_pointer=(_write_pointer+1) % _size;
				return _elems;

			}

			return -1;
		}
	}


	bool next (MBVec<double>& input, const int start, const int end) {
		{
			if(_elems > 0 && _read_pointer >= 0) {
				int old_pointer = _read_pointer;
				_elems--;
				// Go round if we reached the end of the vector
				_read_pointer=(_read_pointer+1) % _size;

				for(int i=start; i<=end; i++) {
					input[i] = _data[old_pointer][i];
				}
				return true;
			}
			else {
				cerr << "Requesting illegal entry." << endl;
				return false;
			}
		}
	}

	void clear() {
		_data.clear();
		_write_pointer=0;
		_read_pointer=0;
		_elems=0;
		_thresval=0.25 * 0;
	}


	void print() {
		cout << "Printing queue:" << endl;
		int i = _read_pointer;
		while(i != _write_pointer) {

			if(i >= 0) {
				cout << _data[i]  << endl;;
			}
			i=(i+1) % _size;
		}
		cout << "End of queue:" << endl;
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

private:
	int _write_pointer;
	int _read_pointer;
	int _size;
	int _elems;
	int _thresval;

	std::vector<MBVec<double> > _data;
};



template <class T> class MBVecQueue2 {
public:
	MBVecQueue2(int size) {

//#ifdef USE_VECQUEUE
	_data.resize(size);
	_data.reserve(size);
/*#else
	_data = new double*[size];
	for(int i=0; i < size; i++) {
		_data[i] = new double[152]; // Only vectors with about 150 entries at the moment
		for(int j=0; j < 152; j++) {
			_data[i][j] = 0.0;
		}
	}
#endif*/
		_write_pointer=0;
		_read_pointer=0;
		_size=size;
		_elems=0;
		_thresval=0.25 * size;

		cout << "Creating MBRingQueue with " << _size << " space." << endl;
	}

	MBVecQueue2() {
//#ifdef USE_VECQUEUE
		// do nothing
/*#else
		_data = NULL;
#endif*/
		_write_pointer=0;
		_read_pointer=0;
		_size=0;
		_elems=0;
		_thresval=0.25 * 0;
	}

	/**
	 * Returns the number of elements after insertion
	 * If vector could not be inserted then -1 indicates that queue is full
	 */
	int push(const MBVec<T>& v, const int start=-1, const int end=-1) {
		{
			if(_elems < _size) {
//#ifdef USE_VECQUEUE
				_data[_write_pointer] = v;
/*#else
				if(start==-1 || end==-1) {
					cerr << "Error: Accessing ArrayQueue without bounds!" << endl;
					for(int i=0; i<=152; i++) {
						_data[_write_pointer][i] = v[i];
					}
				}
				else {
					for(int i=start; i<=end; i++) {
						_data[_write_pointer][i] = v[i];
					}
				}

#endif*/
				_elems++;

				// Go round if we reached the end of the vector
				_write_pointer=(_write_pointer+1) % _size;
				return _elems;

			}

			return -1;
		}
	}

//#ifdef USE_VECQUEUE
	/*MBVec<T>& next () {
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
				cerr << "Requesting illegal entry." << endl;
				return _data[0];
			}
		}
	}*/
//#endif

	/*bool next (MBVec<T>& input, const int start, const int end) {
		{
			if(_elems > 0 && _read_pointer >= 0) {
				int old_pointer = _read_pointer;
				_elems--;
				// Go round if we reached the end of the vector
				_read_pointer=(_read_pointer+1) % _size;

				for(int i=start; i<=end; i++) {
					input[i] = _data[old_pointer][i];
				}
				return true;
			}
			else {
				cerr << "Requesting illegal entry." << endl;
				return false;
			}
		}
	}*/

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

	/*bool checkDuplicates() {
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
					else {
						cout << "Ok" << endl;
						cout << _data[ridx_out] << " / " <<  _data[ridx_in] << endl;
					}
				}
				ridx_in=(ridx_in+1) % _size;
			}

			ridx_out=(ridx_out+1) % _size;
		}
		return dubs;
	}*/

	void clear() {
//#ifdef USE_VECQUEUE
		_data.clear();
/*#else
		for(int i=0; i < _size; i++) {
			for(int j=0; j < 152; j++) {
				_data[i][j] = 0.0;
			}
		}
#endif*/
		_write_pointer=0;
		_read_pointer=0;
		_elems=0;
		_thresval=0.25 * 0;
	}

	/*void print() {
		cout << "Printing queue:" << endl;
		int i = _read_pointer;
		while(i != _write_pointer) {

			if(i >= 0) {
				cout << _data[i]  << endl;;
			}
			i=(i+1) % _size;
		}
		cout << "End of queue:" << endl;
	}

	void print(const int start, const int end, const int jj, const int kk) {
		cout << "Printing queue:" << endl;
		for(int i=start; i<=end;i++) {
			for(int j=jj; j<=kk; j++) {
				cout << _data[i][j] << " ";
			}
			cout << endl;

		}
		cout << "End of queue:" << endl;
	}*/

private:
//#ifdef USE_VECQUEUE
	std::vector<MBVec<T> > _data;
/*#else
	T** _data;
#endif*/


	int _write_pointer;
	int _read_pointer;
	int _size;
	int _elems;
	int _thresval;
};
}

#endif /* SRC_MBQUEUE_HPP_ */
