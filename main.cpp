/*
 *  main.cpp
 *  
 
 main.cpp driver showing how to use pkmMatrix
 row-major floating point matrix utility class
 utilizes Apple Accelerate's vDSP functions for SSE optimizations
 
 Copyright (C) 2011 Parag K. Mital
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 *
 */

#include <iostream>
#include "pkmMatrix.h"

using namespace pkm;
using namespace std;

int main (int argc, char * const argv[]) {
    
	// constructing identity matrix
	Mat idMatrix = pkm::Mat::identity(5);
	idMatrix.print();
	if(1)
	{
		int rows = 3;
		int cols = 4;
		Mat B(rows, cols, 5.0f);
		Mat B2(cols, rows, 10.0f);
		B2.print();
		
		for (int i = 0; i < rows*cols; i++) {
			B.data[i] = i+1;
		}
		B.print();
		
		// Setting a column to random elements:
		B.setTranspose();
		B.print();
		for (int i = 0; i < 3; i++) {
			printf("%f, ", B.row(3)[i]);
		}
		printf("\n");
		Mat B_range = B.rowRange(2,3,false);
		printf("B_range:\n");
		B_range.print();
		B_range.setRand();
		B.setTranspose();
		B.print();

		// creating a diagonal matrix
		B.print();
		B.setTranspose();
		Mat B3_range = B.rowRange(2,3,false);
		Mat B2_range = B2.rowRange(2,3,false);
		B3_range.print();
		B2_range.setNormalize();
		B2_range.print();
		B3_range.copy(B2_range);
		B.setTranspose();
		B.print();
	}
	
	Mat C(5,5);
	
	return 0;
}
