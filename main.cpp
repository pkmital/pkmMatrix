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
	printf("Mat idMatrix = pkm::Mat::identity(5):\n");
	idMatrix.print();
	
	// create a row-vector, and fill with the value 1.0f
	Mat a1 = Mat(1, 5, 1.0f);
	printf("Mat a1 = Mat(1, 5, 1.0f):\n");
	a1.print();
	
	// diagonalize a1, and put in a2
	Mat a2 = pkm::Mat::diag(a1);
	printf("Mat a2 = pkm::Mat::diag(a1):\n");
	a2.print();
	
	// create a 5x5 matrix and clear to 0;
	Mat a3(5, 5, true);
	printf("Mat a3(5, 5, true):\n");
	a3.print();
	
	// create random values between 0-1
	a3.setRand(0.0f, 1.0f);
	printf("a3.setRand(0.0f, 1.0f):\n");
	a3.print();
	
	// normalize each row-vector so the values scale 0-1
	a3.setNormalize();
	printf("a3.setNormalize():\n");
	a3.print();
	
	Mat a4 = Mat(1,5);
	a4.setTo(2.0);
	printf("a4.setTo(2.0):\n");
	a4.print();
	
	Mat a4T = a4.getTranspose();
	printf("a4'\n");
	a4T.print();
	

	Mat a5 = a3.GEMM(a4T);
	printf("Mat a5 = a3.GEMM(a4T):\n");
	a5.print();
	
	a4.data[0] = 10.0f;
	a4.divide(5.0f);
	printf("a4.divide(5.0f):\n");
	a4.print();
	
	Mat a6(5,1);
	a4T.divide(a5, a6);
	printf("a6 = a4T / a5\n");
	a6.print();
	
	a6.setNormalize(false);
	printf("a6.setNormalize():\n");
	a6.print();
	
	Mat a7 = a3.GEMM(idMatrix);
	printf("a4*idMatrix\n");
	a7.print();
	
	
	
	return 0;
}
