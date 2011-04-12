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
	
	
	size_t data_size = 15;
	size_t init_components = 5;
	size_t components = 10;
	
	Mat M_initpw = Mat(data_size, init_components, 1.0f);
	printf("M_initpw:\n");
	M_initpw.print();
	
	Mat M_pw(4,3, 5.0f);
	
	M_pw.reset(data_size, components, true);
	M_pw.setTranspose();
	Mat M_pw1 = M_pw.rowRange(0, M_initpw.cols, false);
	M_pw1.copy(M_initpw.getTranspose());
	M_pw.print();
	printf("Resizing w: %d x %d to %d x %d\n", M_initpw.rows, M_initpw.cols, M_pw.rows, M_pw.cols);
	M_pw.rowRange(M_initpw.cols, components, false).setRand();
	M_pw.setTranspose();
	M_pw.print();
	
	return 0;
}
