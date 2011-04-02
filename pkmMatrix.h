/*
 *  pkmMatrix.h
 *  
 
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

#pragma once

#include <iostream>
#include <assert.h>
#include <Accelerate/Accelerate.h>
//#include <Accelerate/cblas.h>

#ifndef MAX
#define MAX(a,b)  ((a) < (b) ? (b) : (a))
#endif

namespace pkm
{	
	// row-major floating point matrix
	class Mat
	{
		/////////////////////////////////////////
	public:
		// default constructor
		Mat();
		
		// destructor
		~Mat();
		
		// allocate data
		Mat(size_t r, size_t c, bool clear = false);
		
		// pass in existing data
		// non-destructive by default
		Mat(size_t r, size_t c, float *existing_buffer, bool withCopy = true);
		
		// set every element to a value
		Mat(size_t r, size_t c, float val);
		
		// copy-constructor, called during:
		//		pkm::Mat a = rhs;
		//		pkm::Mat a(rhs);
		Mat(const Mat &rhs);
		
		// set every element to a value
		void setTo(float val);		
		
		/////////////////////////////////////////
		
		inline float *row(size_t r);
		
		/////////////////////////////////////////
		
		// element-wise multiplication
		void multiply(Mat rhs, Mat &result);		
		void multiply(Mat rhs);		
		void multiply(float scalar, Mat &result);		
		void multiply(float scalar);
		
		void divide(Mat rhs, Mat &result);		
		void divide(Mat rhs);		
		void divide(float scalar, Mat &result);		
		void divide(float scalar);
		
		void add(Mat rhs, Mat &result);		
		void add(Mat rhs);
		
		void subtract(Mat rhs, Mat &result);		
		void subtract(Mat rhs);
		
		void GEMM(Mat rhs, Mat &result);
		Mat GEMM(Mat rhs);
		
		void setTranspose();		
		Mat getTranspose();
		
		// diagonalize the vector into a square matrix with 
		// the current data vector along the diagonal
		void setDiag();
		Mat getDiag();
		
		// returns a new diagonalized matrix version of A
		static Mat diag(Mat &A);
		
		// get a new identity matrix of size dim x dim
		static Mat identity(size_t dim);
		
		// set every element to a random value between low and high
		void setRand(float low, float high);
		
		// create a random matrix
		static Mat rand(size_t r, size_t c, float low = 0.0, float high = 1.0);
		
		// rescale the values in each row to their maximum
		void setNormalize(bool row_major = true);
		
		// simple print output (be careful with large matrices!)
		void print(bool row_major = true);
		
		/////////////////////////////////////////
		
		size_t rows;
		size_t cols;
		
		float *data;
		float *temp_data;
	};
	
};