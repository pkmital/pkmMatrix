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
		Mat()
		{
			rows = cols = 0;
			data = temp_data = NULL;
		}
		
		// destructor
		~Mat()
		{
			free(data);
			free(temp_data);
		}
		
		// allocate data
		Mat(size_t r, size_t c, bool clear = false)
		{
			rows = r;
			cols = c;
			data = (float *)malloc(rows * cols * sizeof(float));
			
			// sacrifice memory w/ speed, by pre-allocating a temporary buffer
			temp_data = (float *)malloc(rows * cols * sizeof(float));
			
			// set every element to 0
			if(clear)
			{
				vDSP_vclr(data, 1, rows*cols);
			}
		}
		
		// pass in existing data
		// non-destructive by default
		Mat(size_t r, size_t c, float *existing_buffer, bool withCopy = true)
		{
			rows = r;
			cols = c;
			if(withCopy)
			{
				cblas_scopy(rows*cols, existing_buffer, 1, data, 1);
				//memcpy(data, existing_buffer, sizeof(float)*r*c);
			}
			else {
				data = existing_buffer;
			}
			
		}
		
		// set every element to a value
		Mat(size_t r, size_t c, float val)
		{
			rows = r;
			cols = c;
			data = (float *)malloc(rows * cols * sizeof(float));
			// sacrifice memory w/ speed, by pre-allocating a temporary buffer
			temp_data = (float *)malloc(rows * cols * sizeof(float));
			
			// set every element to val
			vDSP_vfill(&val, data, 1, rows * cols);
			
		}
		
		// copy-constructor, called during:
		//		pkm::Mat a = rhs;
		//		pkm::Mat a(rhs);
		Mat(const Mat &rhs)
		{
			if(rhs.data != NULL)
			{
				rows = rhs.rows;
				cols = rhs.cols;
				
				data = (float *)malloc(rows * cols * sizeof(float));
				
				// sacrifice memory w/ speed, by pre-allocating a temporary buffer
				temp_data = (float *)malloc(rows * cols * sizeof(float));
				
				cblas_scopy(rows*cols, rhs.data, 1, data, 1);
				//memcpy(data, rhs.data, sizeof(float)*rows*cols);
			}
			else {
				rows = 0;
				cols = 0;
				data = NULL;
				temp_data = NULL;
			}
		}
		
		
		void setTo(float val)
		{
#ifndef DEBUG
			assert(data != NULL);
#endif	
			vDSP_vfill(&val, data, 1, rows * cols);
		}
		
		
		/////////////////////////////////////////
		
		float *row(size_t r)
		{
#ifndef DEBUG
			assert(data != NULL);
#endif			
			return (data + r*cols);
		}
		
		/////////////////////////////////////////
		
		void multiply(Mat rhs, Mat &result)
		{
#ifndef DEBUG
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(result.data != NULL);
			assert(rows == rhs.rows && 
				   rhs.rows == result.rows &&
				   cols == rhs.cols && 
				   rhs.cols == result.cols);
#endif
			vDSP_vmul(data, 1, rhs.data, 1, result.data, 1, rows*cols);
			
		}
		
		// element-wise multiplication
		// result stored in original matrix
		void multiply(Mat rhs)
		{
#ifndef DEBUG
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows && 
				   cols == rhs.cols);
#endif			
			vDSP_vmul(data, 1, rhs.data, 1, temp_data, 1, rows*cols);
			std::swap(data, temp_data);
		}
		
		
		void multiply(float scalar, Mat &result)
		{
#ifndef DEBUG
			assert(data != NULL);
			assert(result.data != NULL);
			assert(rows == result.rows &&
				   cols == result.cols);
#endif			
			vDSP_vsmul(data, 1, &scalar, result.data, 1, rows*cols);
			
		}
		
		void multiply(float scalar)
		{
#ifndef DEBUG			
			assert(data != NULL);
#endif
			vDSP_vsmul(data, 1, &scalar, data, 1, rows*cols);
		}
		
		void divide(Mat rhs, Mat &result)
		{
#ifndef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(result.data != NULL);
			assert(rows == rhs.rows && 
				   rhs.rows == result.rows &&
				   cols == rhs.cols && 
				   rhs.cols == result.cols);
#endif			
			vDSP_vdiv(rhs.data, 1, data, 1, result.data, 1, rows*cols);
			
		}
		
		void divide(Mat rhs)
		{
#ifndef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows &&
				   cols == rhs.cols);
#endif			
			vDSP_vdiv(rhs.data, 1, data, 1, temp_data, 1, rows*cols);
			std::swap(data, temp_data);
		}
		
		void divide(float scalar, Mat &result)
		{
#ifndef DEBUG			
			assert(data != NULL);
			assert(result.data != NULL);
			assert(rows == result.rows &&
				   cols == result.cols);
#endif	
			
			vDSP_vsdiv(data, 1, &scalar, result.data, 1, rows*cols);
		}
		
		void divide(float scalar)
		{
#ifndef DEBUG			
			assert(data != NULL);
#endif
			vDSP_vsdiv(data, 1, &scalar, data, 1, rows*cols);
		}
		
		void add(Mat rhs, Mat &result)
		{
#ifndef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(result.data != NULL);
			assert(rows == rhs.rows && 
				   rhs.rows == result.rows &&
				   cols == rhs.cols && 
				   rhs.cols == result.cols);
#endif			
			vDSP_vadd(data, 1, rhs.data, 1, result.data, 1, rows*cols);
		}
		
		void add(Mat rhs)
		{
#ifndef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows &&
				   cols == rhs.cols);
#endif			
			vDSP_vadd(data, 1, rhs.data, 1, temp_data, 1, rows*cols);
			std::swap(data, temp_data);
		}
		
		void subtract(Mat rhs, Mat &result)
		{
#ifndef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(result.data != NULL);
			assert(rows == rhs.rows && 
				   rhs.rows == result.rows &&
				   cols == rhs.cols && 
				   rhs.cols == result.cols);
#endif			
			vDSP_vsub(data, 1, rhs.data, 1, data, 1, rows*cols);
			
		}
		
		void subtract(Mat rhs)
		{
#ifndef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows &&
				   cols == rhs.cols);
#endif			
			vDSP_vsub(data, 1, rhs.data, 1, temp_data, 1, rows*cols);
			std::swap(data, temp_data);
		}
		
		
		
		void GEMM(Mat rhs, Mat &result);
		Mat GEMM(Mat rhs);
		
		void setTranspose();
		
		Mat getTranspose();
		
		// diagonalize the vector into a square matrix with 
		// the current data vector along the diagonal
		void setDiag();
		
		// get a diagonalized version of the current vector (non-destructive)
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
		
		// simple output
		void print(bool row_major = true);
		
		/////////////////////////////////////////
		
		size_t rows;
		size_t cols;
		
		float *data;
		float *temp_data;
	};
	
};