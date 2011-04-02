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
		inline void setTo(float val)
		{
#ifndef DEBUG
			assert(data != NULL);
#endif	
			vDSP_vfill(&val, data, 1, rows * cols);
		}	
		
		/////////////////////////////////////////
		
		inline float * row(size_t r)
		{
#ifndef DEBUG
			assert(data != NULL);
#endif			
			return (data + r*cols);
		}
		
		inline Mat rowRange(size_t start, size_t end)
		{
			Mat submat(end-start, cols, row(start));
			return submat;
		}
		
		/////////////////////////////////////////
		
		// element-wise multiplication
		inline void multiply(Mat rhs, Mat &result)
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
		inline void multiply(Mat rhs)
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
		
		inline void multiply(float scalar, Mat &result)
		{
#ifndef DEBUG
			assert(data != NULL);
			assert(result.data != NULL);
			assert(rows == result.rows &&
				   cols == result.cols);
#endif			
			vDSP_vsmul(data, 1, &scalar, result.data, 1, rows*cols);
			
		}
		
		inline void multiply(float scalar)
		{
#ifndef DEBUG			
			assert(data != NULL);
#endif
			vDSP_vsmul(data, 1, &scalar, data, 1, rows*cols);
		}
		
		inline void divide(Mat rhs, Mat &result)
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
		
		inline void divide(Mat rhs)
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
		
		inline void divide(float scalar, Mat &result)
		{
#ifndef DEBUG			
			assert(data != NULL);
			assert(result.data != NULL);
			assert(rows == result.rows &&
				   cols == result.cols);
#endif	
			
			vDSP_vsdiv(data, 1, &scalar, result.data, 1, rows*cols);
		}
		
		inline void divide(float scalar)
		{
#ifndef DEBUG			
			assert(data != NULL);
#endif
			vDSP_vsdiv(data, 1, &scalar, data, 1, rows*cols);
		}
		
		inline void add(Mat rhs, Mat &result)
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
		
		inline void add(Mat rhs)
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
		
		inline void subtract(Mat rhs, Mat &result)
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
		
		inline void subtract(Mat rhs)
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
		
		
		inline void GEMM(Mat rhs, Mat &result)
		{
#ifndef DEBUG
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(result.data != NULL);
			assert(rows == result.rows &&
				   rhs.cols == result.cols &&
				   cols == rhs.rows);
#endif
			
			cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, result.rows, result.cols, cols, 1.0f, data, rows, rhs.data, rhs.rows, 0.0f, result.data, result.cols);
			//vDSP_mmul(data, 1, rhs.data, 1, result.data, 1, result.rows, result.cols, cols);
			
		}
		
		inline Mat GEMM(Mat rhs)
		{
#ifndef DEBUG
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(cols == rhs.rows);
#endif
			
			Mat gemmResult(rows, rhs.cols);
			
			cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, gemmResult.rows, gemmResult.cols, cols, 1.0f, data, cols, rhs.data, rhs.rows, 0.0f, gemmResult.data, gemmResult.cols);
			//vDSP_mmul(data, 1, rhs.data, 1, gemmResult.data, 1, gemmResult.rows, gemmResult.cols, cols);
			return gemmResult;
			
		}
		
		
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