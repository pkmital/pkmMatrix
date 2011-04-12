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

#ifndef DEBUG
#define DEBUG 1
#endif

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
		//		pkm::Mat a(rhs);
		Mat(const Mat &rhs);
		Mat operator=(const Mat &rhs);
		
		void reset(size_t r, size_t c, bool clear = false)
		{
			
			rows = r;
			cols = c;
			
			if (bUserData) {
				data = (float *)malloc(rows * cols * sizeof(float));
				
				// sacrifice memory w/ speed, by pre-allocating a temporary buffer
				temp_data = (float *)malloc(rows * cols * sizeof(float));
			}
			else {
				
				data = (float *)realloc(data, rows * cols * sizeof(float));
				
				// sacrifice memory w/ speed, by pre-allocating a temporary buffer
				temp_data = (float *)realloc(temp_data, rows * cols * sizeof(float));
			}
			bAllocated = true;
			bUserData = false;
			
			// set every element to 0
			if(clear)
			{
				vDSP_vclr(data, 1, rows*cols);
			}
		}
		
		// set every element to a value
		inline void setTo(float val)
		{
#ifdef DEBUG
			assert(data != NULL);
#endif	
			vDSP_vfill(&val, data, 1, rows * cols);
		}	
		
		/////////////////////////////////////////
		
		inline float * row(size_t r)
		{
#ifdef DEBUG
			assert(data != NULL);
#endif			
			return (data + r*cols);
		}
		
		// inclusive of start, exclusive of end
		// can be a copy of the original matrix, or a way of editing the original
		// one by not copying the values (default)
		inline Mat rowRange(size_t start, size_t end, bool withCopy = true)
		{
#ifdef DEBUG
			assert(rows >= end);
#endif
			Mat submat(end-start, cols, row(start), withCopy);
			return submat;
		}
		
		
		inline Mat colRange(size_t start, size_t end, bool withCopy = true)
		{
#ifdef DEBUG
			assert(cols >= end);
#endif
			setTranspose();
			Mat submat(end-start, cols, row(start), withCopy);
			setTranspose();
			submat.setTranspose();
			return submat;
		}
		
		// copy data into the matrix
		inline void copy(Mat rhs)
		{
#ifdef DEBUG
			assert(rhs.rows == rows);
			assert(rhs.cols == cols);
#endif
			cblas_scopy(rows*cols, rhs.data, 1, data, 1);
		}
		
		/////////////////////////////////////////
		
		// element-wise multiplication
		inline void multiply(const Mat &rhs, Mat &result) const 
		{
#ifdef DEBUG
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
		// result stored in newly created matrix
		inline Mat multiply(const Mat &rhs)
		{
#ifdef DEBUG
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows && 
				   cols == rhs.cols);
#endif			
			Mat multiplied_matrix(rows, cols);
			
			vDSP_vmul(data, 1, rhs.data, 1, multiplied_matrix.data, 1, rows*cols);
			return multiplied_matrix;
		}		
		
		inline void multiply(float scalar, Mat &result) const 
		{
#ifdef DEBUG
			assert(data != NULL);
			assert(result.data != NULL);
			assert(rows == result.rows &&
				   cols == result.cols);
#endif			
			vDSP_vsmul(data, 1, &scalar, result.data, 1, rows*cols);
			
		}
		
		inline void multiply(float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif
			vDSP_vsmul(data, 1, &scalar, data, 1, rows*cols);
		}
		
		inline Mat operator/(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows && 
				   cols == rhs.cols);
#endif			
			Mat result(rows, cols);
			vDSP_vdiv(rhs.data, 1, data, 1, result.data, 1, rows*cols);
			return result;
			
		}
		
		inline Mat operator/(float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif			
			Mat result(rows, cols);
			vDSP_vsdiv(data, 1, &scalar, result.data, 1, rows*cols);
			return result;
			
		}
		
		
		inline void divide(const Mat &rhs, Mat &result) const 
		{
#ifdef DEBUG			
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
		
		inline void divide(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows &&
				   cols == rhs.cols);
#endif			
			vDSP_vdiv(rhs.data, 1, data, 1, temp_data, 1, rows*cols);
			std::swap(data, temp_data);
		}
		
		inline void divide(float scalar, Mat &result) const 
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(result.data != NULL);
			assert(rows == result.rows &&
				   cols == result.cols);
#endif	
			
			vDSP_vsdiv(data, 1, &scalar, result.data, 1, rows*cols);
		}
		
		inline void divide(float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif
			vDSP_vsdiv(data, 1, &scalar, data, 1, rows*cols);
		}
		
		inline void add(const Mat &rhs, Mat &result) const 
		{
#ifdef DEBUG			
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
		
		inline void add(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows &&
				   cols == rhs.cols);
#endif			
			vDSP_vadd(data, 1, rhs.data, 1, temp_data, 1, rows*cols);
			std::swap(data, temp_data);
		}
		
		inline void subtract(const Mat &rhs, Mat &result) const 
		{
#ifdef DEBUG			
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
		
		inline void subtract(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows &&
				   cols == rhs.cols);
#endif			
			vDSP_vsub(data, 1, rhs.data, 1, temp_data, 1, rows*cols);
			std::swap(data, temp_data);
		}
		
		
		inline void GEMM(const Mat &rhs, Mat &result) const 
		{
#ifdef DEBUG
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
		
		inline Mat GEMM(const pkm::Mat &rhs)
		{
#ifdef DEBUG
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(cols == rhs.rows);
#endif
			
			Mat gemmResult(rows, rhs.cols);
			
			printf("lda: %d\nldb: %d\nldc: %d\n", rows, rhs.rows, gemmResult.rows); 
			cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, gemmResult.rows, gemmResult.cols, cols, 1.0f, data, cols, rhs.data, rhs.cols, 0.0f, gemmResult.data, gemmResult.cols);
			//vDSP_mmul(data, 1, rhs.data, 1, gemmResult.data, 1, gemmResult.rows, gemmResult.cols, cols);
			return gemmResult;
			
		}
		
		inline Mat operator*(const pkm::Mat &rhs)
		{
#ifdef DEBUG
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(cols == rhs.rows);
#endif
			
			Mat gemmResult(rows, rhs.cols);
			//ldb must be >= MAX(N,1): ldb=30 N=3533Parameter 11 to routine cblas_sgemm was incorrect
			cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, gemmResult.rows, gemmResult.cols, cols, 1.0f, data, cols, rhs.data, rhs.cols, 0.0f, gemmResult.data, gemmResult.cols);
			//vDSP_mmul(data, 1, rhs.data, 1, gemmResult.data, 1, gemmResult.rows, gemmResult.cols, cols);
			return gemmResult;
		}
		
		inline Mat operator*(float scalar)
		{
#ifdef DEBUG
			assert(data != NULL);
#endif
			
			Mat gemmResult(rows, cols);
			vDSP_vsmul(data, 1, &scalar, gemmResult.data, 1, rows*cols);

			return gemmResult;
		}
		
		inline void setTranspose()
		{
#ifdef DEBUG      
			assert(data != NULL);
#endif      
			vDSP_mtrans(data, 1, temp_data, 1, cols, rows);
			//cblas_scopy(rows*cols, temp_data, 1, data, 1);
			std::swap(data, temp_data);					// swap will break certain operations for col/row range as their pointers will have changed. :(
			std::swap(rows, cols);
		}
		
		Mat getTranspose();
		
		
		// diagonalize the vector into a square matrix with 
		// the current data vector along the diagonal
		inline void setDiag()
		{
#ifdef DEBUG
			assert(data != NULL);
#endif	
			if(rows == 1 && cols > 1 || cols == 1 && rows > 1)
			{
				int diagonal_elements = MAX(rows,cols);
				
				// create a square matrix
				temp_data = (float *)realloc(temp_data, diagonal_elements*diagonal_elements*sizeof(float));
				
				// set values to 0
				vDSP_vclr(temp_data, 1, diagonal_elements*diagonal_elements);
				
				// set diagonal elements to the current vector in data
				for (int i = 0; i < diagonal_elements; i++) {
					temp_data[i*diagonal_elements+i] = data[i];
				}
				
				// store in data
				std::swap(data, temp_data);
				
				// reallocate temp data for future processing
				temp_data = (float *)realloc(temp_data, diagonal_elements*diagonal_elements*sizeof(float));
				
				// save dimensions
				rows = cols = diagonal_elements;
			}
		}
		Mat getDiag();
		
		// returns a new diagonalized matrix version of A
		static Mat diag(Mat &A);
		
		// get a new identity matrix of size dim x dim
		static Mat identity(size_t dim);
		
		static Mat zeros(size_t rows, size_t cols)
		{
			return Mat(rows, cols, true);
		}
		
		// set every element to a random value between low and high
		void setRand(float low = 0.0, float high = 1.0);
		
		// create a random matrix
		static Mat rand(size_t r, size_t c, float low = 0.0, float high = 1.0);
		
		// sum across rows or columns creating a vector from a matrix, or a scalar from a vector
		Mat sum(bool across_rows = true);
		
		// repeat a vector for size times
		static Mat repeat(Mat &m, size_t size)
		{
			// repeat a column vector across cols
			if(m.rows > 1 && m.cols == 1 && size > 1)
			{
				Mat repeated_matrix(size, m.rows);
				for (int i = 0; i < size; i++) {
					cblas_scopy(m.rows, m.data, 1, repeated_matrix.data + (i*m.rows), 1);
				}
				repeated_matrix.setTranspose();
				return repeated_matrix;
			}
			else if( m.rows == 1 && m.cols > 1 && size > 1)
			{
				Mat repeated_matrix(size, m.cols, 5.0f);
				
				for (int i = 0; i < size; i++) {
					cblas_scopy(m.cols, m.data, 1, repeated_matrix.data + (i*m.cols), 1);
				}
				return repeated_matrix;
			}
			else {
				printf("[ERROR]: repeat requires a vector and a size to repeat on.");
				Mat a;
				return a;
			}

		}
		
		// repeat a vector for size times
		static void repeat(Mat &dst, const Mat &m, size_t size)
		{
			// repeat a column vector across cols
			if(m.rows > 1 && m.cols == 1 && size > 1)
			{
				dst.reset(size, m.rows);
				for (int i = 0; i < size; i++) {
					cblas_scopy(m.rows, m.data, 1, dst.data + (i*m.rows), 1);
				}
				dst.setTranspose();
			}
			else if( m.rows == 1 && m.cols > 1 && size > 1)
			{
				dst.reset(size, m.cols);
				
				for (int i = 0; i < size; i++) {
					cblas_scopy(m.cols, m.data, 1, dst.data + (i*m.cols), 1);
				}
			}
			else {
				printf("[ERROR]: repeat requires a vector and a size to repeat on.");
				
			}
			
		}
		
		// rescale the values in each row to their maximum
		void setNormalize(bool row_major = true);
		
		// simple print output (be careful with large matrices!)
		void print(bool row_major = true);
		
		/////////////////////////////////////////
		
		size_t rows;
		size_t cols;
		
		float *data;
		float *temp_data;
		
		bool bAllocated;
		bool bUserData;
	};
	
};