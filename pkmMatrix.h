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
#include <vector>
using namespace std;
#ifndef DEBUG
#define DEBUG 1
#endif

#ifndef MAX
#define MAX(a,b)  ((a) < (b) ? (b) : (a))
#endif

#ifndef MIN
#define MIN(a,b)  ((a) > (b) ? (b) : (a))
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
		
        Mat(vector<float> m);
        
		// allocate data
		Mat(int r, int c, bool clear = false);
		
		// pass in existing data
		// non-destructive by default
		Mat(int r, int c, float *existing_buffer, bool withCopy);
		
		// set every element to a value
		Mat(int r, int c, float val);
		
		// copy-constructor, called during:
		//		pkm::Mat a(rhs);
		Mat(const Mat &rhs);
		Mat & operator=(const Mat &rhs);
		Mat & operator=(const vector<float> &rhs);
		
		
		inline Mat operator+(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows &&
				   cols == rhs.cols);
#endif			
			Mat newMat(rows, cols);
			vDSP_vadd(data, 1, rhs.data, 1, newMat.data, 1, rows*cols);
			//cblas_scopy(rows*cols, temp_data, 1, data, 1);
			return newMat;
		}
		
		
		inline Mat operator+(float rhs)
		{	
#ifdef DEBUG			
			assert(data != NULL);
#endif			
			Mat newMat(rows, cols);
			vDSP_vsadd(data, 1, &rhs, newMat.data, 1, rows*cols);
			return newMat;
		}
		
		inline Mat operator-(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows &&
				   cols == rhs.cols);
#endif			
			Mat newMat(rows, cols);
			vDSP_vsub(data, 1, rhs.data, 1, newMat.data, 1, rows*cols);
			//cblas_scopy(rows*cols, temp_data, 1, data, 1);
			return newMat;
		}
		

		
		inline Mat operator-(const float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif			
			Mat newMat(rows, cols);
			float rhs = -scalar;
			vDSP_vsadd(data, 1, &rhs, newMat.data, 1, rows*cols);
			return newMat;
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
		
		inline Mat operator>(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows && 
				   cols == rhs.cols);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] > rhs.data[i];
			return result;
		}
		
		inline Mat operator>(float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] > scalar;
			return result;
		}
		
		inline Mat operator>=(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows && 
				   cols == rhs.cols);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] >= rhs.data[i];
			return result;
		}
		
		inline Mat operator>=(float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] >= scalar;
			return result;
		}
		
		inline Mat operator<(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows && 
				   cols == rhs.cols);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] < rhs.data[i];
			return result;
		}
		
		inline Mat operator<(float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] < scalar;
			return result;
		}
		
		inline Mat operator<=(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows && 
				   cols == rhs.cols);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] <= rhs.data[i];
			return result;
		}
		
		inline Mat operator<=(float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] <= scalar;
			return result;
		}
		
		inline Mat operator==(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows && 
				   cols == rhs.cols);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] == rhs.data[i];
			return result;
		}
		
		inline Mat operator==(float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] == scalar;
			return result;
		}
		
		inline Mat operator!=(const Mat &rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows && 
				   cols == rhs.cols);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] != rhs.data[i];
			return result;
		}
		
		inline Mat operator!=(float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif			
			Mat result(rows, cols);
			for(int i = 0; i < rows*cols; i++)
				result.data[i] = data[i] != scalar;
			return result;
		}
		
		inline float & operator[](int idx)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rows*cols >= idx);
#endif	
			return data[idx];
		}
		
		// return a vector composed on non-zero indices of logicalMat
		inline Mat operator[](Mat rhs)
		{
#ifdef DEBUG			
			assert(data != NULL);
			assert(rhs.data != NULL);
			assert(rows == rhs.rows && 
				   cols == rhs.cols);
#endif	
			std::vector<float> newMat;
			int count = 0;
			for(int i = 0; i < rows*cols; i++)
			{
				if (rhs.data[i] > 0) {
					newMat.push_back(data[i]);
				}
			}
			if (newMat.size() > 0) {
				Mat result(1,newMat.size());
				for(int i = 0; i < newMat.size(); i++)
				{
					result.data[i] = newMat[i];
				}
				return result;
			}
			else {
				Mat empty;
				return empty;
			}
		}
		

		
		friend Mat operator-(float lhs, const Mat &rhs)
		{
#ifdef DEBUG			
			assert(rhs.data != NULL);
#endif			
			Mat newMat(rhs.rows, rhs.cols);
			float scalar = -lhs;
			vDSP_vsadd(rhs.data, 1, &scalar, newMat.data, 1, rhs.rows*rhs.cols);
			return newMat;
		}
		
		friend Mat operator*(float lhs, const Mat &rhs)
		{
#ifdef DEBUG
			assert(rhs.data != NULL);
#endif
			
			Mat gemmResult(rhs.rows, rhs.cols);
			vDSP_vsmul(rhs.data, 1, &lhs, gemmResult.data, 1, rhs.rows*rhs.cols);
			
			return gemmResult;
		}
		friend Mat operator+(float lhs, const Mat &rhs)
		{
#ifdef DEBUG			
			assert(rhs.data != NULL);
#endif			
			Mat newMat(rhs.rows, rhs.cols);
			vDSP_vsadd(rhs.data, 1, &lhs, newMat.data, 1, rhs.rows*rhs.cols);
			//cblas_scopy(rows*cols, temp_data, 1, data, 1);
			return newMat;
		}
		
		
		// can be used to create an already declared matrix without a copy constructor
		void reset(int r, int c, bool clear = false)
		{
			rows = r;
			cols = c;
			current_row = 0;
			bCircularInsertionFull = false;
			
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
		
		inline void clear()
		{
			if (rows == 0 || cols == 0) {
				return;
			}

			vDSP_vclr(data, 1, rows * cols);
		}
		
		/////////////////////////////////////////
		
		inline float * row(int r)
		{
#ifdef DEBUG
			assert(data != NULL);
#endif			
			return (data + r*cols);
		}
		
		inline void insertRow(float *buf, int row_idx)
		{
			float * rowData = row(row_idx);
			cblas_scopy(cols, buf, 1, rowData, 1);
		}
		
        inline void push_back(Mat m)
        {
            if (rows > 0 && cols > 0) {
                if (m.cols != cols) {
                    printf("[ERROR]: pkm::Mat push_back(Mat m) requires same number of columns!\n");
                    return;
                }
                cblas_scopy(rows*cols, data, 1, temp_data, 1);
                data = (float *)realloc(data, (rows+m.rows)*cols*sizeof(float));
                cblas_scopy(rows*cols, temp_data, 1, data, 1);
                cblas_scopy(m.rows*cols, m.data, 1, data + (rows*cols), 1);
                rows+=m.rows;
                temp_data = (float *)realloc(temp_data, rows*cols*sizeof(float));
            }
            else {
                *this = m;
            }

        }
        
        inline void push_back(vector<float> m)
        {
            if (rows > 0 && cols > 0) {
                if (m.size() != cols) {
                    printf("[ERROR]: pkm::Mat push_back(vector<float> m) requires same number of columns in Mat as length of vector!\n");
                    return;
                }
                cblas_scopy(rows*cols, data, 1, temp_data, 1);
                data = (float *)realloc(data, (rows+1)*cols*sizeof(float));
                cblas_scopy(rows*cols, temp_data, 1, data, 1);
                cblas_scopy(cols, &(m[0]), 1, data + (rows*cols), 1);
                rows++;
                temp_data = (float *)realloc(temp_data, rows*cols*sizeof(float));
            }
            else {
                *this = m;
            }
            
        }
        
		inline void insertRowCircularly(float *buf)
		{
			insertRow(buf, current_row);
			current_row = (current_row + 1) % rows;
			if (current_row == 0) {
				bCircularInsertionFull = true;
			}
		}
		
		// inclusive of start, exclusive of end
		// can be a copy of the original matrix, or a way of editing the original
		// one by not copying the values (default)
		inline Mat rowRange(int start, int end, bool withCopy = true)
		{
#ifdef DEBUG
			assert(rows >= end);
#endif
			Mat submat(end-start, cols, row(start), withCopy);
			return submat;
		}
		
		
		inline Mat colRange(int start, int end, bool withCopy = true)
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
		void copy(Mat rhs)
		{
#ifdef DEBUG
			assert(rhs.rows == rows);
			assert(rhs.cols == cols);
#endif
			cblas_scopy(rows*cols, rhs.data, 1, data, 1);
			
		}
		
		void copy(Mat &rhs, Mat &indx)
		{
#ifdef DEBUG
			assert(indx.rows == rows);
			assert(indx.cols == cols);
#endif
			int idx = 0;
			for(int i = 0; i < rows; i++)
			{
				for(int j = 0; j < cols; j++)
				{
					if (indx.data[i*cols + j]) {
						data[i*cols + j] = rhs[idx];
						idx++;
					}
				}
				
			}
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
			//cblas_scopy(rows*cols, temp_data, 1, data, 1);
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
			//cblas_scopy(rows*cols, temp_data, 1, data, 1);
			std::swap(data, temp_data);
		}
		
		inline void add(float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif			
			vDSP_vsadd(data, 1, &scalar, data, 1, rows*cols);
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
			//cblas_scopy(rows*cols, temp_data, 1, data, 1);
			std::swap(data, temp_data);
		}
		
		inline void subtract(float scalar)
		{
#ifdef DEBUG			
			assert(data != NULL);
#endif			
			float rhs = -scalar;
			vDSP_vsadd(data, 1, &rhs, data, 1, rows*cols);
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
				rows = cols = diagonal_elements;
				//free(data);
				//data = (float *)malloc(rows*cols*sizeof(float));
				//cblas_scopy(rows*cols, temp_data, 1, data, 1);
				std::swap(data, temp_data);
				
				// reallocate temp data for future processing
				temp_data = (float *)realloc(temp_data, diagonal_elements*diagonal_elements*sizeof(float));
				
				// save dimensions
			}
		}
		Mat getDiag();
		
		// returns a new matrix with each el the log(el)
		static Mat log(Mat &A);
		
		// returns a new matrix with each el the exp(el)
		static Mat exp(Mat &A);
		
		// returns a new diagonalized matrix version of A
		static Mat diag(Mat &A);
		
		// get a new identity matrix of size dim x dim
		static Mat identity(int dim);
		
		static Mat zeros(int rows, int cols)
		{
			return Mat(rows, cols, true);
		}
		
		// set every element to a random value between low and high
		void setRand(float low = 0.0, float high = 1.0);
		
		// create a random matrix
		static Mat rand(int r, int c, float low = 0.0, float high = 1.0);
		
		// sum across rows or columns creating a vector from a matrix, or a scalar from a vector
		Mat sum(bool across_rows = true);
		
		// repeat a vector for size times
		static Mat repeat(Mat &m, int size)
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
		static void repeat(Mat &dst, const Mat &m, int size)
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
		
		static float meanMagnitude(float *buf, int size)
		{
			float mean;
			vDSP_meamgv(buf, 1, &mean, size);
			return mean;
		}
		
		static float sumOfAbsoluteDifferences(float *buf1, float *buf2, int size)
		{
			int a = size;
			float diff = 0;
			float *p1 = buf1, *p2 = buf2;
			while (a) {
				diff += fabs(*p1++ - *p2++);
				a--;
			}
			return diff/(float)size;
		}
		
		static float mean(float *buf, int size)
		{
			float val;
			vDSP_meanv(buf, 1, &val, size);
			return val;
		}
		
		static float var(float *buf, int size)
		{
			float m = mean(buf, size);
			float v = 0;
			float sqr = 0;
			float *p = buf;
			int a = size;
			while (a) {
				sqr = (*p++ - m);
				v += sqr*sqr;
				a--;
			}
			return v/(float)size;
		}
		
		static float rms(float *buf, int size)
		{
			float val;
			vDSP_rmsqv(buf, 1, &val, size);
			return val;
		}
		
		static float min(Mat &A)
		{
			float minval;
			vDSP_minv(A.data, 1, &minval, A.rows*A.cols);
			return minval;
		}
		
		static float max(Mat &A)
		{
			float maxval;
			vDSP_maxv(A.data, 1, &maxval, A.rows*A.cols);
			return maxval;
		}
		
		static float sum(Mat &A)
		{
			float sumval;
			vDSP_sve(A.data, 1, &sumval, A.rows*A.cols);
			return sumval;
		}
		
		
		// rescale the values in each row to their maximum
		void setNormalize(bool row_major = true);
		
		void normalizeRow(int r)
		{
			float min, max;
			vDSP_minv(&(data[r*cols]), 1, &min, cols);
			vDSP_maxv(&(data[r*cols]), 1, &max, cols);
			float height = max-min;
			min = -min;
			vDSP_vsadd(&(data[r*cols]), 1, &min, &(data[r*cols]), 1, cols);
			if (height != 0) {
				vDSP_vsdiv(&(data[r*cols]), 1, &height, &(data[r*cols]), 1, cols);	
			}
		}
		
		void divideEachVecByMaxVecElement(bool row_major);
		void divideEachVecBySum(bool row_major);
        
        
        bool save(string filename)
        {
            FILE *fp;
            fp = fopen(filename.c_str(), "w");
            fprintf(fp, "%d %d\n", rows, cols);
            for(int i = 0; i < rows; i++)
            {
                for(int j = 0; j < cols; j++)
                {
                    fprintf(fp, "%f, ", data[i*cols + j]);
                }
            }
            fclose(fp);
            return true;
        }
        
        bool load(string filename)
        {
            if (bAllocated) {
                free(data);
                free(temp_data);
                rows = cols = 0;
            }
            FILE *fp;
            fp = fopen(filename.c_str(), "r");
            fscanf(fp, "%d %d\n", &rows, &cols);
            data = (float *)malloc(sizeof(float) * rows * cols);
            temp_data = (float *)malloc(sizeof(float) * rows * cols);
            for(int i = 0; i < rows; i++)
            {
                for(int j = 0; j < cols; j++)
                {
                    fscanf(fp, "%f, ", &(data[i*cols + j]));
                }
            }
            fclose(fp);
            return true;
        }
		
		// simple print output (be careful with large matrices!)
		void print(bool row_major = true);
		// only prints maximum of 5 rows/cols
		void printAbbrev(bool row_major = true);
		
		/////////////////////////////////////////
		
		int current_row;	// for circular insertion
		bool bCircularInsertionFull;
		int rows;
		int cols;
		
		float *data;
		float *temp_data;
		
		bool bAllocated;
		bool bUserData;
	};
	
};