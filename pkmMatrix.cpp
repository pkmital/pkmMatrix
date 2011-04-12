/*
 *  pkmMatrix.cpp
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

#include "pkmMatrix.h"

using namespace pkm;

Mat::Mat()
{
	bUserData = false;
	rows = cols = 0;
	data = temp_data = NULL;
	bAllocated = false;
	
}

// destructor
Mat::~Mat()
{
	if(bAllocated && !bUserData && data != NULL)
	{
		free(data);
	}
	if(bAllocated && temp_data != NULL) 
	{
		free(temp_data);	
	}
	rows = cols = 0;
	data = temp_data = NULL;
	bAllocated = false;
}

// allocate data
Mat::Mat(size_t r, size_t c, bool clear)
{
	data = temp_data = NULL;
	
	bUserData = false;
	rows = r;
	cols = c;
	data = (float *)malloc(rows * cols * sizeof(float));
	
	// sacrifice memory w/ speed, by pre-allocating a temporary buffer
	temp_data = (float *)malloc(rows * cols * sizeof(float));

	bAllocated = true;
	
	// set every element to 0
	if(clear)
	{
		vDSP_vclr(data, 1, rows*cols);
	}
}

// pass in existing data
// non-destructive by default
// this WILL destroy the passed in data when object leaves scope if
// with copy is not true
Mat::Mat(size_t r, size_t c, float *existing_buffer, bool withCopy)
{
	data = temp_data = NULL;
	
	bUserData = false;
	rows = r;
	cols = c;
	
	// sacrifice memory w/ speed, by pre-allocating a temporary buffer
	temp_data = (float *)malloc(rows * cols * sizeof(float));
	
	if(withCopy)
	{
		data = (float *)malloc(rows * cols * sizeof(float));
		
		cblas_scopy(rows*cols, existing_buffer, 1, data, 1);
		//memcpy(data, existing_buffer, sizeof(float)*r*c);
	}
	else {
		// user gave us data, don't free it.
		bUserData = true;
		data = existing_buffer;
	}
	
	bAllocated = true;
}

// set every element to a value
Mat::Mat(size_t r, size_t c, float val)
{
	data = temp_data = NULL;
	
	bUserData = false;
	rows = r;
	cols = c;
	data = (float *)malloc(rows * cols * sizeof(float));
	// sacrifice memory w/ speed, by pre-allocating a temporary buffer
	temp_data = (float *)malloc(rows * cols * sizeof(float));
	
	bAllocated = true;
	
	// set every element to val
	vDSP_vfill(&val, data, 1, rows * cols);
	
}

// copy-constructor, called during:
//		pkm::Mat a = rhs;
//		pkm::Mat a(rhs);
Mat::Mat(const Mat &rhs)
{		
	if(rhs.data != NULL)
	{
		bUserData = false;
		
		rows = rhs.rows;
		cols = rhs.cols;
		
		data = (float *)malloc(rows * cols * sizeof(float));
		
		// sacrifice memory w/ speed, by pre-allocating a temporary buffer
		temp_data = (float *)malloc(rows * cols * sizeof(float));
		
		bAllocated = true;
		
		cblas_scopy(rows*cols, rhs.data, 1, data, 1);
		//memcpy(data, rhs.data, sizeof(float)*rows*cols);
	}
	else {
		rows = 0;
		cols = 0;
		data = NULL;
		temp_data = NULL;
		bUserData = false;
		bAllocated = false;
	}
}

Mat Mat::operator=(const Mat &rhs)
{
	if(bAllocated && !bUserData && data != NULL)
	{
		free(data);
	}
	if(bAllocated && temp_data != NULL) 
	{
		free(temp_data);	
	}
	
	if(rhs.data != NULL)
	{
		bUserData = false;
		rows = rhs.rows;
		cols = rhs.cols;
		
		data = (float *)malloc(rows * cols * sizeof(float));
		
		// sacrifice memory w/ speed, by pre-allocating a temporary buffer
		temp_data = (float *)malloc(rows * cols * sizeof(float));
		
		bAllocated = true;
		
		if (data != rhs.data) 
			cblas_scopy(rows*cols, rhs.data, 1, data, 1);
		//memcpy(data, rhs.data, sizeof(float)*rows*cols);
		
		return *this;
	}
	else {
		bUserData = false;
		rows = 0;
		cols = 0;
		data = NULL;
		temp_data = NULL;
		
		bAllocated = false;
		return *this;
	}			
}



/////////////////////////////////////////



/////////////////////////////////////////


Mat Mat::getTranspose()
{
#ifndef DEBUG			
	assert(data != NULL);
#endif	
	Mat transposedMatrix(cols, rows);
	
	if (rows == 1 || cols == 1) {
		cblas_scopy(rows*cols, data, 1, transposedMatrix.data, 1);
		//memcpy(transposedMatrix.data, data, sizeof(float)*rows*cols);
	}
	else {
		vDSP_mtrans(data, 1, transposedMatrix.data, 1, cols, rows);
	}
	
	return transposedMatrix;
}


// get a diagonalized version of the current vector (non-destructive)
Mat Mat::getDiag()
{
#ifndef DEBUG
	assert(data != NULL);
#endif	
	if(rows == 1 && cols > 1 || cols == 1 && rows > 1)
	{
		int diagonal_elements = MAX(rows,cols);
		
		// create a square matrix
		Mat diagonalMatrix(diagonal_elements,diagonal_elements, true);
		
		// set diagonal elements to the current vector in data
		for (int i = 0; i < diagonal_elements; i++) {
			diagonalMatrix.data[i*diagonal_elements+i] = data[i];
		}
		return diagonalMatrix;
	}
	else {
		printf("[ERROR]: Cannot diagonalize a matrix. Either rows or cols must be == 1.");
		Mat A;
		return A;
	}
}

Mat Mat::diag(Mat &A)
{
	if(A.rows == 1 && A.cols > 1 || A.cols == 1 && A.rows > 1)
	{
		int diagonal_elements = MAX(A.rows,A.cols);
		
		// create a square matrix
		Mat diagonalMatrix(diagonal_elements,diagonal_elements, true);
		
		// set diagonal elements to the current vector in data
		for (int i = 0; i < diagonal_elements; i++) {
			diagonalMatrix.data[i*diagonal_elements+i] = A.data[i];
		}
		return diagonalMatrix;
	}
	else {
		printf("[ERROR]: Cannot diagonalize a matrix. Either rows or cols must be == 1.");
		Mat A;
		return A;
	}
}

Mat Mat::identity(size_t dim)
{
	
	// create a square matrix
	Mat identityMatrix(dim,dim, true);
	
	// set diagonal elements to the current vector in data
	for (size_t i = 0; i < dim; i++) {
		identityMatrix.data[i*dim+i] = 1;
	}
	
	return identityMatrix;
}


// set every element to a random value between low and high
void Mat::setRand(float low, float high)
{
	float width = (high-low);
	float *ptr = data;
	for (int i = 0; i < rows*cols; i++) {
		*ptr = low + (float(::random())/float(RAND_MAX))*width;
		++ptr;
	}
}

// create a random matrix
Mat Mat::rand(size_t r, size_t c, float low, float high)
{
	Mat randomMatrix(r, c);
	randomMatrix.setRand(low, high);
	return randomMatrix;
}

Mat Mat::sum(bool across_rows)
{
	// sum across rows
	if(across_rows)
	{
		Mat result(1, cols);
		for (int i = 0; i < cols; i++) {
			vDSP_sve(data+i, cols, result.data+i, rows);
		}				
		return result;
	}
	// cols
	else
	{
		Mat result(rows, 1);
		for (int i = 0; i < rows; i++) {
			vDSP_sve(data+(i*cols), 1, result.data+i, cols);
		}
		return result;
	}
	
}

// normalize the values for each row-vector
void Mat::setNormalize(bool row_major)
{
	if (row_major) {
		for (int r = 0; r < rows; r++) {
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
	}
	// or for each column
	else {
		for (int c = 0; c < cols; c++) {
			float min, max;
			vDSP_minv(&(data[c]), cols, &min, rows);
			vDSP_maxv(&(data[c]), cols, &max, rows);
			float height = max-min;
			min = -min;
			vDSP_vsadd(&(data[c]), cols, &min, &(data[c]), cols, rows);
			if (height != 0) {
				vDSP_vsdiv(&(data[c]), cols, &height, &(data[c]), cols, rows);	
			}
		}
	}
}

void Mat::print(bool row_major)
{
	if(row_major)
	{
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				printf("%8.2f ", data[r*cols + c]);
			}
			printf("\n");
		}
		printf("\n");
	}
	else {
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				printf("%8.2f ", data[c*rows + r]);
			}
			printf("\n");
		}
		printf("\n");
	}
	
}



