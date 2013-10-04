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
#include "pkmDTW.h"
#include <vector>

using namespace pkm;
using namespace std;

int main (int argc, char * const argv[]) {
    
    //srand ( time(NULL) );
    srand ( 0 );
    
    int rows = 1;
    int cols = 5;
    
    Mat a(rows, cols);
    for (int i = 0; i < rows * cols; i++) {
        a[i] = i * 10;
    }
    a.print();
    
    Mat c = pkm::Mat::resize(a, rows * cols * 10);
    c.print();
    
    
//    pkmDTW dtw;
//    unsigned long t = clock();
//    int numCandidates = 5;
//    for (int i =0; i < numCandidates; i++) {
//        //Mat e(rows + rand() % (rows * 2), cols);
//        Mat e(rows, cols);
//        e.setRand(-5,5);
//        dtw.addToDatabase(e);
//        //e.print();
//    }
//    t = clock() - t;
//    cout << "ticks: " << t << " (" << (float)t / CLOCKS_PER_SEC << " seconds)" << endl;
//    float distance;
//    int subscript;
//    vector<int> bestPathI, bestPathJ;
//    dtw.getNearestCandidate(d, distance, subscript, bestPathI, bestPathJ);
//    t = clock() - t;
//    cout << "best: " << subscript << " dist:  " << distance << endl;
//    cout << "ticks: " << t << " (" << (float)t / CLOCKS_PER_SEC << " seconds)" << endl;
//    dtw.getNearestCandidateEuclidean(d, distance, subscript);
//    t = clock() - t;
//    cout << "best: " << subscript << " dist:  " << distance << endl;
//    cout << "ticks: " << t << " (" << (float)t / CLOCKS_PER_SEC << " seconds)" << endl;
//    
    /*
    Mat b = pkm::Mat::sgn(a);
    b.print();
    
    float data[] = {1,1,7,1,2,1,1,1,3};
    Mat c(3,3,data,false);
    c.print();
    c.inv();
    c.print();
    */
    
    
//    unsigned long t = clock();
//    for (int i = 0; i < 100000; i++) {
//        float f = sqrtf(i);
//    }
//    t = clock() - t;
//    cout << "ticks: " << t << endl;
//    
//    Mat a(1,100000,false);
//    
//    for (int i = 0; i < 100000; i++) {
//        a[i] = i;
//    }
//    a.sqrt();
//    t = clock() - t;
//    cout << "ticks: " << t << endl;
    
    
//    int i, j , c1, c2, pivot[size], ok;
//    float A[size][size], b[size], AT[size*size];	/* single precision!!! */
//    
//    
//    A[0][0]=3.1;  A[0][1]=1.3;  A[0][2]=-5.7;	/* matrix A */
//    A[1][0]=1.0;  A[1][1]=-6.9; A[1][2]=5.8;	
//    A[2][0]=3.4;  A[2][1]=7.2;  A[2][2]=-8.8;	
//    
//    b[0]=-1.3;			/* if you define b as a matrix then you */
//    b[1]=-0.1;			/* can solve multiple equations with */
//    b[2]=1.8;			/* the same A but different b */ 	
//    
//    for (i=0; i<size; i++)		/* to call a Fortran routine from C we */
//    {				/* have to transform the matrix */
//        for(j=0; j<size; j++) AT[j+size*i]=A[j][i];		
//    }						
//    
//    c1=size;			/* and put all numbers we want to pass */
//    c2=1;    			/* to the routine in variables */
//    
//    /* find solution using LAPACK routine SGESV, all the arguments have to */
//    /* be pointers and you have to add an underscore to the routine name */
//    sgesv_(&c1, &c2, AT, &c1, pivot, b, &c1, &ok);      
//    
//    /*
//     parameters in the order as they appear in the function call
//     order of matrix A, number of right hand sides (b), matrix A,
//     leading dimension of A, array that records pivoting, 
//     result vector b on entry, x on exit, leading dimension of b
//     return value */ 
//    
//    for (j=0; j<size; j++) printf("%e\n", b[j]);	/* print vector x */


    /*
    Mat a(1,8,5.0f);
    a.print();
    
    Mat b;
    b = a.getDiag();
    b.print();
    
    b.removeRow(3);
    b.print();
    b.removeRow(6);
    b.print();
    b.removeRow(1);
    b.print();
    
    Mat c = Mat::sqrt(b);
    c.print();
    */
     
     
    
    
    /*
    float means[] = {0.5, 0.5};
    float sigma[] = {0.15915, 0.0, 0.0, 0.15915};
    float input[] = {0.5, 0.5};
    
    Mat m(1,2,means,true);
    m.print();
    Mat s(2,2,sigma,true);
    s.print();
    Mat i(1,2,input,true);
    i.print();
    
    float size = 10;
    for (int x = 0; x < size; x++) {
        i[0] = x / size;
        for (int y = 0; y < size; y++) {
            i[1] = y / size;
            cout << Mat::gaussianPosterior(i, m, s) << " ";
        }
        cout << endl;
    }
    */
    
    
	/*
	Mat C = A[B];
	C.print();
	
	C.setRand(10.0, 20.0);
	A.copy(C, B);
	A.print();
	
	
	if(0)
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
	*/
	return 0;
}
