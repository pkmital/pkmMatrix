//
//  pkmDTW.cpp
//  pkmMatrix
//
//  Created by Parag Mital on 10/19/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "pkmDTW.h"


// -------------------------------------------------------------------------
void pkmDTW::calculateBounds(Mat &input, Mat &upperBound, Mat &lowerBound)
{
    // compute how much we allow the time to stretch in terms of the query's subscripts
    // Sakoe-Chiba uses fixed range
    int subscriptRange = range * input.rows;
    
    upperBound = Mat(input.rows, input.cols);
    lowerBound = Mat(input.rows, input.cols);
    
    for (int i = 0; i < input.rows; i++) {
        int startRow = max(0, i - subscriptRange);
        int endRow = min(i + subscriptRange, input.rows - 1);
        for (int j = 0; j < input.cols; j++) {
            vDSP_Length len = endRow - startRow + 1;
            vDSP_maxv(input.row(startRow) + j, input.cols, upperBound.row(i) + j, len);
            vDSP_minv(input.row(startRow) + j, input.cols, lowerBound.row(i) + j, len);
        }
    }
}
// -------------------------------------------------------------------------
