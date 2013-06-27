// -----------------------------------------------------------------------------
//  pkmDTW.h
//  pkmMatrix
//
//  Created by Parag Mital on 10/19/12.
//  Copyright (c) 2012 Parag K Mital. All rights reserved.
//
// -----------------------------------------------------------------------------

#pragma once

#include "pkmMatrix.h"

#define WITH_OF
//#define WITH_FEATURE_WEIGHTING

#ifdef WITH_OF
#include "ofMain.h"
#endif

using namespace pkm;

// -----------------------------------------------------------------------------
class pkmDTW
{
public:
    // -------------------------------------------------------------------------
    pkmDTW()
    {
        bSetQuery = false;
        bHaveCandidates = false;
        bUseZNormalize = false;
        bUseCosineDistance = false;
        
        range = 0.5;
#ifdef WITH_FEATURE_WEIGHTING
        featureWeights = Mat(1,36,1.0f);
        featureWeights[0] = 20.0f;
#endif
        bestSoFar = INFINITY;
    }
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    //  Change the possible range of the warping envelope
    //
    //  'r' is a floating point value within (0, 1) 
    //  This value determines how much the warping is allowed to move from the diagonal
    // -------------------------------------------------------------------------
    void setRange(float r)
    {
        range = r;
    }
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    //  Add elements to the database of possible candidates
    //
    //  'candidate': size is frames x dimensions
    // -------------------------------------------------------------------------
    void addToDatabase(Mat &el)
    {
        candidates.push_back(el);
        
        bHaveCandidates = true;
    }
    // -------------------------------------------------------------------------
    
    
    // -------------------------------------------------------------------------
    void getNearestCandidate(Mat &q, 
                             float &distance, 
                             int &subscript, 
                             vector<int> &bestPathI,  // candidate's frame   (source)
                             vector<int> &bestPathJ)  // query's frame       (target)
    {
        if (!bHaveCandidates) {
            cout << "[ERROR::pkmDTW]: Add sequences to the database first using pkmDTW::addToDatabase(el)!" << endl;
            return;
        }
        // establish the query
        setQuery(q);
        subscript = 0;
        // search all candidates linearly
        for (int i = 0; i < candidates.size(); i++) 
        {
            vector<int> pathI, pathJ;    
            Mat differenceMatrix, dtwDistance;
            differenceMatrix = computeDifferenceMatrix(candidates[i]);
            float thisDistance = dtw(q, candidates[i], differenceMatrix, dtwDistance, pathI, pathJ);

            //cout << "i: " << i << endl;
            //dtwDistance.print();
            if (thisDistance < bestSoFar) 
            {
                bestSoFar = thisDistance;
                bestPathI = pathI;
                bestPathJ = pathJ;
                subscript = i;
            }
        }
        distance = bestSoFar;
        bestSoFar = INFINITY;
        
        
    }
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    void getNearestCandidateEuclidean(const Mat &q, 
                                      float &distance, 
                                      int &subscript) 
    {
        if (!bHaveCandidates) {
            cout << "[ERROR::pkmDTW]: Add sequences to the database first using pkmDTW::addToDatabase(el)!" << endl;
            return;
        }
        subscript = 0;
        Mat query = q;
        Mat distanceMatrix = Mat(q.rows, q.cols);
        // search all candidates linearly
        for (int i = 0; i < candidates.size(); i++) 
        {
            query.subtract(candidates[i], distanceMatrix);
            distanceMatrix.abs();
            Mat distance2 = distanceMatrix.sum(false);
            float thisDistance = Mat::sum(distance2);
            if (thisDistance < bestSoFar) 
            {
                bestSoFar = thisDistance;
                subscript = i;
            }
        }
        distance = bestSoFar;
        bestSoFar = INFINITY;
    }
    // -------------------------------------------------------------------------
  
    
    // -------------------------------------------------------------------------
    void save()
    {
        FILE *fp;
#ifdef WITH_OF
        fp = fopen(ofToDataPath("dtwDatabase/dtw.txt").c_str(), "w");
#else
        fp = fopen("dtwDatabase/dtw.txt", "w");
#endif
        fprintf(fp, "%d", candidates.size());
        fclose(fp);
        for (int i = 0; i < candidates.size(); i++) {
            char buf[256];
            sprintf(buf, "dtwDatabase/candidate%08d.txt", i);
#ifdef WITH_OF
            candidates[i].save(ofToDataPath(buf).c_str());
#else
            candidates[i].save(buf);
#endif
        }
    }
    // -------------------------------------------------------------------------
    
    
    // -------------------------------------------------------------------------
    void load()
    {
        int numCamdidates = 0;
        FILE *fp;
#ifdef WITH_OF
        fp = fopen(ofToDataPath("dtwDatabase/dtw.txt").c_str(), "r");
#else
        fp = fopen("dtwDatabase/dtw.txt", "r");
#endif
        fscanf(fp, "%d", &numCamdidates);
        fclose(fp);
        candidates.resize(numCamdidates);
        for (int i = 0; i < numCamdidates; i++) {
            char buf[256];
            sprintf(buf, "dtwDatabase/candidate%08d.txt", i);
#ifdef WITH_OF
            candidates[i].load(ofToDataPath(buf).c_str());
#else
            candidates[i].load(buf);
#endif
        }
        
        if (numCamdidates > 0) {
            bHaveCandidates = true;
        }
    }
    // -------------------------------------------------------------------------
    
protected:
    
    // -------------------------------------------------------------------------
    // Establish the query to compare against all candidates
    //
    //  'q' is a query matrix of size T x D,
    //  T is the number of frames or time-steps
    //  D is the dimension of each data point
    // -------------------------------------------------------------------------
    void setQuery(Mat &q)
    {
        query = q;
        if(bUseZNormalize)
        {
            query.zNormalizeEachCol();
        }
        queryTransposed = query;
        queryTransposed.setTranspose();
        
        Mat temp = query;
        temp.sqr();
        queryNormalization = temp.sum(false);
        queryNormalization.sqrt();
        queryNormalization.setTranspose();
        
        bSetQuery = true;
    }
    // -------------------------------------------------------------------------
    
    
    // -------------------------------------------------------------------------
    // Compare stored query with the incoming matrix candidate
    // 
    //  'candidate': size is frames x dimensions
    //  'differenceMatrix': size will be candidate's rows x query's rows,
    //      i.e. differenceMatrix(i,j) indexes [1 - cosine distance of (candidate(i,:), query(j,:))]
    // -------------------------------------------------------------------------
    Mat computeDifferenceMatrix(Mat &candidate)
    {
        Mat differenceMatrix;
        if (bSetQuery) {
            if(bUseCosineDistance)
            {
                Mat temp(candidate.rows, candidate.cols);
                temp.copy(candidate);
                if (bUseZNormalize) {
                    temp.zNormalizeEachCol();
                }
                temp.sqr();
                Mat candidateNormalization = temp.sum(false);
                candidateNormalization.sqrt();
                Mat normalization = candidateNormalization.GEMM(queryNormalization);
                differenceMatrix = candidate.GEMM(queryTransposed);
                differenceMatrix.divide(normalization);
                
                // remove these next 3 lines for a similarity matrix instead
                float factor = -1;
                float term = 1;
                vDSP_vsmsa(differenceMatrix.data, 1, &factor, &term, differenceMatrix.data, 1, differenceMatrix.size());
            }
            else
            {
                int padding = query.rows * range;
                differenceMatrix = Mat(candidate.rows, query.rows, HUGE_VALF);
                
                Mat ssd(1, candidate.cols);
                for (int i = 0; i < candidate.rows; i++) 
                {
                    Mat p1(1, candidate.cols, candidate.row(i), false);
                    for (int j = max(0, i - padding); j < min(query.rows, i + padding - 1); j++) 
                    {
                        Mat p2(1, query.cols, query.row(j), false);
                        p1.subtract(p2, ssd);
                        ssd.sqr();
#ifdef WITH_FEATURE_WEIGHTING
                        ssd.multiply(featureWeights, ssd);
#endif
                        differenceMatrix.data[query.rows*i + j] = sqrtf(Mat::sum(ssd));
                        
                        /*
                        float sum = 0;
                        for (int k = 0; k < candidate.cols; k++) {
                            //sum += sqrtf(powf(candidate.row(i)[k] - query.row(j)[k],2));
                            sum += fabs(candidate.row(i)[k] - query.row(j)[k]);                            
                        }
                        //sum = (sum / (float)candidate.cols);
                        differenceMatrix.data[query.rows*i + j] = sum;
                         */
                    }
                }
                //differenceMatrix.print();
            }
        }
        return differenceMatrix;
    }
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    float dtw(const Mat &query, 
             const Mat &candidate, 
             Mat &differenceMatrix,
             Mat &dtwDistance,
             vector<int> &pathI,
             vector<int> &pathJ)
    {
        // calculate the dtw distance matrix
        Mat traceBack(differenceMatrix.rows, differenceMatrix.cols);
        dtwDistance = differenceMatrix;
        int subscriptRange = differenceMatrix.cols * range;
        float x, y, z;
        int i, j;
        for (i = 0; i < differenceMatrix.rows; i++) 
        {
            float *dist = dtwDistance.row(i);
            float *tb = traceBack.row(i);
            float minCost = INFINITY;
            //int k = max(0, subscriptRange - i);
            //for (j = max(0, i - subscriptRange); j < min(i + subscriptRange - 1, differenceMatrix.cols); j++, k++) 
            for (j = 0; j < differenceMatrix.cols; j++) 
            {
                if (i == 0 && j == 0) {
                    *dist = *(dtwDistance.data);
                    minCost = *dist;
                    continue;
                }
                
                /*
                // get distance for all branches
                if ((j - 1 < 0) || (k - 1 < 0))                     x = INFINITY;                     // horizontal
                else                                                x = dtwDistance.row(i)[j-1]; 
                if ((i - 1 < 0 )|| (k + 1 > 2 * subscriptRange))    y = INFINITY;                     // veritcal
                else                                                y = dtwDistance.row(i-1)[j];      
                if ((i - 1 < 0) || (j - 1 < 0))                     z = INFINITY;                     // diagonal
                else                                                z = dtwDistance.row(i-1)[j-1];
                */
                
                // get distance for all branches
                if (j - 1 < 0)                x = INFINITY;                     // horizontal
                else                          x = dtwDistance.row(i)[j-1]; 
                if (i - 1 < 0)                y = INFINITY;                     // veritcal
                else                          y = dtwDistance.row(i-1)[j];      
                if (i - 1 < 0 || j - 1 < 0)   z = INFINITY;                     // diagonal
                else                          z = dtwDistance.row(i-1)[j-1];
                
                
                // find minimum branch and store path
                float val;
                if (x < y) {        // horizontal
                    val = x;
                    tb[j] = 0;
                }
                else {              // vertical
                    val = y;
                    tb[j] = 1;
                }
                if (z < val) {      // diagonal
                    val = z;
                    tb[j] = 2;
                }
                
                // aggregate distance
                dist[j] = val + dist[j];
                
                if (dist[j] < minCost) {
                    minCost = dist[j];
                }
            }
            
            // abandon early
            if (minCost > bestSoFar) {
                return INFINITY;
            }
        }
        
        // calculate path
        i--;
        j--;
        while(i >= 0 && j >= 0) 
        {
            pathI.push_back(i);
            pathJ.push_back(j);
            float t = traceBack.row(i)[j];
            if (t == 0) {                   // horizontal
                j--;
            }
            else if (t == 1) {              // vertical
                i--;
            }
            else {                          // diagonal
                i--;
                j--;
            }
        }
        return *(dtwDistance.last());
    }
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Calculates the Sakoe-Chiba Band for a multidimensional input T x D
    //
    //  'input' is a T x D matrix (untouched, though not declared const since
    //  vDSP's max/min functions are not const) with:
    // 
    //  T frames, or time-steps
    //  D dimensions
    //
    //  [ d1 d2 .  . dn ]  t1
    //  [ .  .  .  .  . ]  t2
    //  [ .  .  .  .  . ]   .  time
    //  [ .  .  .  .  . ]   .
    //  [ .  .  .  .  . ]  tm
    //         dims
    //
    // 'upperBound' is computed as: UW_i = max(C_max(1,i-r), . . . , C_min(i+r,n)) and
    // 'lowerBound' is computed as: LW_i = min(C_max(1,i-r), . . . , C_min(i+r,n))
    // -------------------------------------------------------------------------
    void calculateBounds(Mat &input, 
                         Mat &upperBound, 
                         Mat &lowerBound);

    
    
    
private:
    // -------------------------------------------------------------------------
    float           bestSoFar;
    float           range;
    Mat             query, queryTransposed, queryNormalization;
    vector<Mat>     candidates;
    Mat             queryLB;
    Mat             queryUB;
#ifdef WITH_FEATURE_WEIGHTING
    Mat             featureWeights;
#endif
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    bool            bUseZNormalize, bSetQuery, bHaveCandidates, bUseCosineDistance;
    // -------------------------------------------------------------------------
};