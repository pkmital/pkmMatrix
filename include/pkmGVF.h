// -----------------------------------------------------------------------------
//  pkmDTW.h
//  pkmMatrix
//
//  Created by Parag Mital on 10/19/12.
//  Copyright (c) 2012 Parag K Mital. All rights reserved.
//
/*
Copyright (C) 2011 Parag K. Mital

 This program is free software: you can redistribute it and/or modify  
 it under the terms of the GNU General Public License as published by  
 the Free Software Foundation, version 3.
 
 This program is distributed in the hope that it will be useful, but 
 WITHOUT ANY WARRANTY; without even the implied warranty of 
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License 
 along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

// -----------------------------------------------------------------------------

#pragma once

#include "pkmMatrix.h"
#include "GestureVariationFollower.h"
#include <Eigen/Core>

#define WITH_OF

#ifdef WITH_OF
#include "ofMain.h"
#endif

using namespace pkm;

// -----------------------------------------------------------------------------
class pkmGVF
{
public:
    // -------------------------------------------------------------------------
    pkmGVF()
    {
        numCandidates = 0;
        bHaveCandidates = false;
        
        Eigen::VectorXf sigs(3);
        sigs(0) = 0.00001;
        sigs(1) = 0.001;
        sigs(2) = 0.00001;
        int numParticles = 1e6;
        float tolerance = 1.0;
        gvf = new GestureVariationFollower(numParticles, sigs, 1.0 / (tolerance*tolerance), numParticles / 100);

    }
    // -------------------------------------------------------------------------

    
    // -------------------------------------------------------------------------
    //  Add elements to the database of possible candidates
    //
    //  'candidate': size is frames x dimensions
    // -------------------------------------------------------------------------
    void addToDatabase(Mat &el)
    {
        vector<float> lut_el;
        lut_el.push_back(allFeatures.rows);
        allFeatures.push_back(el);
        lut_el.push_back(allFeatures.rows - lut_el[0]);
        lut.push_back(lut_el);
        
        numCandidates++;
        bHaveCandidates = true;
    }
    // -------------------------------------------------------------------------
    
    
    // -------------------------------------------------------------------------
    void normalizeDatabase()
    {
        meanFeature = allFeatures.mean();
        stdFeature = allFeatures.stddev();
        
        allFeatures.zNormalizeEachCol();
        
        allFeatures.save(ofToDataPath("all-features-normalized.txt"));
    }
    // -------------------------------------------------------------------------
    
    
    void normalizeQuery(Mat &query)
    {
        query.subtract(meanFeature);
        query.divide(stdFeature);
    }
    
    // -------------------------------------------------------------------------
    void buildDatabase()
    {
        
        for(int i = 0; i < numCandidates; i++)
        {
            gvf->addTemplate();
            
            int idxOfGesture = lut.row(i)[0];
            int numFramesInGesture = lut.row(i)[1];
            
            vector<float> query;
            query.resize(allFeatures.cols);
            
            for(int f = 0; f < numFramesInGesture; f++)
            {
                std::memcpy(&query[0], allFeatures.row(idxOfGesture + f), sizeof(float)*allFeatures.cols);
//
//                cout << "gesture: " << i;
//                for(int m = 0; m < allFeatures.cols; m++)
//                    cout << " " << query[m];
//                cout << endl;
                
                gvf->fillTemplate(i, query);
            }
        }
        
        Eigen::VectorXf meanPVRS(3);
        Eigen::VectorXf rangePVRS(3);
        
        meanPVRS << 0.0 , 1.0 , 1.0; // phase, speed, scale
        rangePVRS << 0.1 , 0.1 , 0.1; // how much to spread uniformly
        
        gvf->spreadParticles(meanPVRS, rangePVRS);
        
        gvf->saveTemplates(ofToDataPath("gvf-database.txt"));
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
            cout << "[ERROR::pkmGVF]: Add sequences to the database first using pkmDTW::addToDatabase(el)!" << endl;
            return;
        }
        
        normalizeQuery(q);
        
        vector<float> query;
        query.resize(q.cols);
        std::memcpy(&query[0], q.data, sizeof(float)*q.cols);
        
        
//        Eigen::VectorXf meanPVRS(3);
//        Eigen::VectorXf rangePVRS(3);
//        
//        meanPVRS << 0.0 , 1.0 , 1.0; // phase, speed, scale
//        rangePVRS << 0.1 , 0.1 , 0.1; // how much to spread uniformly
//        
//        gvf->spreadParticles(meanPVRS, rangePVRS);

        gvf->infer(query);
        int c;
        float i;
        gvf->getEstimatedStatus(c, i);
        subscript = c;
        cout << "gvf | " << c << " " << i << endl;
        bestPathI.push_back(i);
    }
    // -------------------------------------------------------------------------
    
    
    // -------------------------------------------------------------------------
    void getNearestCandidate(float *q, int numFeatures,
                             float &distance,
                             int &subscript,
                             vector<int> &bestPathI,  // candidate's frame   (source)
                             vector<int> &bestPathJ)  // query's frame       (target)
    {
        if (!bHaveCandidates) {
            cout << "[ERROR::pkmGVF]: Add sequences to the database first using pkmDTW::addToDatabase(el)!" << endl;
            return;
        }
        
        Mat qmat(1,numFeatures,q,false);
        normalizeQuery(qmat);
        
        vector<float> query;
        query.resize(numFeatures);
        std::memcpy(&query[0], q, sizeof(float)*numFeatures);
        
        
        //        Eigen::VectorXf meanPVRS(3);
        //        Eigen::VectorXf rangePVRS(3);
        //
        //        meanPVRS << 0.0 , 1.0 , 1.0; // phase, speed, scale
        //        rangePVRS << 0.1 , 0.1 , 0.1; // how much to spread uniformly
        //
        //        gvf->spreadParticles(meanPVRS, rangePVRS);
        
        gvf->infer(query);
        float i;
        gvf->getEstimatedStatus(subscript, i);
        cout << "gvf | " << subscript << " " << roundf(i*lut.row(subscript)[1]) << "/" << lut.row(subscript)[1] << endl;
        bestPathI.push_back(roundf(i*lut.row(subscript)[1]));
    }
    //
    
    // -------------------------------------------------------------------------
    void save()
    {
        allFeatures.save(ofToDataPath("all-features.txt"));
        lut.save(ofToDataPath("all-features-lut.txt"));
    }
    // -------------------------------------------------------------------------
    
    
    // -------------------------------------------------------------------------
    void load()
    {
        allFeatures.load(ofToDataPath("all-features.txt"));
        lut.load(ofToDataPath("all-features-lut.txt"));
        
        numCandidates = lut.rows;
        
        normalizeDatabase();
        
        buildDatabase();
        
        if (numCandidates > 0) {
            bHaveCandidates = true;
        }
        
    }
    // -------------------------------------------------------------------------
    
protected:
    

private:
    
    Mat lut;                // 0, what index in allFeatures is the gesture
                            // 1, how many rows in allFeatures the gesture has
    Mat allFeatures;
    
    Mat meanFeature, stdFeature;
        
    GestureVariationFollower *gvf;
    int numCandidates;
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    bool            bHaveCandidates;
    // -------------------------------------------------------------------------
};