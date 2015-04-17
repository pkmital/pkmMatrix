// -----------------------------------------------------------------------------
//  pkmDTW.h
//  pkmMatrix
//
//  Created by Parag Mital on 10/19/12.
//  Copyright (c) 2012 Parag K Mital. All rights reserved.
//
/*
Copyright (C) 2011 Parag K. Mital

The Software is and remains the property of Parag K Mital
("pkmital") The Licensee will ensure that the Copyright Notice set
out above appears prominently wherever the Software is used.

The Software is distributed under this Licence:

- on a non-exclusive basis,

- solely for non-commercial use in the hope that it will be useful,

- "AS-IS" and in order for the benefit of its educational and research
purposes, pkmital makes clear that no condition is made or to be
implied, nor is any representation or warranty given or to be
implied, as to (i) the quality, accuracy or reliability of the
Software; (ii) the suitability of the Software for any particular
use or for use under any specific conditions; and (iii) whether use
of the Software will infringe third-party rights.

pkmital disclaims:

- all responsibility for the use which is made of the Software; and

- any liability for the outcomes arising from using the Software.

The Licensee may make public, results or data obtained from, dependent
on or arising out of the use of the Software provided that any such
publication includes a prominent statement identifying the Software as
the source of the results or the data, including the Copyright Notice
and stating that the Software has been made available for use by the
Licensee under licence from pkmital and the Licensee provides a copy of
any such publication to pkmital.

The Licensee agrees to indemnify pkmital and hold them
harmless from and against any and all claims, damages and liabilities
asserted by third parties (including claims for negligence) which
arise directly or indirectly from the use of the Software or any
derivative of it or the sale of any products based on the
Software. The Licensee undertakes to make no liability claim against
any employee, student, agent or appointee of pkmital, in connection
with this Licence or the Software.


No part of the Software may be reproduced, modified, transmitted or
transferred in any form or by any means, electronic or mechanical,
without the express permission of pkmital. pkmital's permission is not
required if the said reproduction, modification, transmission or
transference is done without financial return, the conditions of this
Licence are imposed upon the receiver of the product, and all original
and amended source code is included in any transmitted product. You
may be held legally responsible for any copyright infringement that is
caused or encouraged by your failure to abide by these terms and
conditions.

You are not permitted under this Licence to use this Software
commercially. Use for which any financial return is received shall be
defined as commercial use, and includes (1) integration of all or part
of the source code or the Software into a product for sale or license
by or on behalf of Licensee to third parties or (2) use of the
Software or any derivative of it for research with the final aim of
developing software products for sale or license to a third party or
(3) use of the Software or any derivative of it for research with the
final aim of developing non-software products for sale or license to a
third party, or (4) use of the Software to provide any service to an
external organisation for which payment is received. If you are
interested in using the Software commercially, please contact pkmital to
negotiate a licence. Contact details are: parag@pkmital.com
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