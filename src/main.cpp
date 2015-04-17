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
#include <vector>

using namespace pkm;
using namespace std;


int main (int argc, char * const argv[]) {
    
    size_t n_observations = 10000;
    size_t n_features = 500;
    
    pkm::Mat data(n_observations, n_features);
//    data.printAbbrev();
    
    for (int i = 0; i < 100; i++)
    {
        auto start = std::chrono::steady_clock::now();
        pkm::Mat mean_data = data.mean();
        auto end = std::chrono::steady_clock::now();
        std::cout << "Mean calculated in " << double((end-start).count())/double(std::chrono::steady_clock::period::den) << "s" << std::endl;
    }

    
	return 0;
}
