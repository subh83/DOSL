/** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal search Library (DOSL)                                              *
*    A template-based C++ library for discrete search                                    *
*    Version 3.x                                                                         *
*    ----------------------------------------------------------                          *
*    Copyright (C) 2017  Subhrajit Bhattacharya                                          *
*                                                                                        *
*    This program is free software: you can redistribute it and/or modify                *
*    it under the terms of the GNU General Public License as published by                *
*    the Free Software Foundation, either version 3 of the License, or                   *
*    (at your option) any later version.                                                 *
*                                                                                        *
*    This program is distributed in the hope that it will be useful,                     *
*    but WITHOUT ANY WARRANTY; without even the implied warranty of                      *
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                       *
*    GNU General Public License for more details <http://www.gnu.org/licenses/>.         *
*                                                                                        *
*                                                                                        *
*    Contact:  subhrajit@gmail.com                                                       *
*              https://www.lehigh.edu/~sub216/ , http://subhrajit.net/                   *
*                                                                                        *
*                                                                                        *
*************************************************************************************** **/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>

// Open CV:
#include <opencv2/opencv.hpp> 
#include <opencv2/highgui.hpp> 

// DOSL
#include <dosl/encapsulations/cvMulticlassPathPlanner.tcc>
#include <dosl/aux-utils/string_utils.hpp> // compute_program_path

#ifndef _DOSL_ALGORITHM // can pass at command line during compilation: -D_DOSL_ALGORITHM=AStar
    #define _DOSL_ALGORITHM  AStar
#endif


int main(int argc, char *argv[])
{
    // Call syntax 1: a.out nPaths
    // Call syntax 2:  a.out "startx,starty" "goalx,goaly" nPaths
    // Call syntax 3:  a.out map_file.png "startx,starty" "goalx,goaly" nPaths
    
    compute_program_path();
    
    int nPaths = 3;
    std::string map_file = program_path + "../files/expt/L457.png";
    cv::Point start (60,60);
    cv::Point goal (200,300);
    
    if (argc>1) {
        nPaths = atoi(argv[argc-1]); // last parameter
        if (argc>2) {
            sscanf (argv[argc-2], "%d,%d", &goal.x, &goal.y); // second-to-last parameter
            sscanf (argv[argc-3], "%d,%d", &start.x, &start.y); // third-to-last parameter
            if (argc>4)
                map_file = argv[1];
        }
    }
    
    printf (_BOLD _YELLOW "Note: " YELLOW_ BOLD_ "Using algorithm " _YELLOW  MAKESTR(_DOSL_ALGORITHM)  YELLOW_ 
                        ". Run 'make' to recompile with a different algorithm.\n"
                        "Note: This program currently does not support 'ThetaStar' algorithm.\n");
    
    // ------------------------------------------------------------------
    
    cv::Mat obs_map = cv::imread (map_file, CV_LOAD_IMAGE_GRAYSCALE);
    
    // find paths
    cvMulticlassPathPlanner<_DOSL_ALGORITHM> path_planner (obs_map, start, goal, nPaths, true);
    path_planner.find_paths ();
    
    // draw paths and display
    obs_map = path_planner.draw_paths ({}, 2);
    cv::imshow ("Final Paths", obs_map);
    cv::waitKey(0);
}

