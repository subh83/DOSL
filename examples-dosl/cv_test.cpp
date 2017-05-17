/** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal Search Library (DOSL)                                              *
*    A template-based C++ library for discrete search                                    *
*    Version 3.x                                                                         *
*    ----------------------------------------------------------                          *
*    Copyright (C) 2016  Subhrajit Bhattacharya                                          *
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

// Other libraries:
// Open CV:
#include <opencv/cv.h>
#include <opencv/cvaux.h>
#include <opencv/highgui.h>
#include "cvParseMap2d.h"

int main(int argc, char *argv[])
{
    std::string program_fName (argv[0]);
    std::string program_folderName = program_fName.substr(0, program_fName.find_last_of("/\\")+1);
    
    std::string imagefName = program_folderName + "exptfiles/L457.png";
    cvParseMap2d my_map = cvParseMap2d (imagefName, true); // setting 'true' for the second parameter computes representative points
    
    // Make changes to the map
    // -----------------------
    
    for (int a=2; a<10; ++a)
        for (int b=2; b<15; ++b)
            my_map.getPixel(a, b) = OBSTACLE_PIXEL;
    
    my_map.update(); // need to update since the map was changed (will recompute representative points)
    
    // Plotting:
    // --------
    my_map.print_info ();
    
    // Get the un-occupied map
    cv::Mat image_to_display = my_map.getCvMat (COLOR_MAP);
    
    // plot representative points
    for (int i=0; i<my_map.repPts.size(); ++i) {
        cv::circle (image_to_display, my_map.repPts[i], 3, cv::Scalar(0,0,255), -1);
        image_to_display.at<cv::Vec3b> (my_map.repPts[i].y, my_map.repPts[i].x) = cv::Vec3b(0,255,0);
        printf("reprsentative point %d: (%d, %d).\n", i, my_map.repPts[i].x, my_map.repPts[i].y);
    }
    
    // display
    printf("Map: %s. width = %d, high = %d.\n", imagefName.c_str(), my_map.width(), my_map.height());
    cv::resize (image_to_display, image_to_display, cv::Size(), 2, 2);
    cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE);
    cv::imshow("Display window", image_to_display);
    cv::waitKey(0);
}


