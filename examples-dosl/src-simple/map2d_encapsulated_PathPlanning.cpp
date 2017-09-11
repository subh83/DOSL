#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>

// Open CV:
#include <opencv2/opencv.hpp> 
#include <opencv2/highgui.hpp> 

// DOSL
#include <dosl/encapsulations/cvPathPlanner.tcc>
#include <dosl/aux-utils/string_utils.hpp> // compute_program_path

#ifndef _DOSL_ALGORITHM // can pass at command line during compilation: -D_DOSL_ALGORITHM=AStar
    #define _DOSL_ALGORITHM  AStar
#endif


int main(int argc, char *argv[])
{
    // a.out map_file.png "startx,starty" "goalx,goaly"
    // a.out "startx,starty" "goalx,goaly"
    
    compute_program_path();
    
    std::string map_file = program_path + "../files-expt/L457.png";
    cv::Point start (60,60);
    cv::Point goal (200,300);
    
    if (argc>1) {
        sscanf (argv[argc-1], "%d,%d", &goal.x, &goal.y); // last parameter
        sscanf (argv[argc-2], "%d,%d", &start.x, &start.y); // second-to-last parameter
        if (argc>3)
            map_file = argv[1];
    }
    
    printf (_BOLD _YELLOW "Note: " YELLOW_ BOLD_ "Using algorithm " _YELLOW  MAKESTR(_DOSL_ALGORITHM)  YELLOW_ 
                        ". Run 'make' to recompile with a different algorithm.\n");
    
    // ------------------------------------------------------------------
    
    cv::Mat obs_map = cv::imread (map_file, CV_LOAD_IMAGE_GRAYSCALE);
    
    // find paths
    cvPathPlanner<_DOSL_ALGORITHM> path_planner (obs_map, start, goal, false);
    
    // draw paths and display
    obs_map = path_planner.draw_path (cvScalar(0.0,0.0,255.0), 2);
    cv::imshow ("Final Path", obs_map);
    cv::waitKey(0);
}

