/** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal Search Library (DOSL)                                              *
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


#ifndef __DOSL_PARSEMAP2D_H
#define __DOSL_PARSEMAP2D_H

#include <vector>
#include <cstring>
#include <string>
#include <ostream>
#include <fstream>

// OpenCV
#include <opencv2/opencv.hpp> 
#include <opencv2/highgui.hpp> 

// debugging macro
#define _map_info(m)  printf(#m ": width = %d, high = %d.\n", m.cols, m.rows); std::cout << std::flush;


enum PARSE_MAP_TYPE { ORIGINAL_MAP, FREE_MAP, OBSTACLE_MAP, OBSTACLE_LABEL_MAP };
enum MAP_COLOR_TYPE { GRAYSCALE_MAP, COLOR_MAP };
enum MAP_PIXEL_VALS { OBSTACLE_PIXEL=0, DEFAULT_THRESHOLD=100, FREE_PIXEL=255 };

class cvParseMap2d {
public:
    cv::Mat map;
    int obsThreshold, repPtObsSizeThresh;
    bool repPtsComputed;
    
    // void update(void) -->
    void update (const cv::Mat& new_map=cv::Mat(), int computeRepPts=-1, int obsSizeThresh=-1, int threshold=-1) { 
        if (new_map.empty()) { /* use current 'map' */ }
        else if (new_map.channels() == 1)  new_map.copyTo (map);
        else  cv::cvtColor (new_map, map, cv::COLOR_RGB2GRAY);
        
        if (obsSizeThresh>=0)  repPtObsSizeThresh = obsSizeThresh;
        if (threshold>0)  obsThreshold = threshold;
        repPtsComputed = false;
        if ( computeRepPts>0 || (computeRepPts<0 && repPtsComputed>0) ) {
            computeRepresentativePoints();
            repPtsComputed = 1;
        }
    }
    
    void update (std::string fName, int computeRepPts=0, int obsSizeThresh=0, int threshold=DEFAULT_THRESHOLD) {
        // std::cout << "Reading image file: " << fName << std::endl;
        map = cv::imread (fName, CV_LOAD_IMAGE_GRAYSCALE);
        update (cv::Mat(), computeRepPts, obsSizeThresh, threshold);
    }
    
    // --------------------------------------
    
    cvParseMap2d() : repPtsComputed(false) { }
    cvParseMap2d(const cv::Mat& new_map, int computeRepPts=0, int obsSizeThresh=0, int threshold=DEFAULT_THRESHOLD) 
                                                    { update (new_map, computeRepPts, obsSizeThresh, threshold); }
    cvParseMap2d(std::string fName, int computeRepPts=0, int obsSizeThresh=0, int threshold=DEFAULT_THRESHOLD) 
                                                    { update (fName, computeRepPts, obsSizeThresh, threshold); }
    
    // --------------------------------------
    
    int width(void) { return (map.cols); }
    int height(void) { return (map.rows); }
    void print_info() {
        printf("cvParseMap2d (%x): width = %d, height = %d, channels = %d;\n"
               "\tobsThreshold = %d, repPtsComputed = %s, #repPts = %d.\n", 
                    this, map.cols, map.rows, map.channels(), 
                    obsThreshold, (repPtsComputed?"true":"false"), repPts.size());
       std::cout << std::flush;
   }; 
    
    uchar& getPixel (int x, int y) { return (map.at<uchar>(y,x)); }
    bool isFree (int x, int y) { return ((bool)(map.at<uchar>(y,x)>=obsThreshold)); }
    bool isInFrame (int x, int y) { return (x>=0 && y>=0 && x<map.cols && y<map.rows); }
    bool isObstacle (int x, int y) { return ((bool)(map.at<uchar>(y,x)<obsThreshold)); }
    bool isOutsideFrame (int x, int y) { return (x<0 || y<0 || x>=map.cols || y>=map.rows); }
    
    // --------------------------------------
    
    cv::Mat getCvMat (MAP_COLOR_TYPE color_type = GRAYSCALE_MAP, PARSE_MAP_TYPE map_type = FREE_MAP, int draw_rep_pts=1) {
        cv::Mat retMat; // return a copy
        
        if (map_type==FREE_MAP && !map.empty())
            cv::threshold (map, retMat, obsThreshold, FREE_PIXEL, cv::THRESH_BINARY);
        else if (map_type==OBSTACLE_MAP && !map.empty())
            cv::threshold (map, retMat, obsThreshold, FREE_PIXEL, cv::THRESH_BINARY_INV);
        else if (map_type==OBSTACLE_LABEL_MAP && !obsLabelMap.empty())
            obsLabelMap.copyTo (retMat);
        else if (map_type==ORIGINAL_MAP && !map.empty())
            map.copyTo (retMat);
        
        if (color_type == COLOR_MAP)
            cv::cvtColor (retMat, retMat, cv::COLOR_GRAY2RGB);
        
        if (draw_rep_pts && repPtsComputed)
            for (int i=0; i<repPts.size(); ++i) {
                cv::circle (retMat, repPts[i], 3, cv::Scalar(0,0,255), -1);
                retMat.at<cv::Vec3b> (repPts[i].y, repPts[i].x) = cv::Vec3b(0,255,0);
            }
        
        return (retMat);
    }
    
    // ===================================================
    
    cv::Mat obsLabelMap;
    std::vector<cv::Point> repPts;
    
    void computeRepresentativePoints (void) {
        repPts.clear();
        if (map.empty()) return;
        
        cv::threshold (map, obsLabelMap, obsThreshold, 1, cv::THRESH_BINARY_INV);
        int label_count = 2; // starts at 2 because 0, 1 are used already

        for (int y=0; y<obsLabelMap.rows; y++) {
            for (int x=0; x<obsLabelMap.cols; x++) {
                if (obsLabelMap.at<uchar>(y,x) != 1)
                    continue;

                cv::Rect rect;
                cv::floodFill (obsLabelMap, cv::Point(x,y), label_count, &rect, 0, 0, 4);
                
                if (MAX(rect.width,rect.height) < repPtObsSizeThresh) continue;
                
                std::vector<int> dists (rect.width, 0);
                int j = rect.y + rect.height/2;
                int ii, dd = 0;
                for(ii=0; ii<rect.width; ii++) {
                    if (obsLabelMap.at<uchar>(j,ii+rect.x) == label_count)  dd++;
                    else dd = 0;
                    dists[ii] = dd;
                }
                int bestdist=-1, bestii=-1;
                dd = 0;
                for(ii=rect.width-1; ii>=0; ii--) {
                    if (obsLabelMap.at<uchar>(j,ii+rect.x) == label_count)  dd++;
                    else dd = 0;
                    dists[ii] = MIN(dd, dists[ii]);
                    if (dd>0 && dists[ii]>bestdist) { bestdist = dists[ii]; bestii = ii; }
                }
                
                repPts.push_back (cv::Point(bestii+rect.x,j));
                label_count++;
            }
        }
    }
    
};

#endif
