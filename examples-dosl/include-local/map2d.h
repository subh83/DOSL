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

#ifndef __DOSL_MAP_2D_H
#define __DOSL_MAP_2D_H

#include <string>
#include <unordered_map>
#include <cassert>
#include <opencv/cv.h>
#include <opencv/cvaux.h>
#include <opencv/highgui.h> 

#include "mathevalAux.h"
#include <dosl/aux-utils/double_utils.hpp>

// =======================================

enum chart_limits_t { PIXEL_BOUNDARY, PIXEL_CENTER };

class PixelChartTransformation {
public:
    int width, height;
    double x_min, x_max, y_min, y_max;
    double dx, dy; // pixel size
    bool isEmpty;
    
    PixelChartTransformation (int w, int h, double cMinX, double cMaxX, double cMinY, double cMaxY, 
                                    int chart_limits=PIXEL_BOUNDARY) : 
            isEmpty(false), width(w), height(h), 
            x_min(cMinX), x_max(cMaxX), y_min(cMinY), y_max(cMaxY) {
        if (chart_limits == PIXEL_BOUNDARY) {
            dx = (x_max - x_min) / ((double)width);
            dy = (y_max - y_min) / ((double)height);
        }
        else {  // "pixel_center"
            dx = (x_max - x_min) / ((double)width-1.0);
            dy = (y_max - y_min) / ((double)height-1.0);
            x_min -= (dx/2.0);
            y_min -= (dy/2.0);
        }
    }
    
    PixelChartTransformation (int w, int h) {
        *this = PixelChartTransformation (w, h, -0.5, (double)w-0.5, -0.5, (double)h-0.5, PIXEL_BOUNDARY);
    }
    
    PixelChartTransformation () : isEmpty (true) { }
    
    // -------------------
    
    bool empty(void) const { return (isEmpty); }
    
    void print (std::string prefix) {
        std::cout << prefix << " width=" << width 
                            << ", height=" << height
                            << ", x_min=" << x_min
                            << ", x_max=" << x_max
                            << ", y_min=" << y_min
                            << ", y_max=" << y_max
                            << ", dx=" << dx
                            << ", dy=" << dy << std::endl;
    }
};

// =======================================
// Declarations

class ChartPoint;
class MapPoint;
class Pixel;


class ChartPoint {
public:
    double x, y;
    
    ChartPoint () { }
    template<class UserCoordType> ChartPoint (UserCoordType xx, UserCoordType yy) : x((double)xx), y((double)yy) { }
    template<class UserPointClass> ChartPoint (const UserPointClass& p) : x((double)p.x), y((double)p.y) { }
    
    MapPoint toMapPoint (const PixelChartTransformation& pct) const;
    Pixel toPixel (const PixelChartTransformation& pct) const;
    
    bool inDomain (const PixelChartTransformation& pct, bool include_max_val=true) const;
};

class MapPoint {
public:
    double x, y;
    
    MapPoint () { }
    template<class UserCoordType> MapPoint (UserCoordType xx, UserCoordType yy) : x((double)xx), y((double)yy) { }
    template<class UserPointClass> MapPoint (const UserPointClass& p) : x((double)p.x), y((double)p.y) { }
    
    ChartPoint toChartPoint (const PixelChartTransformation& pct) const;
    Pixel toPixel (const PixelChartTransformation& pct) const;
    
    bool inDomain (const PixelChartTransformation& pct, bool include_max_val=true) const;
};

class Pixel {
public:
    int x, y;
    
    Pixel () { }
    template<class UserCoordType> Pixel (UserCoordType xx, UserCoordType yy) : x((int)xx), y((int)yy) { }
    template<class UserPointClass> Pixel (const UserPointClass& p) : x((int)p.x), y((int)p.y) { }
    
    ChartPoint toChartPoint (const PixelChartTransformation& pct) const;
    MapPoint toMapPoint (const PixelChartTransformation& pct) const;
    
    bool inDomain (const PixelChartTransformation& pct, bool include_max_val=false) const;
};

// --------------------------------------
// Definitions
/*
  |-----|-----|-----|
  ^                     x_min of Chart
  ^                     -0.5 of Map
     ^                  0 of Map
        ^               0.5 of Map
  |-----|               0 of Pixel
*/

MapPoint ChartPoint::toMapPoint (const PixelChartTransformation& pct) const {
    MapPoint ret;
    ret.x = (x-pct.x_min)/pct.dx - 0.5;
    ret.y = (y-pct.y_min)/pct.dy - 0.5;
    return (ret);
}

ChartPoint MapPoint::toChartPoint (const PixelChartTransformation& pct) const {
    ChartPoint ret;
    ret.x = (x+0.5)*pct.dx + pct.x_min;
    ret.y = (y+0.5)*pct.dy + pct.y_min;
    return (ret);
}

Pixel MapPoint::toPixel (const PixelChartTransformation& pct) const {
    Pixel ret;
    
    ret.x = (int) round (x);
    if ( ret.x==-1 && isEqual_d(x,-0.5) ) ret.x = 0;
    else if ( ret.x==pct.width && isEqual_d(x,((double)pct.width)-0.5) ) ret.x = pct.width - 1;
    
    ret.y = (int) round (y);
    if ( ret.y==-1 && isEqual_d(y,-0.5) ) ret.y = 0;
    else if ( ret.y==pct.height && isEqual_d(y,((double)pct.height)-0.5) ) ret.y = pct.height - 1;
    
    return (ret);
}

Pixel ChartPoint::toPixel (const PixelChartTransformation& pct) const {
    return (toMapPoint(pct).toPixel(pct));
}

MapPoint Pixel::toMapPoint (const PixelChartTransformation& pct) const {
    MapPoint ret;
    ret.x = (double)x;
    ret.y = (double)y;
    return (ret);
}

ChartPoint Pixel::toChartPoint (const PixelChartTransformation& pct) const {
    return (toMapPoint(pct).toChartPoint(pct));
}

// ---------------------------------------

bool ChartPoint::inDomain (const PixelChartTransformation& pct, bool include_max_val) const {
    if (include_max_val)
        return ( isGreaterEq_d(x, pct.x_min) && isLessEq_d(x, pct.x_max) &&
                        isGreaterEq_d(y, pct.y_min) && isLessEq_d(y, pct.y_max) );
    else
        return ( isGreaterEq_d(x, pct.x_min) && isLess_d(x, pct.x_max) &&
                        isGreaterEq_d(y, pct.y_min) && isLess_d(y, pct.y_max) );
}

bool MapPoint::inDomain (const PixelChartTransformation& pct, bool include_max_val) const {
    if (include_max_val)
        return ( isGreaterEq_d(x, -0.5) && isLessEq_d(x, ((double)pct.width)-0.5) &&
                        isGreaterEq_d(y, -0.5) && isLessEq_d(y, ((double)pct.height)-0.5) );
    else
        return ( isGreaterEq_d(x, -0.5) && isLess_d(x, ((double)pct.width)-0.5) &&
                        isGreaterEq_d(y, -0.5) && isLess_d(y, ((double)pct.height)-0.5) );
}

bool Pixel::inDomain (const PixelChartTransformation& pct, bool include_max_val) const {
    if (include_max_val)
        return ( isGreaterEq_i(x,0) && isLessEq_i(x,pct.width) && isGreaterEq_i(y,0) && isLessEq_i(y,pct.height) );
    else
        return ( isGreaterEq_i(x,0) && isLess_i(x,pct.width) && isGreaterEq_i(y,0) && isLess_i(y,pct.height) );
}


// =======================================
// =======================================

class map2d {
public:
    PixelChartTransformation pc_transform;
    
    std::vector< std::vector< double > > map;
    MathEvaluator fun;
    double val;
    
    // ---------------
    map2d () { }
    
    map2d (int width, int height, const PixelChartTransformation& pct=PixelChartTransformation(), double default_val=0.0 ) {
        if (pct.empty())  pc_transform = PixelChartTransformation (width, height);
        else  pc_transform = pct;
        assert ((pc_transform.width==width && pc_transform.height==height));
        
        map.resize (width);
        for (int x=0; x<width; ++x) {
            map[x].resize (height);
            for (int y=0; y<height; ++y)
                map[x][y] = default_val;
        }
    }
    
    map2d (const PixelChartTransformation& pct, double default_val=0.0 ) {
        *this = map2d (pct.width, pct.height, pct, default_val);
    }
        
    map2d (double uniform_val) : val(uniform_val) { }
    
    // pixel map
    map2d (std::string imagefname, const PixelChartTransformation& pct=PixelChartTransformation(), 
                                                            double pix0val=0.0, double pix255val=1.0) { 
        cv::Mat image = cv::imread (imagefname, CV_LOAD_IMAGE_GRAYSCALE);
        
        if (pct.empty())  pc_transform = PixelChartTransformation (image.cols, image.rows);
        else  pc_transform = pct;
        
        assert ((pc_transform.width==image.cols && pc_transform.height==image.rows));
        
        map.resize (pc_transform.width);
        for (int x=0; x<pc_transform.width; ++x) {
            map[x].resize (pc_transform.height);
            for (int y=0; y<pc_transform.height; ++y)
                map[x][y] = pix0val + (pix255val-pix0val)*((double)image.at<uchar>(y,x))/255.0;
        }
    }
    
    // function that uses chart coodinates
    map2d (const MathEvaluator& ff, const PixelChartTransformation& pct=PixelChartTransformation(), 
                    bool reduce_to_pix_map=true) { 
        fun = ff;
        pc_transform = pct;
        
        if (reduce_to_pix_map && !pct.empty()) {
            ChartPoint coords;
            map.resize (pct.width);
            for (int x=0; x<pct.width; ++x) {
                map[x].resize (pct.height);
                for (int y=0; y<pct.height; ++y) {
                    coords = Pixel(x,y).toChartPoint (pc_transform); //pc_transform.pix2chart(x,y);
                    map[x][y] = fun.eval( {{"x",coords.x}, {"y",coords.y}} );
                    //std::cout << "("<<x<<","<<y<<"): " << map[x][y] << "; ";
                }
            }
        }
    }
    
    
    // ------------------------
    
    double evalChartPoint (const ChartPoint& p) { // use:  M.evalChartPoint (ChartPoint (x,y))
        if (map.size()>0) {
            Pixel pix = p.toPixel (pc_transform);
            return (map[pix.x][pix.y]); // TODO: interpolate
        }
        else if (!fun.empty())
            return ( fun.eval ( {{"x",p.x},{"y",p.y}} ) );
        else
            return (val);
    }
    
    double evalMapPoint (const MapPoint& p) { // use:  M.MapPoint (MapPoint (x,y))
        if (map.size()>0) {
            Pixel pix = p.toPixel (pc_transform);
            return (map[pix.x][pix.y]); // TODO: interpolate.
        }
        else if (!fun.empty()) {
            ChartPoint pc = p.toChartPoint (pc_transform);
            return ( fun.eval ( {{"x",pc.x},{"y",pc.y}} ) );
        }
        else
            return (val);
    }
    
    double evalPixel (const Pixel& p) { // use:  M.Pixel (Pixel (x,y))
        if (map.size()>0) {
            return (map[p.x][p.y]);
        }
        else if (!fun.empty()) {
            ChartPoint pc = p.toChartPoint (pc_transform);
            return ( fun.eval ( {{"x",pc.x},{"y",pc.y}} ) );
        }
        else
            return (val);
    }
    
    /*double evalAtCoord (double px, double py) {
        if (map.size()>0) { // TODO: interpolate.
            MapPoint<int> pix = pc_transform.chart2pix_i (px, py);
            return (map[pix.x][pix.y]);
        }
        else if (!fun.empty())
            return ( fun.eval ( {{"x",px},{"y",py}} ) );
        else
            return (val);
    }
    
    template<class PointType>
    double evalAtCoord (PointType p) {
        return (evalAtCoord(p.x, p.y));
    }
        
    // ------------------------
    
    template<class PointCoordType>
    double evalAtPix (PointCoordType px, PointCoordType py) {
        if (map.size()>0) { // TODO: interpolate.
            
            return (map[(int)px][(int)py]);
        }
        else if (!fun.empty()) {
            MapPoint<double> p = pc_transform.pix2chart (px, py);
            return ( fun.eval ( {{"x",p.x},{"y",p.y}} ) );
        }
        else
            return (val);
    }
    
    template<class PointType>
    double evalAtPix (PointType pix) {
        return (evalAtPix(pix.x, pix.y));
    }*/
    
    // ------------------------
    
    double& at (int x, int y) {
        assert (map.size()>0);
        return (map[x][y]);
    }
    
};


#endif
