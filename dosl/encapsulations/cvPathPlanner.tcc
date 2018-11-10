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
#ifndef __CV_PATH_PLANNER_TCC
#define __CV_PATH_PLANNER_TCC

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

// Other libraries:
// Open CV:
#include <opencv2/opencv.hpp> 
#include <opencv2/highgui.hpp> 
// DOSL library
#include "../dosl"
// local
#include "../aux-utils/cvParseMap2d.hpp"
#include "../aux-utils/double_utils.hpp"


#ifndef __DOSL_CV_ENUMS
#define SQUARE_GRID 4  // o1o
#define TRIANGULATION_CONNECTIVITY 8 // 1oo
enum cv_graph_connectivity_type { graph_connectivity_UNDEFIND = 0, 
                                  FOUR_CONNECTED = SQUARE_GRID, EIGHT_CONNECTED = SQUARE_GRID + 1, // 010, 011
                                  UNIFORM_TRIANGULATION = TRIANGULATION_CONNECTIVITY, // 100 (unsupported)
                                  SQUARE_TRIANGULATION = (SQUARE_GRID | TRIANGULATION_CONNECTIVITY) // 110
                                  };
enum cv_vetrex_location_type {vertex_location_UNDEFIND=0, PIXEL_CORNER=1, PIXEL_CENTER=2};
enum cv_heuristic_type { heuristic_UNDEFINED = 0,
                         ZERO_HEURISTIC, // Dijkstra's
                         EUCLIDEAN_HEURISTIC, 
                         EIGHT_CONNECTED_HEURISTIC,
                         MANHATTAN_HEURISTIC // 
                         };
#define __DOSL_CV_ENUMS
#endif

// ==============================================================================

// A node of the graph
template <class doslAlgorithm>
class cvPPnode : public doslAlgorithm::template Node< cvPPnode<doslAlgorithm>, double>
{
public:
    int x, y;
    
    void put_in_grid (void) { }
    
    bool isCoordsEqual(const cvPPnode& n) const {
        return (isEqual_i(x,n.x) && isEqual_i(y,n.y));
    }
    
    // *** This must be defined for the node
    bool operator==(const cvPPnode& n) const {
        if (!isCoordsEqual(n)) return (false);
        return (true);
    }
    
    // constructor
    cvPPnode () { }
    cvPPnode (int xx, int yy) : x(xx), y(yy) { put_in_grid(); }
    
    
    // Inherited functions being overwritten
    int getHashBin (void) const {
        return (abs(x));
    }
    
    // print
    void print (std::string head="", std::string tail="\n") const {
        /*_dosl_cout << _GREEN + head << " (" << this << ")" GREEN_ " x=" << x << ", y=" << y << "; ";
        (this->G == std::numeric_limits<double>::max())? printf("INF") : printf("G = %0.8f \n", this->G);
        _dosl_cout << "Successors: ";
        for (auto it=this->Successors.begin(); it!=this->Successors.end(); ++it)
            printf ("%x (%f), ", it->first, it->second);
        std::cout << tail << _dosl_endl;*/
    }
};

// ==============================================================================

template <class doslAlgorithm>
class cvPathPlanner : public doslAlgorithm::template Algorithm< cvPPnode<doslAlgorithm>, double>
{
public:
    typedef cvPPnode<doslAlgorithm> cvPPnode_;
    
    cvParseMap2d  my_map;
    cv_graph_connectivity_type GRAPH_TYPE; // 4, 6 or 8
    cv_vetrex_location_type VERTEX_LOCATION; // PIXEL_CORNER or PIXEL_CENTER
    cv_heuristic_type HEURISTIC;
    
    // Image display variables / parameters
    bool visualize;
    cv::Mat image_to_display;
    double PLOT_SCALE;
    double VERTEX_SIZE, LINE_THICKNESS;
    int VIS_INTERVAL, VERTEX_COLORS;
    
    // variables for saving image
    int frameno;
    std::ostringstream imgPrefix;
    
    // variables decsribing problem
    int WIDTH, HEIGHT;
    cvPPnode_ startNode, goalNode, lastExpanded;
    
    // =========================================
    
    template<typename T>
    cv::Point cv_plot_coord(T x, T y) { return( cv::Point( round(PLOT_SCALE*(x-0)), round(PLOT_SCALE*(y-0)) ) ); }
    
    void cvPlotPoint (cv::Point pt, CvScalar colr, int size=1) {
        for (int i=pt.x; i<=pt.x+(size-1); i++)
            for (int j=pt.y; j<=pt.y+(size-1); j++) {
                if (i<0 || i>=image_to_display.cols || j<0 || j>=image_to_display.rows) continue;
                image_to_display.at<cv::Vec3b>(j,i) = cv::Vec3b ((uchar)colr.val[0], (uchar)colr.val[1], (uchar)colr.val[2]);
            }
    }
    
    // =========================================
    
    std::vector<cv::Point> path;
    double cost;
    
    
    void find_path (cv::Point start, cv::Point goal, bool vis=false) {
        // set start and goal
        startNode = cvPPnode_ (start.x, start.y);
        goalNode = cvPPnode_ (goal.x, goal.y);
        startNode.put_in_grid(); goalNode.put_in_grid();
        visualize = vis;
        
        if (visualize && image_to_display.empty()) {
            image_to_display = my_map.getCvMat (COLOR_MAP);
            cv::resize (image_to_display, image_to_display, cv::Size(), PLOT_SCALE , PLOT_SCALE );
            cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE);
            cv::imshow("Display window", image_to_display);
            //cv::waitKey(0);
        }
        
        // ------------
        this->search();
        // ------------
        
        auto dosl_path = this->reconstructPath (goalNode);
        double cost = 0.0;
        CvScalar colr = cvScalar (50.0 + rand()%150, 50.0 + rand()%150, 50.0 + rand()%150);
        
        cv::Point2d lastPt, thisPt(startNode.x,startNode.y);
        for (int a=dosl_path.size()-1; a>=0; --a) {
            lastPt = thisPt;
            thisPt = cv::Point2d(0.0, 0.0);
            for (auto it=dosl_path[a].begin(); it!=dosl_path[a].end(); ++it) {
                thisPt.x += it->second * it->first->x;
                thisPt.y += it->second * it->first->y;
            }
            path.push_back (cv::Point((int)round(thisPt.x),(int)round(thisPt.y)));
            if (visualize)
                cv::line (image_to_display, cv_plot_coord(thisPt.x,thisPt.y), cv_plot_coord(lastPt.x,lastPt.y), 
                                        colr, LINE_THICKNESS*PLOT_SCALE);
            cost += sqrt ((thisPt.x-lastPt.x)*(thisPt.x-lastPt.x) + (thisPt.y-lastPt.y)*(thisPt.y-lastPt.y));
        }
        
        if (visualize) {
            cv::imshow("Display window", image_to_display);
            cv::waitKey(0);
        }
    }
    
    // -----------------------------------------------------------
    
    // Constructor
    cvPathPlanner ( cv::Mat obs_map, 
                    cv_graph_connectivity_type graph_type=graph_connectivity_UNDEFIND, 
                    cv_vetrex_location_type vertex_location=vertex_location_UNDEFIND,
                    cv_heuristic_type heuristic_type=heuristic_UNDEFINED )
                            : GRAPH_TYPE(graph_type), VERTEX_LOCATION(vertex_location), HEURISTIC(heuristic_type)
    {
        // read data for planning
        my_map = cvParseMap2d (obs_map, false);
        WIDTH = my_map.width(); HEIGHT = my_map.height();
        
        
        // set graph type
        if (GRAPH_TYPE == UNIFORM_TRIANGULATION) {
            _dosl_info("GRAPH_TYPE=UNIFORM_TRIANGULATION is currently unsupported in pixel-based map planning.");
            GRAPH_TYPE = graph_connectivity_UNDEFIND;
        }
        
        // enforce planner requirements
        if (this->algorithm_name()=="SStar") {
            if (!(GRAPH_TYPE & TRIANGULATION_CONNECTIVITY)) {
                if (GRAPH_TYPE != graph_connectivity_UNDEFIND)
                    _dosl_info("SStar Algorithm: Setting GRAPH_TYPE=SQUARE_TRIANGULATION (SStar do not support FOUR_CONNECTED graphs).");
                GRAPH_TYPE = SQUARE_TRIANGULATION;
                            /* TODO: path reconstruction does not work with EIGHT_CONNECTED because of degenerate 3-simplices being rejected
                                     in computation of 'getAllAttachedMaximalSimplices' (during path reconstruction only?).
                                     Need to fix this by allowing faces of the degenerate simplex to be returned instead.
                                     TODO: Create function 'getAllAttachedNondegenerateMaximalSimplices'. */
             }
             HEURISTIC = ZERO_HEURISTIC; // Non-zero heuistic such as EUCLIDEAN_HEURISTIC is not fully supported by SStar yet.
        }
        
        if (GRAPH_TYPE == 0) GRAPH_TYPE=EIGHT_CONNECTED;
        
        if (HEURISTIC == 0) {
            if (GRAPH_TYPE == FOUR_CONNECTED)  HEURISTIC = MANHATTAN_HEURISTIC;
            else if (GRAPH_TYPE & SQUARE_GRID) HEURISTIC = EIGHT_CONNECTED_HEURISTIC;
            else HEURISTIC = EUCLIDEAN_HEURISTIC;
        }
        
        if (VERTEX_LOCATION == 0) VERTEX_LOCATION = PIXEL_CENTER;
        
        
        // display options
        PLOT_SCALE = 1.0;
        VERTEX_SIZE = 1.0;
        LINE_THICKNESS = 2.0;
        VERTEX_COLORS = 1; // 1 or 0
        VIS_INTERVAL = 100;
        
        // saving options
        //frameno = 0;
        //imgPrefix << MAKESTR(_DOSL_ALGORITHM) << GRAPH_TYPE << "homotopy2d_";
        
        //visualize = vis;
        
        // Set planner variables
        this->AllNodesSet.HashTableSize = ceil(WIDTH - 0 + 1);
        
        // compute path
        // find_path();
    }
    
    // -----------------------------------------------------------
    
    bool isNodeInWorkspace (const cvPPnode_& tn) {
        if (VERTEX_LOCATION == PIXEL_CORNER)
            return (tn.x>=0 && tn.y>=0 && tn.x<=WIDTH && tn.y<=HEIGHT);
        else if (VERTEX_LOCATION == PIXEL_CENTER)
            return (tn.x>=0 && tn.y>=0 && tn.x<WIDTH && tn.y<HEIGHT);
        return (true);
    }
    
    bool isEdgeAccessible (const cvPPnode_& tn1, const cvPPnode_& tn2) {
        if ( (!isNodeInWorkspace(tn1)) || !(isNodeInWorkspace(tn2)) )  return (false);
        
        if (VERTEX_LOCATION==PIXEL_CORNER && (GRAPH_TYPE & SQUARE_GRID)) {
            if (tn1.x!=tn2.x && tn1.y!=tn2.y) // diagonal edge. 8-connected grid.
                return ( my_map.isFree( MIN(tn1.x,tn2.x), MIN(tn1.y,tn2.y) ) );
            else if (tn1.x==tn2.x) // need to check two cells
                /* return (!( (tn1.x==WIDTH || my_map.isObstacle ( round(tn1.x), round(MIN(tn1.y,tn2.y)) ) ) &
                            ( tn1.x==0 || my_map.isObstacle ( round(tn1.x-1), round(MIN(tn1.y,tn2.y)) ) ) )); */
                /*return ( (tn1.x==WIDTH || my_map.isFree(tn1.x, MIN(tn1.y,tn2.y))) &&    // right cell is free
                            (tn1.x==0 || my_map.isFree(tn1.x-1, MIN(tn1.y,tn2.y)))  ); // left cell is free */
                /*return ( ( my_map.isFree(tn1.x, MIN(tn1.y,tn2.y)) ) ||    // right cell is free
                            ( my_map.isFree(tn1.x-1, MIN(tn1.y,tn2.y)) )  ); // left cell is free */
                return ( ( tn1.x>0 && my_map.isFree(tn1.x-1,MIN(tn1.y,tn2.y)) ) ||
                            ( tn1.x<WIDTH && my_map.isFree(tn1.x,MIN(tn1.y,tn2.y)) ) );
            else if (tn1.y==tn2.y) // need to check two cells (above and below)
                /*return (!( ( tn1.y==HEIGHT || my_map.isObstacle ( round(MIN(tn1.x,tn2.x)), round(tn1.y) ) ) &
                            ( tn1.y==0 || my_map.isObstacle ( round(MIN(tn1.x,tn2.x)), round(tn1.y-1) ) ) ));*/
                return ( ( tn1.y>0 && my_map.isFree(MIN(tn1.x,tn2.x),tn1.y-1) ) ||
                            ( tn1.y<HEIGHT && my_map.isFree(MIN(tn1.x,tn2.x),tn1.y) ) );
        }
        else if (GRAPH_TYPE & SQUARE_GRID)
            return ( my_map.isFree(tn1.x,tn1.y)  && my_map.isFree(tn2.x,tn2.y) );
        else
            // TODO: make more precise for UNIFORM_TRIANGULATION
            return ( my_map.isFree(approx_floor(tn1.x),approx_floor(tn1.y))  && 
                        my_map.isFree(approx_floor(tn2.x),approx_floor(tn2.y)) );
        
        return (true);
    }
    
    // -----------------------------------------------------------
    
    void getSuccessors (cvPPnode_ &n, std::vector<cvPPnode_>* s, std::vector<double>* c) // *** This must be defined
    {
        // This function should account for obstacles and size of environment.
        cvPPnode_ tn;
        
        if (GRAPH_TYPE & SQUARE_GRID) {
            for (int a=-1; a<=1; ++a)
                for (int b=-1; b<=1; ++b) {
                    if (a==0 && b==0) continue;
                    
                    if (GRAPH_TYPE==FOUR_CONNECTED  &&  (a!=0 && b!=0)) continue; // diagonal
                    
                    if (GRAPH_TYPE & TRIANGULATION_CONNECTIVITY) { // use triangulation for 'SStar'
                        int xParity = ((int)round(fabs(n.x))) % 2;
                        if (xParity==0 && (a!=0 && b==-1)) continue;
                        if (xParity==1 && (a!=0 && b==1)) continue;
                    }
                    
                    tn.x = n.x + a;
                    tn.y = n.y + b;
                    
                    if (!isEdgeAccessible(tn,n)) continue;
                    
                    //updateHSignature(n,tn);
                    
                    s->push_back(tn);
                    double dx=tn.x-n.x, dy=tn.y-n.y;
                    c->push_back(sqrt(dx*dx+dy*dy)); 
                }
        }
        /*else if (GRAPH_TYPE == UNIFORM_TRIANGULATION) {
            double th;
            for (int a=0; a<6; ++a) {
                th = a * PI_BY_3;
                tn.x = n.x + 1.0*cos(th);
                tn.y = n.y + 1.0*sin(th);
                
                if (!isEdgeAccessible(tn,n)) continue;
                
                //updateHSignature(n,tn);
                    
                s->push_back(tn);
                c->push_back(1.0);
            }
        
        }*/
        
    }
    
    // -----------------------------------------------------------
    // 'isSegmentFree' is required by ThetaStar
    
    bool isSegmentFree (cvPPnode_ &n1, cvPPnode_ &n2, double* c)
    {
        double dx=(double)(n2.x-n1.x), dy=(double)(n2.y-n1.y);
        *c = sqrt(dx*dx + dy*dy);
        
        double step_count = 2.0 * ceil(*c);
        double xstep=dx/step_count, ystep=dy/step_count;
        
        double xx, yy;
        for (int a=0; a<step_count; ++a) {
            xx = n1.x + a*xstep; yy = n1.y + a*ystep;
            if (my_map.isObstacle ( (int)round(xx), (int)round(yy) ) )
                return (false);
        }
        return (true);
    }
    
    // -----------------------------------------------------------
    
    double getHeuristics (cvPPnode_& n)
    {
        if (HEURISTIC == ZERO_HEURISTIC)
            return (0.0);
        
        double dx = fabs(goalNode.x - n.x);
        double dy = fabs(goalNode.y - n.y);
        
        if (HEURISTIC == EUCLIDEAN_HEURISTIC)
            return (sqrt(dx*dx + dy*dy));
        
        if (HEURISTIC == EIGHT_CONNECTED_HEURISTIC)
            return (fabs(dx-dy) + SQRT2*MIN(dx,dy));
        
        if (HEURISTIC == MANHATTAN_HEURISTIC)
            return (dx + dy);
    }
    
    // -----------------------------------------------------------
    
    std::vector<cvPPnode_> getStartNodes (void) 
    {
        std::vector<cvPPnode_> startNodes;
        
        startNodes.push_back (startNode);
        if (visualize)
            cv::circle (image_to_display, cv_plot_coord(startNode.x,startNode.y), VERTEX_SIZE*PLOT_SCALE, 
                                                                cvScalar (200.0, 150.0, 150.0), -1, 8);
        
        return (startNodes);
    }
    
    // -----------------------------------------------------------
    
    void nodeEvent (cvPPnode_ &n, unsigned int e) 
    {
        if (!visualize) return;
        
        CvScalar col = cvScalar(0.0, 0.0, 0.0);;
        int thickness = -1; //lineThickness;
        double radFactor = 1.0;
        
        bool pauseForVis=false, drawVertex=true;
        
        // --------------------------------------------
        if (e & this->EXPANDED) {
            lastExpanded = n;
            if (VERTEX_COLORS)
                col = cvScalar(200.0, 255.0, 200.0);
            else
                col = cvScalar(255.0, 255.0, 255.0);
        }
        
        else if ((e & this->HEAP) == this->PUSHED) {
            if (VERTEX_COLORS)
                col = cvScalar(255.0, 0.0, 0.0); // blue
            else
                col = cvScalar(200.0, 200.0, 200.0);
        }
        
        else if (e & this->UNEXPANDED) {
            col = cvScalar(255.0, 0.0, 150.0); // purple
            printf ("Backtrack started!!\n");
        }
        
        else
            return;
                
        //-------------------------------------------
        if (drawVertex) {
            double nodeRad = radFactor*VERTEX_SIZE;
            if (n.x<WIDTH  &&  n.y<HEIGHT  &&  my_map.isFree(round(n.x), round(n.y)) )
                cvPlotPoint (cv_plot_coord(n.x,n.y), col, PLOT_SCALE);
        }
        
        if (this->ExpandCount % VIS_INTERVAL == 0  ||  this->NodeHeap.size() == 0) {
            cv::imshow("Display window", image_to_display);
            std::cout << std::flush;
            /* if (SAVE_IMG_INTERVAL>0 && ExpandCount%SAVE_IMG_INTERVAL == 0) {
                char imgFname[1024];
                sprintf(imgFname, "%s%s%05d.png", out_folderName.c_str(), imgPrefix.str().c_str(), ExpandCount);
                cv::imwrite(imgFname, image_to_display);
            } */
            cvWaitKey(1); //(10);
        }
        
        if (pauseForVis) {
            cvWaitKey();
        }
    }
    
    // ---------------------------------------
    
    bool stopSearch (cvPPnode_ &n) {
        return (n == goalNode);
    }
    
    // =========================================
    
    void draw_path (cv::Mat& in_map, cv::Mat& out_map, CvScalar color=cvScalar(0.0,0.0,255.0),
                                            int thickness=1, int lineType=8, int shift=0) {
        int nChannels = out_map.channels();
        
        if (&in_map != &out_map) {
            in_map.copyTo (out_map);
            if (out_map.channels()<3 && nChannels==3)
                cv::cvtColor (out_map, out_map, CV_GRAY2RGB);
        }
        
        for (int j=0; j<path.size()-1; ++j)
            cv::line (out_map, path[j], path[j+1], color, thickness, lineType, shift);
    }
    
    
    cv::Mat draw_path (CvScalar color=cvScalar(0.0,0.0,255.0), int thickness=1, int lineType=8, int shift=0) {
        cv::Mat ret = my_map.getCvMat (COLOR_MAP, FREE_MAP, 0);
        draw_path (ret, ret, color, thickness, lineType, shift);
        return (ret);
    }
};

// ==============================================================================


#endif
