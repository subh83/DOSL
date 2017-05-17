/** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal search Library (DOSL)                                              *
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
#include <fstream>
#include <limits>
#include <csignal>
#include <vector>

// Other libraries:
// Open CV:
#include <opencv/cv.h>
#include <opencv/cvaux.h>
#include <opencv/highgui.h>

// DOSL library

// Optional parameters -- to be declared before including dosl: 
// #define _DOSL_DEBUG 1                // 0,1, or 2
// #define _DOSL_AUTOCORRECT 1          // 0 or 1
// #define _DOSL_VERBOSE_LEVEL 0        // 0 to 5
// #define _DOSL_VERBOSE_ITEMS "NONE"   // "fun1,fun2" or "NONE" or "ALL"
           // Common items: "reconstructPath,findCamefromPoint,getAllAttachedMaximalSimplices,getAllMaximalSimplicesFromSet"
// #define _DOSL_PRINT_COLORS 1
// #define _DOSL_EVENTHANDLER 1

// Include main DOSL header
#include <dosl/dosl>

// The following allows us to use the 'DOSL_CLASS' macro.
#ifndef _DOSL_ALGORITHM // Allows us to pass at command line during compilation: g++ ... -D_DOSL_ALGORITHM=AStar
    #define _DOSL_ALGORITHM  AStar
#endif

// =======================

#define PI       3.14159265359
#define PI_BY_3  1.0471975512
#define SQRT3    1.73205080757
#define SQRT3BY2 0.86602540378

#define SMALL_EPS  1e-8
// relaxed comparisons (appropriate if x and y take discrete values)
#define isEqual_d(x,y)  ( fabs((x)-(y)) < SMALL_EPS ) // x == y
#define isLess_d(x,y)  ( (x)+SMALL_EPS < (y) ) // x < y
#define isGreater_d(x,y)  ( (x) > (y)+SMALL_EPS ) // x < y

#define sign(x) ((x>0.0)?1.0:((x<0.0)?-1.0:0.0))

#define _VIS 1
#define VIS_INTERVAL 10
#define SAVE_IMG_INTERVAL -1 // 0 to not save at all. -1 to save last frame only.
#define VERTEX_COLORS 1
#define DRAW_INITIAL_GRID 1

#define GRAPH_TYPE 6 // 6 or 8

// GRAPH_TYPE == 8
#define GRAPH8_X_WIGGLE 0  // 0 or 1
#define GRAPH8_X_WIGGLE_AMOUNT 0.1  // should be <0.5

// GRAPH_TYPE == 6
#define ENABLE_MESH_REFINEMENT  1 // 0 or 1
#define REFINE_DISK_RAD  6.0
#define REFINE_N_LAYERS  4.0

// ==============================================================================
// globals

double MAX_X=300.0, MIN_X=-300.0, MAX_Y=300.0, MIN_Y = -300.0;
std::vector< std::vector<double> > OBS_RECT = { {69.99, -70.01, 170.01, 70.01} };
                                    // Rectangles: {x1, y1, x2, y2}, x1<x2, y1<y2.  add little offset/wiggle

// ==============================================================================

// A node of the graph
class myNode : public DOSL_CLASS(Node)<myNode,double>  // Expands to 'SStarNode' or 'AStarNode'
{
public:
    double x, y;
    
    #if GRAPH_TYPE == 8
        #if GRAPH8_X_WIGGLE == 0
        void put_in_grid (void) { x = (double)round(x); y = (double)round(y); }
        #else
        void put_in_grid (void) {
            x = (double)round(x); y = (double)round(y);
            if ( ((int)fabs(y))%2 == 1) 
                x += GRAPH8_X_WIGGLE_AMOUNT;
        }
        #endif
    #elif GRAPH_TYPE == 6
    void put_in_grid (void) {
        int yLevel = round (y / SQRT3BY2);
        y = SQRT3BY2*yLevel;
        if (yLevel%2 == 0) // even
            x = round(x+1e-6);
        else
            x = round(x+0.5+1e-6) - 0.5;
    }
    
    #endif
    
    myNode (double xx, double yy) : x(xx), y(yy) { put_in_grid(); }
    bool operator==(const myNode& n) const { return ( isEqual_d(x,n.x) && isEqual_d(y,n.y) ); }; // This must be defined for the node
    
    // constructor
    myNode () { }
    
    // Functions being overwritten
    int getHashBin (void) const { return ( MAX(round(fabs(x)+1e-6), round(fabs(x)-1e-6)) ); }
    
    double dist (void) const { /* { return (sqrt(x*x+y*y)); } */
        // Assume single rectangular obstacle symmetric about x-axis
        double xMin = MIN(OBS_RECT[0][0], OBS_RECT[0][2]);
        double xMax = MAX(OBS_RECT[0][0], OBS_RECT[0][2]);
        double yAbs = fabs(OBS_RECT[0][1]);
        
        #if GRAPH_TYPE == 8  // &&  GRAPH8_X_WIGGLE == 0
        xMin = (double)floor(xMin); xMax = (double)ceil(xMax);
        yAbs = (double)ceil(yAbs);
        #elif GRAPH_TYPE == 6
        myNode top_left (xMin-0.5, yAbs+SQRT3BY2/2); top_left.put_in_grid();
        myNode top_right (xMax+0.5, yAbs+SQRT3BY2/2); top_right.put_in_grid();
        xMin = top_left.x+1e-8; yAbs = top_left.y-1e-8; xMax = top_right.x-1e-8;
        #endif
        
        if (x > xMin  &&  fabs(y) < yAbs*x/xMin) {
            double px, py, ret = (sqrt(xMin*xMin + yAbs*yAbs));
            if (fabs(y) < yAbs) {
                ret += xMax - xMin;
                px = xMax; py = sign(y)*yAbs;
            }
            else {
                px = xMin; py = sign(y)*yAbs;
            }
            return (ret + sqrt((x-px)*(x-px) + (y-py)*(y-py)));
        }
        else
            return (sqrt(x*x+y*y));
    }
    
    void print (std::string head="", std::string tail="") const {
        _dosl_cout << _GREEN + head << " (" << this << ")" GREEN_ "x=" << x << ", y=" << y << ", "; //dist=" << dist() << ". G=" << G; 
        printf("dist=%0.8f. G=", dist());
        (G==std::numeric_limits<double>::max())? printf("INF") : printf("%0.8f", G);
        printf(" (diff=%e)", G - dist());
        _dosl_cout << _dosl_endl;
        _dosl_printf_nobreak("Successors: ");
        for (auto it=Successors.begin(); it!=Successors.end(); ++it)
            printf ("%x, ", it->first);
        _dosl_cout << tail << _dosl_endl;
    }
};

// ==============================================================================

class inversionDiskInfo {
public:
    double cx, cy; // center
    double R, numLayers; // disk radius is R
    // derived
    double rInvThresh, rThresh;
    // parameters
    double k;
    
    // -------------------------------
    
    inversionDiskInfo (myNode c, double RR, double nn) : cx(c.x), cy(c.y), R(RR), numLayers(nn) {
        k = 0.5;
        rInvThresh = R + numLayers;
        rThresh = f2inv (rInvThresh);
    }
    
    // -------------------------------
    
    double f2 (double hexDist) {
        return (R + (R - hexDist)/k);
    }
    
    double f2inv (double hexDist) {
        return (R - k*(hexDist - R));
    }
    
    // -------------------------------
    
    double g1 (double hexDist) {
        double newHexDist;
        if ( isLess_d(hexDist, rThresh) )
            newHexDist = rInvThresh * hexDist / rThresh; // f1
        else if ( isLess_d(hexDist, R) ) 
            newHexDist = f2 (hexDist); // f2
        else if ( !isGreater_d(hexDist, rInvThresh) ) // <=
            newHexDist = f2inv (hexDist); // f2 inverse
        else
            newHexDist = std::numeric_limits<double>::max();
        
        return (newHexDist);
    }
    
    double g2 (double hexDist) {
        double newHexDist;
        if ( !isGreater_d(hexDist, rInvThresh) ) // <=
            newHexDist = rThresh * hexDist / rInvThresh; // f1 inverse
        else
            newHexDist = std::numeric_limits<double>::max();
        
        return (newHexDist);
    }
    
    // -------------------------------
    
    double getHexDist (myNode n) {  // hexagonal metric.
        double adx=fabs(n.x-cx), ady=fabs(n.y-cy);
        return (MAX(adx + ady/SQRT3, ady/SQRT3BY2));
    }
    
    std::vector<myNode> get_inverted_pt (myNode n, double* hexDist_p = NULL, std::vector<double>* newHexDists_p = NULL) {
        
        myNode invpt;
        invpt.x = n.x; invpt.y = n.y; // copy coordinates of n
        double scaleFactor, newHexDist, hexDist = getHexDist(n);
        std::vector<myNode> ret;
        
        if (hexDist_p)  *hexDist_p = hexDist;
        if (newHexDists_p)  newHexDists_p->clear();
        
        newHexDist = g1 (hexDist);
        scaleFactor = newHexDist / hexDist;
        invpt.x = cx + (n.x-cx)*scaleFactor; invpt.y = cy + (n.y-cy)*scaleFactor;
        ret.push_back (invpt);
        if (newHexDists_p) newHexDists_p->push_back (newHexDist);
        
        newHexDist = g2 (hexDist);
        scaleFactor = newHexDist / hexDist;
        invpt.x = cx + (n.x-cx)*scaleFactor; invpt.y = cy + (n.y-cy)*scaleFactor;
        ret.push_back (invpt);
        if (newHexDists_p) newHexDists_p->push_back (newHexDist);
        
        return (ret);
    }
};

// ---------------------------------------------------

class SearchProblem : public DOSL_CLASS(Problem)<myNode,double>  // SStarProblem or AStarProblem
{
public:
    // Image display variables / parameters
    cv::Mat image_to_display;
    double PLOT_SCALE, PROBLEM_SCALE;
    double VERTEX_SIZE, LINE_THICKNESS;
    
    // variables for saving image
    int frameno;
    std::ostringstream imgPrefix;
    
    // variables decsribing problem
    myNode startNode, goalNode, lastExpanded;
    
    // refinement options
    std::vector<inversionDiskInfo> refiners;
    
    // --------------------------------------------
    
    template<typename T>
    CvPoint cv_plot_coord(T x, T y) { return( cvPoint( round(PLOT_SCALE*(x-MIN_X)), round(PLOT_SCALE*(y-MIN_Y)) ) ); }
    
    // Constructors
    SearchProblem ()
    {
        // variables
        PLOT_SCALE = 10.0; //20.0; //10.0; //10.0; //8.0; //2.0; // 20;
        PROBLEM_SCALE = 0.1; //0.05; //0.1; //0.2; //0.4; //1.0; // 0.1;
                    // ^^ must have: PLOT_SCALE * PROBLEM_SCALE = 1.0 for window size to be as desired.
        
        MAX_X*=PROBLEM_SCALE; MIN_X*=PROBLEM_SCALE;
        MAX_Y*=PROBLEM_SCALE; MIN_Y*=PROBLEM_SCALE;
        for (int a=0; a<OBS_RECT.size(); a++) {
            OBS_RECT[a][0]*=PROBLEM_SCALE; OBS_RECT[a][1]*=PROBLEM_SCALE; 
            OBS_RECT[a][2]*=PROBLEM_SCALE; OBS_RECT[a][3]*=PROBLEM_SCALE; 
            #if ENABLE_MESH_REFINEMENT
            // select one corner of obstacle for mesh refinement
            refiners.push_back ( inversionDiskInfo (myNode (OBS_RECT[a][0], OBS_RECT[a][1]), REFINE_DISK_RAD, REFINE_N_LAYERS) );
            #endif
        }
        
        LINE_THICKNESS = 1; // CV_FILLED
        VERTEX_SIZE = 0.2; // 0.3; //0.5;
        
        // Initiation
        frameno = 0;
        imgPrefix << MAKESTR(_DOSL_ALGORITHM) << GRAPH_TYPE << "demo_";
        
        // Draw
        #if _VIS
        image_to_display.create (ceil( PLOT_SCALE * (MAX_Y - MIN_Y) ), ceil( PLOT_SCALE * (MAX_X - MIN_X) ), CV_8UC3 );
        image_to_display = cv::Scalar(255.0, 255.0, 255.0);
        #endif

        #if GRAPH_TYPE == 8
        const double ystep=1.0, thstep=PI/4.0;
        #elif GRAPH_TYPE == 6
        const double ystep=SQRT3BY2, thstep=PI/3.0;
        #endif
        
        #if _VIS
        #if DRAW_INITIAL_GRID
        // all vertices
        cv::Scalar col (255.0, 240.0, 240.0);
        for (double xx=MIN_X; xx<MAX_X; ++xx) 
            for (double yy=MIN_Y; yy<MAX_Y; yy+=ystep) {
                myNode tn(xx, yy); tn.put_in_grid();
                if (!isNodeAccessible(tn)) continue;
                cv::circle (image_to_display, cv_plot_coord(tn.x,tn.y), VERTEX_SIZE*PLOT_SCALE, col, -1, 8);
                for (double th=0.0; th<2*PI-0.0001; th+=thstep)
                    cv::line (image_to_display, cv_plot_coord(tn.x,tn.y), cv_plot_coord(tn.x+1.0*cos(th),tn.y+1.0*sin(th)), 
                            cvScalar(240.0,240.0,240.0), 1);
            }
        #endif
        // obstacles
        for (int a=0; a<OBS_RECT.size(); a++) 
            cv::rectangle(image_to_display, cv_plot_coord(OBS_RECT[a][0],OBS_RECT[a][1]), cv_plot_coord(OBS_RECT[a][2],OBS_RECT[a][3]), 
                        cvScalar(100.0,100.0,100.0), -1 );
        
        cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE);
        cv::imshow("Display window", image_to_display);
        cv::waitKey(0);
        #endif
        
        // Set planner variables
        AllNodesSet.HashTableSize = ceil(MAX_X - MIN_X + 1);
    }
    
    
    bool isNodeAccessible (const myNode& tn) {
        if ( tn.x<MIN_X || tn.x>MAX_X || tn.y<MIN_Y || tn.y>MAX_Y )  return (false);
        bool insideObstacle = false;
        for (int a=0; a<OBS_RECT.size(); ++a)
            if ( OBS_RECT[a][0]<=tn.x && OBS_RECT[a][1]<=tn.y && OBS_RECT[a][2]>=tn.x && OBS_RECT[a][3]>=tn.y ) {
                insideObstacle = true;
                break;
            }
        if (insideObstacle) return (false);
        return (true);
    }
    
    bool isEdgeAccessible (const myNode& tn1, const myNode& tn2) {
        if ( (!isNodeAccessible(tn1)) || !(isNodeAccessible(tn2)) )  return (false);
        double cx=tn2.y-tn1.y, cy=-(tn2.x-tn1.x), cc=(tn2.x*tn1.y-tn1.x*tn2.y);
        for (int a=0; a<OBS_RECT.size(); ++a) {
            if ( (tn1.x<=OBS_RECT[a][0] && tn2.x<=OBS_RECT[a][0])  ||  (tn1.x>=OBS_RECT[a][2] && tn2.x>=OBS_RECT[a][2])  || 
                    (tn1.y<=OBS_RECT[a][1] && tn2.y<=OBS_RECT[a][1])  ||  (tn1.y>=OBS_RECT[a][3] && tn2.y>=OBS_RECT[a][3])  )
                continue; 
            // all corners of the obstacle should lie on the same side of the line
            double last_sign=-3.0, this_sign;
            for (int ix=0; ix<=2; ix+=2)
                for (int iy=1; iy<=3; iy+=2) {
                    this_sign = sign(cx*OBS_RECT[a][ix] + cy*OBS_RECT[a][iy] + cc);
                    if (last_sign>-2.0  &&  last_sign != this_sign)
                        return (false);
                    last_sign = this_sign;
                }
        }
        return (true);
    }
    
    // -----------------------------------------------------------
    
    #if ENABLE_MESH_REFINEMENT && GRAPH_TYPE == 6
    void insertNodeIntoCSvector (myNode& tn, myNode &n, std::vector<myNode>* s, std::vector<double>* c) {
        bool distlessthanone = true;
        double dx, dy, dist;
        if (isEdgeAccessible(tn,n)) {
            s->push_back(tn);
            dx=n.x-tn.x, dy=n.y-tn.y;
            dist = sqrt (dx*dx + dy*dy);
            if (isGreater_d(dist,1.0)) distlessthanone=false;
            c->push_back(dist);
        }
    }
    #endif
    
    void getSuccessors (myNode &n, std::vector<myNode>* s, std::vector<double>* c) 
    {
        // This function should account for obstacles and size of environment.
        
        #if GRAPH_TYPE == 8
        myNode tn;
        
        for (int a=-1; a<=1; ++a)
            for (int b=-1; b<=1; ++b) {
                int xParity = ((int)round(fabs(n.x))) % 2;
                if (a==0 && b==0) continue;
                if (xParity==0 && (a!=0 && b==-1)) continue;
                if (xParity==1 && (a!=0 && b==1)) continue;
                
                tn.x = n.x + a;
                tn.y = n.y + b;
                
                #if GRAPH8_X_WIGGLE
                tn.x = round (tn.x);
                if (((int)round(fabs(tn.y)))%2==1) 
                    tn.x += GRAPH8_X_WIGGLE_AMOUNT;
                #endif
                
                //if (!isNodeAccessible(tn)) continue;
                if (!isEdgeAccessible(tn,n)) continue;
                
                s->push_back(tn);
                double dx=tn.x-n.x, dy=tn.y-n.y;
                c->push_back(sqrt(dx*dx+dy*dy)); 
            }
        
        #elif GRAPH_TYPE == 6
        myNode tin, in;
        tin.x=n.x; tin.y=n.y; in.x=n.x; in.y=n.y;
        
        #if ENABLE_MESH_REFINEMENT
        std::vector<myNode> ins, tiins;
        double nHexDist=0.0, tinHexDist;
        
        int refinerIdx = 0;
        // check if n is inside or on bounday of a disk. If yes, invert it.
        while (refinerIdx < refiners.size()) {
            ins = refiners[refinerIdx].get_inverted_pt (n, &nHexDist);
            if ( isLess_d (nHexDist, refiners[refinerIdx].R) ) { // 'n' is strictly inside the disk -- work with inverted point
                in = ins[0];
                break;
            }
            else if ( isEqual_d (nHexDist, refiners[refinerIdx].R) ) { // 'n' is on boundary of a disk
                break;
            }
            ++refinerIdx;
        }
        #endif
        
        double th;
        for (int a=0; a<6; ++a) {
            th = a * PI_BY_3;
            tin.x = in.x + 1.0*cos(th);
            tin.y = in.y + 1.0*sin(th);
            
            #if ENABLE_MESH_REFINEMENT
            
            if (refinerIdx < refiners.size()) {
                
                tiins = refiners[refinerIdx].get_inverted_pt (tin, &tinHexDist);
                
                // ---------------------------------------------------
                if ( isEqual_d (nHexDist, refiners[refinerIdx].R) ) { // 'n' was on outer boundary
                    // NOTE: n==in is a outer boundary point. tin==tn is a neighbor.
                    
                    if ( isLess_d (tinHexDist, refiners[refinerIdx].R) )  // the neighbor, tin==tn, is inside outer boundary.
                        continue; // do not insert this
                    else { // the neighbor, tin==tn, is outside outr boundary. Will insert 'tin' and its inverse
                        insertNodeIntoCSvector (tin, n, s, c);
                        // chek if 'tin' itsef is not on boundary
                        if (!isEqual_d (tinHexDist, refiners[refinerIdx].R)  ) 
                            insertNodeIntoCSvector (tiins[0], n, s, c); // g1 branch
                        continue;
                    }
                } // ---------------------------------------------------
                
                else if ( isLess_d (nHexDist, refiners[refinerIdx].R) ) { // 'n' was strictly inside outer boundary
                    //NOTE: 'in' was the inverse of 'n' (and distinct from it). Thus, tin != tn.
                    //      Need to use tiin (inverse of 'tin')
                    
                    if ( isGreater_d (nHexDist, refiners[refinerIdx].rThresh) )  { // 'n' is strictly between the inner and outer boundary
                        insertNodeIntoCSvector (tiins[0], n, s, c); // branch g1, function f (since nHexDist < R)
                        continue;
                    }
                    else if ( isEqual_d (nHexDist, refiners[refinerIdx].rThresh) ) { // 'n' is on the inner boundary
                        
                        if ( isGreater_d (tinHexDist, refiners[refinerIdx].rInvThresh) ) // 'tin' is outside inverted inner disk
                            continue; // dont insert 'tiin'
                        else { // insert both inverses of 'tin'
                            insertNodeIntoCSvector (tiins[1], n, s, c); // branch g2
                            // Check if 'tin' is itself not on the boundary of inverted inner disk
                            if (!isEqual_d (tinHexDist, refiners[refinerIdx].rInvThresh) ) 
                                insertNodeIntoCSvector (tiins[0], n, s, c); // branch g1, function f (since nHexDist < R)
                            continue;
                        }
                    }
                    else { // 'n' is inside inner disk boundary
                        insertNodeIntoCSvector (tiins[1], n, s, c); // branch g2
                        continue;
                    }
                        
                }
                // ---------------------------------------------------
                // else ...
            }
            
            
            #endif
            
            // in == n   &&   tin == tn
            if (!isEdgeAccessible(tin,in)) continue;
            s->push_back(tin);
            c->push_back(1.0);
        }
        
        #endif
    }
    
    // -----------------------------------------------------------
    
    double getHeuristics (myNode& n)
    {
        // printf("in GetHeuristics: %f\n", sqrt((double)(n.x*n.x + n.y*n.y)));
        /* double dx = goalNode.x - n.x;
        double dy = goalNode.y - n.y;
        return (sqrt(dx*dx + dy*dy)); */
        return (0.0);
    }
    
    // -----------------------------------------------------------
    
    std::vector<myNode> getStartNodes (void) 
    {
        std::vector<myNode> startNodes;
        myNode tn;
        
        printf("boundaries: %f, %f, %f, %f\n", MAX_X, MIN_X, MAX_Y, MIN_Y);
        
        tn = myNode(0.0, 0.0);
        tn.put_in_grid();
        if (isNodeAccessible(tn)) { 
            startNodes.push_back (tn);
            #if _VIS
            cv::circle (image_to_display, cv_plot_coord(tn.x,tn.y), VERTEX_SIZE*PLOT_SCALE, getExpandedNodeeColor(tn), -1, 8);
            #endif
        }
        
        return (startNodes);
    }
    
    
    // =================================================================
    
    CvScalar getExpandedNodeeColor (myNode &n) {
        double gScale = 1.0; // exp(-0.01*n.G); //1.0;
        double r=0.0, g=0.0, b=0.0;
        while (abs(r-b)<50 && g-b<100) {
            r = 60.0 + (rand()%151);
            g = 60.0 + (rand()%149);
            b = 60.0 + (rand()%157);
        } 
        r *= 0.5; //1.0 + sin(20*n.G/(MAX_X-MIN_X)); // / (MAX_X-MIN_X);
        return (cvScalar(b, g, r));
    }
    
    // -------------------------------------------
    
    
    
    void nodeEvent (myNode &n, unsigned int e) 
    {
        CvScalar col = cvScalar(0.0, 0.0, 0.0);;
        int thickness = -1; //LINE_THICKNESS;
        double radFactor = 1.0;
        
        bool pauseForVis=false, drawVertex=true;;
        
        // --------------------------------------------
        if (e & EXPANDED) {
            // n.print ("expanded");
            lastExpanded = n;
            bool cameFromNull = false;
            #if defined(S_STAR_ALGORITHM)
            cameFromNull = (n.CameFromSimplex==NULL);
            #endif
            if (!cameFromNull) {
                
                if ( !isEqual_d (n.G, n.dist()) ) {
                    //printf("Distance mismatch! Expanded:%d, event:%d-%d-%d; dist=%f\n", self_p->Expanded, e,EXPANDED,(e|EXPANDED), fabs(self_p->G - dist));
                    #if defined(S_STAR_ALGORITHM) && _YAGSBPL_DEBUG>1
                    printf("*** Distance mismatch! diff=%f\n", fabs(n.G - n.dist()));
                    n.CameFromSimplex->print ("Came-from simplex: ");
                    pauseForVis = true;
                    #endif
                    #if VERTEX_COLORS
                    double red = MIN(255.0, MIN(255.0, 10*255.0*fabs(n.G - n.dist())/n.dist() ) );
                    col = cvScalar(0.0, 255.0-red, red); // red
                    #else
                    col = cvScalar(200.0, 200.0, 200.0);
                    #endif
                }
                else
                    #if VERTEX_COLORS
                    col = cvScalar(0.0, 255.0, 0.0);
                    #else
                    col = cvScalar(200.0, 200.0, 200.0);
                    #endif
            }
            else {
                printf ("Expanded, but came-from is NULL!!\n");
                pauseForVis = true;
            }
        }
        
        else if ((e & HEAP) == PUSHED) {
            col = cvScalar(255.0, 0.0, 0.0); // blue
        }
        
        else if (e & UNEXPANDED) {
            col = cvScalar(255.0, 0.0, 150.0); // purple
            if ( !isEqual_d (n.G, n.dist()) ) {
                #if defined(S_STAR_ALGORITHM) && _YAGSBPL_DEBUG>1
                printf("*** Distance mismatch even with backtrace correction backtrace! diff=%f\n", fabs(n.G - n.dist()));
                n.CameFromSimplex->print ("Came-from simplex: ");
                pauseForVis = true;
                #endif
                col = cvScalar(255.0, 0.0, 255.0); // red-purple
            }
            pauseForVis = true;
        }
        
        else
            return;
                
        //-------------------------------------------
        #if _VIS
        if (drawVertex) {
            double nodeRad = radFactor*VERTEX_SIZE;
            cv::circle (image_to_display, cv_plot_coord(n.x,n.y), nodeRad*PLOT_SCALE, col, thickness, 8);
            for (auto its=n.Successors.begin(); its!=n.Successors.end(); ++its) {
                double Dx = (*its).first->x - n.x, Dy = (*its).first->y - n.y;
                double Dd = sqrt(Dx*Dx + Dy*Dy);
                double nx=Dx/Dd, ny=Dy/Dd;
                cv::line(image_to_display, 
                    cv_plot_coord(n.x+nodeRad*nx, n.y+nodeRad*ny), 
                    cv_plot_coord((*its).first->x-nodeRad*nx, (*its).first->y-nodeRad*ny), 
                        cvScalar(200.0,200.0,200.0), 0.1*PLOT_SCALE);
            }
        }
        
        if (ExpandCount%VIS_INTERVAL == 0  ||  NodeHeap.size() == 0) {
            cv::imshow("Display window", image_to_display);
            std::cout << std::flush;
            if (SAVE_IMG_INTERVAL>0 && ExpandCount%SAVE_IMG_INTERVAL == 0) {
                char imgFname[1024];
                sprintf(imgFname, "outfiles/%s%05d.png", imgPrefix.str().c_str(), ExpandCount);
                cv::imwrite (imgFname, image_to_display);
            }
            cvWaitKey(1); //(10);
        }
        
        if (pauseForVis) {
            // cvWaitKey();
        }
        #endif
    }
    
    // ---------------------------------------
    
    bool stopSearch (myNode &n) {
        // return (false);
        // return (sqrt(n.x*n.x +n.y*n.y) > 30);
        return (n==goalNode);
    }
};

// ==============================================================================

int main(int argc, char *argv[])
{
    SearchProblem test_search_problem;
    test_search_problem.goalNode = myNode(250*test_search_problem.PROBLEM_SCALE, -150*test_search_problem.PROBLEM_SCALE);
    // Run
    test_search_problem.search();
    
    // Get path
    auto path = test_search_problem.reconstructPath (test_search_problem.goalNode);
    
    printf("\nPath: ");
    #if VERTEX_COLORS
    CvScalar path_color = cvScalar(250.0,100.0,100.0);
    #else
    CvScalar path_color = cvScalar(200.0,100.0,250.0);
    #endif
    myNode thisPt= myNode(0.0,0.0), lastPt;
    std::vector<myNode> allPts;
    for (int a=path.size()-1; a>=0; --a) {
        lastPt = thisPt;
        thisPt = myNode(0.0, 0.0);
        for (auto it=path[a].begin(); it!=path[a].end(); ++it) {
            thisPt.x += it->second * it->first->x;
            thisPt.y += it->second * it->first->y;
        }
        allPts.push_back (thisPt);
        printf ("[%f,%f]; ", thisPt.x, thisPt.y);
        #if _VIS
        cv::line (test_search_problem.image_to_display, 
                test_search_problem.cv_plot_coord(thisPt.x,thisPt.y), test_search_problem.cv_plot_coord(lastPt.x,lastPt.y), 
                            path_color, 2);
        #endif
    }
    #if _VIS
    for (int a=allPts.size()-1; a>=0; --a)
        cv::circle (test_search_problem.image_to_display, 
                test_search_problem.cv_plot_coord(allPts[a].x,allPts[a].y),
                    test_search_problem.VERTEX_SIZE*test_search_problem.PLOT_SCALE*0.5, cvScalar(255.0,0.0,0.0), -1, 8);
    #endif
    
    printf("\n");
    cv::imshow("Display window", test_search_problem.image_to_display);
    
    test_search_problem.clear();
    
    if (SAVE_IMG_INTERVAL != 0) {
        char imgFname[1024];
        sprintf(imgFname, "outfiles/%s_path.png", test_search_problem.imgPrefix.str().c_str());
        cv::imwrite(imgFname, test_search_problem.image_to_display);
    }
    cvWaitKey();
}

