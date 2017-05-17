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

#ifndef _DOSL_ALGORITHM // can pass at command line during compilation: -D_DOSL_ALGORITHM=AStar
    #define _DOSL_ALGORITHM  AStar
#endif
#include <dosl/dosl>

// Local libraries/headers
#include "cvParseMap2d.h"
#include "JSONparser.tcc"

// =======================
// Other parameters

#define GRAPH_TYPE 6 // 6 or 8
#define PROB_HMTPY 1

// math macros / constants
#define sign(x) ((x>0.0)?1.0:((x<0.0)?-1.0:0.0))
#define PI       3.14159265359
#define PI_BY_3  1.0471975512
#define SQRT3BY2 0.86602540378
#define INFINITESIMAL_DOUBLE  1e-6

// output options
#define _STAT 0

#define _VIS 1
#define VERTEX_COLORS 1

#define SAVE_IMAGE -1
#define SAVE_IMG_INTERVAL 1000 // 1 //10000

// ---------------------------------------------------
// derived values
#if GRAPH_TYPE==8
    #define COORD_TYPE int
#else
    #define COORD_TYPE double
#endif

// ==============================================================================
// helper functions

int approx_floor (double x, double tol=INFINITESIMAL_DOUBLE) {
    int ret = floor (x);
    if (ret+1-x<tol) ++ret;
    return (ret);
}

double wiggle = 1e-3;

// ==============================================================================

// A node of the graph
class myNode : public DOSL_CLASS(Node)<myNode,double>
{
public:
    COORD_TYPE x, y;
    std::vector<int> h;
    
    #if GRAPH_TYPE == 8
        void put_in_grid (void) { }
    #elif GRAPH_TYPE == 6
        void put_in_grid (void) {
            int yLevel = round (y / SQRT3BY2);
            y = SQRT3BY2*yLevel;
            if (yLevel%2 == 0) // even
                x = round(x+INFINITESIMAL_DOUBLE);
            else
                x = round(x+0.5+INFINITESIMAL_DOUBLE) - 0.5;
        }
    #endif
    
    bool isCoordsEqual(const myNode& n) const {
        #if GRAPH_TYPE == 8 // COORD_TYPE int
        return ((x==n.x) && (y==n.y));
        #else
        return (fabs(x-n.x)<INFINITESIMAL_DOUBLE  &&  fabs(y-n.y)<INFINITESIMAL_DOUBLE);
        #endif
    }
    
    // *** This must be defined for the node
    bool operator==(const myNode& n) const {
        if (!isCoordsEqual(n)) return (false);
        if (h.size()!=n.h.size()) return (false);
	    for (int a=0; a<h.size(); a++)
	        if (h[a]!=n.h[a]) return (false);
        return (true);
    }
    
    // constructor
    myNode () { }
    myNode (COORD_TYPE xx, COORD_TYPE yy) : x(xx), y(yy) { put_in_grid(); }
    
    
    // Inherited functions being overwritten
    int getHashBin (void) const {
        #if GRAPH_TYPE == 8 // COORD_TYPE int
        return (abs(x));
        #else
        return ( MAX(round(fabs(x)+INFINITESIMAL_DOUBLE), round(fabs(x)-INFINITESIMAL_DOUBLE)) );
        #endif
    }
    
    // print
    void print (std::string head="", std::string tail="\n") const {
        _dosl_cout << _GREEN + head << " (" << this << ")" GREEN_ " x=" << x << ", y=" << y << "; ";
        (G==std::numeric_limits<double>::max())? printf("INF") : printf("G = %0.8f; h = [", G);
        for (int a=0; a<h.size(); ++a){
            printf("%d", a);
            if (a!=h.size()-1) printf(", ");
        } 
        printf("]\n");
        _dosl_cout << "Successors: ";
        for (auto it=Successors.begin(); it!=Successors.end(); ++it)
            printf ("%x (%f), ", it->first, it->second);
        std::cout << tail << _dosl_endl;
    }
};

// ==============================================================================

class searchProblem : public DOSL_CLASS(Problem)<myNode,double>
{
public:
    // Fime names and JSON objects
    std::string   map_image_fName, expt_fName, expt_folderName, expt_Name;
    JSONcontainer expt_container;
    cvParseMap2d  my_map;
    int nClassesToFind;
    
    // Image display variables / parameters
    cv::Mat image_to_display;
    double PLOT_SCALE;
    double VERTEX_SIZE, LINE_THICKNESS;
    
    // variables for saving image
    int frameno;
    std::ostringstream imgPrefix;
    
    // variables decsribing problem
    COORD_TYPE MAX_X, MIN_X, MAX_Y, MIN_Y;
    myNode startNode, goalNode, lastExpanded;
    
    // homotopy classes
    int nClasses;
    std::vector<myNode> homotopyGoals;
    
    // -----------------------------------------------------------
    
    template<typename T>
    CvPoint cv_plot_coord(T x, T y) { return( cvPoint( round(PLOT_SCALE*(x-MIN_X)), round(PLOT_SCALE*(y-MIN_Y)) ) ); }
    
    void cvPlotPoint (CvPoint pt, CvScalar colr, int size=1) {
        for (int i=pt.x; i<=pt.x+(size-1); i++)
            for (int j=pt.y; j<=pt.y+(size-1); j++) {
                if (i<0 || i>=image_to_display.cols || j<0 || j>=image_to_display.rows) continue;
                image_to_display.at<cv::Vec3b>(j,i) = cv::Vec3b ((uchar)colr.val[0], (uchar)colr.val[1], (uchar)colr.val[2]);
            }
    }
    
    // Constructor
    searchProblem (std::string expt_f_name, std::string expt_name)
    {
        expt_fName = expt_f_name; expt_Name = expt_name;
        expt_folderName = expt_fName.substr(0, expt_fName.find_last_of("/\\")+1);
        
        // Read from file
        std::ifstream my_fstream (expt_fName);
        expt_container = JSONcontainer (my_fstream)[expt_Name];
        map_image_fName = expt_folderName + expt_container["map_name"].as<std::string>();
        my_map = cvParseMap2d (map_image_fName, true); // computes representative points
        
        // read data for planning
        MAX_X=my_map.width(); MIN_X=0; MAX_Y=my_map.height(); MIN_Y=0;
        startNode = myNode (expt_container["start"][0].as<int>(), expt_container["start"][1].as<int>());
        goalNode = myNode (expt_container["goal"][0].as<int>(), expt_container["goal"][1].as<int>());
        startNode.put_in_grid(); goalNode.put_in_grid(); 
        nClassesToFind = expt_container["top_class"].as<int>();
        nClasses = 0;
        
        // display options
        PLOT_SCALE = 2.0;
        VERTEX_SIZE = 1.0;
        LINE_THICKNESS = 4.0; // CV_FILLED
        
        // saving options
        frameno = 0;
        imgPrefix << MAKESTR(ALGORITHM) << GRAPH_TYPE << "_";
        
        #if _VIS
        image_to_display = my_map.getCvMat (COLOR_MAP);
        cv::resize (image_to_display, image_to_display, cv::Size(), PLOT_SCALE , PLOT_SCALE );
        cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE);
        cv::imshow("Display window", image_to_display);
        cv::waitKey(0);
        #endif
        
        // Set planner variables
        AllNodesSet.HashTableSize = ceil(MAX_X - MIN_X + 1);
    }
    
    // -----------------------------------------------------------
    
    bool isNodeInWorkspace (const myNode& tn) {
        if ( tn.x<MIN_X || tn.x>MAX_X || tn.y<MIN_Y || tn.y>MAX_Y )  return (false);
        return (true);
    }
    
    bool isEdgeAccessible (const myNode& tn1, const myNode& tn2) {
        if ( (!isNodeInWorkspace(tn1)) || !(isNodeInWorkspace(tn2)) )  return (false);
        
        #if GRAPH_TYPE == 8 // COORD_TYPE int
        // the following works only for 8-connected grid!!
        if (tn1.x!=tn2.x && tn1.y!=tn2.y) // diagonal edge
            return ( my_map.isFree ( round(MIN(tn1.x,tn2.x)), round(MIN(tn1.y,tn2.y)) ) );
        else if (tn1.x==tn2.x) // need to check two cells
            return (!( (tn1.x==MAX_X || my_map.isObstacle ( round(tn1.x), round(MIN(tn1.y,tn2.y)) ) ) &
                        ( tn1.x==MIN_X || my_map.isObstacle ( round(tn1.x-1), round(MIN(tn1.y,tn2.y)) ) ) ));
        else if (tn1.y==tn2.y) // need to check two cells
            return (!( ( tn1.y==MAX_Y || my_map.isObstacle ( round(MIN(tn1.x,tn2.x)), round(tn1.y) ) ) &
                        ( tn1.y==MIN_Y || my_map.isObstacle ( round(MIN(tn1.x,tn2.x)), round(tn1.y-1) ) ) ));
        #else
        // TODO: make more precise
        return ( my_map.isFree(approx_floor(tn1.x),approx_floor(tn1.y))  && 
                    my_map.isFree(approx_floor(tn2.x),approx_floor(tn2.y)) );
        #endif
        
        return (true);
    }
    
    // -----------------------------------------------------------
    
    void updateHSignature (myNode &n, myNode &tn) { // updates h/H-signature of tn
        int pm;
        tn.h = n.h;
        if (n.x!=tn.x) {
            int p_start, p_end, p_step;
            if (n.x<tn.x) { p_start = 0; p_end = my_map.repPts.size()-1; p_step = 1; }
            else { p_start = my_map.repPts.size()-1; p_end = 0; p_step = -1; }
        
            for (int p=p_start; p!=p_end+p_step; p+=p_step) {
                std::vector<double> this_rep_pt = { my_map.repPts[p].x+(p+1)*wiggle, my_map.repPts[p].y+(p+1)*wiggle };
	            pm = 0;
	            if (n.y>this_rep_pt[1] && tn.y>this_rep_pt[1]) {
	                if (n.x<=this_rep_pt[0] && tn.x>this_rep_pt[0]) 
	                    pm = 1;
	                else if (n.x>this_rep_pt[0] && tn.x<=this_rep_pt[0]) 
	                    pm = -1;
	            }
	            if (pm)
	            {
                    #if	PROB_HMTPY
                    if (tn.h.size()>0 && tn.h[tn.h.size()-1]==-pm*(p+1))
                        tn.h.pop_back();
                    else
                        tn.h.push_back(pm*(p+1));
                    #else // homology
                    tn.h[p] += pm;
                    #endif
	            }
            }
        }
    }
    
    void getSuccessors (myNode &n, std::vector<myNode>* s, std::vector<double>* c) // *** This must be defined
    {
        // This function should account for obstacles and size of environment.
        myNode tn;
        
        #if GRAPH_TYPE == 8
        for (int a=-1; a<=1; ++a)
            for (int b=-1; b<=1; ++b) {
                if (a==0 && b==0) continue;
                
                #ifdef DOSL_ALGORITHM_SStar
                int xParity = ((int)round(fabs(n.x))) % 2;
                if (xParity==0 && (a!=0 && b==-1)) continue;
                if (xParity==1 && (a!=0 && b==1)) continue;
                #endif
                
                tn.x = n.x + a;
                tn.y = n.y + b;
                
                if (!isEdgeAccessible(tn,n)) continue;
                
                updateHSignature(n,tn);
                
                s->push_back(tn);
                double dx=tn.x-n.x, dy=tn.y-n.y;
                c->push_back(sqrt(dx*dx+dy*dy)); 
            }
        
        #elif GRAPH_TYPE == 6
        double th;
        for (int a=0; a<6; ++a) {
            th = a * PI_BY_3;
            tn.x = n.x + 1.0*cos(th);
            tn.y = n.y + 1.0*sin(th);
            
            if (!isEdgeAccessible(tn,n)) continue;
            
            updateHSignature(n,tn);
                
            s->push_back(tn);
            c->push_back(1.0);
        }
        
        #endif
        
    }
    
    // -----------------------------------------------------------
    
    double getHeuristics (myNode& n)
    {
        /* double dx = goalNode.x - n.x;
        double dy = goalNode.y - n.y;
        return (sqrt(dx*dx + dy*dy)); */
        return (0.0);
    }
    
    // -----------------------------------------------------------
    
    std::vector<myNode> getStartNodes (void) 
    {
        std::vector<myNode> startNodes;
        
        startNodes.push_back (startNode);
        #if _VIS
        cv::circle (image_to_display, cv_plot_coord(startNode.x,startNode.y), VERTEX_SIZE*PLOT_SCALE, 
                                                                cvScalar (200.0, 150.0, 150.0), -1, 8);
        #endif
        
        return (startNodes);
    }
    
    // -----------------------------------------------------------
    
    void nodeEvent (myNode &n, unsigned int e) 
    {
        #if _VIS
        
        CvScalar col = cvScalar(0.0, 0.0, 0.0);;
        int thickness = -1; //lineThickness;
        double radFactor = 1.0;
        
        bool pauseForVis=false, drawVertex=true;
        
        // --------------------------------------------
        if (e & EXPANDED) {
            // ++++++++++++++++++++++++++++++++++++++++++++++++++++
            // Testing specific vertices
            /* if (fabs(n.x-172.0)<0.01 && fabs(n.y-98.7269)<0.01)
                RUNTIME_VERBOSE_SWITCH = 1; */
            // ++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            lastExpanded = n;
            bool cameFromNull = false;
            #ifdef DOSL_ALGORITHM_SStar
            cameFromNull = (n.CameFromSimplex==NULL);
            #endif
            if (!cameFromNull) {
                #if VERTEX_COLORS
                col = cvScalar(150.0, 255.0, 150.0);
                #else
                col = cvScalar(255.0, 255.0, 255.0);
                #endif
            }
            else {
                printf ("Expanded, but came-from is NULL!!\n");
                // pauseForVis = true;
            }
        }
        
        else if ((e & HEAP) == PUSHED) {
            #if VERTEX_COLORS
            col = cvScalar(255.0, 0.0, 0.0); // blue
            #else
            col = cvScalar(200.0, 200.0, 200.0);
            #endif
        }
        
        else if (e & UNEXPANDED) {
            col = cvScalar(255.0, 0.0, 150.0); // purple
            printf ("Backtrack started!!\n");
        }
        
        else
            return;
                
        //-------------------------------------------
        if (drawVertex) {
            double nodeRad = radFactor*VERTEX_SIZE;
            if (n.x<MAX_X  &&  n.y<MAX_Y  &&  my_map.isFree(round(n.x), round(n.y)) )
                cvPlotPoint (cv_plot_coord(n.x,n.y), col, PLOT_SCALE);
        }
        
        if (ExpandCount%SAVE_IMG_INTERVAL == 0  ||  NodeHeap.size() == 0) {
            cv::imshow("Display window", image_to_display);
            std::cout << std::flush;
            #if SAVE_IMAGE>0
            char imgFname[1024];
            sprintf(imgFname, "outfiles/%s%05d.png", imgPrefix.str().c_str(), ExpandCount);
            cvSaveImage(imgFname, ipl_image_p);
            #endif
            cvWaitKey(1); //(10);
        }
        
        if (pauseForVis) {
            cvWaitKey();
        }
        #endif
    }
    
    // ---------------------------------------
    
    bool stopsearch (myNode &n) {
        if (n.isCoordsEqual(goalNode)) {
            homotopyGoals.push_back (n);
            ++nClasses;
            n.print ("Found a path to ");
            if (nClasses>=nClassesToFind)
                return (true);
        }
        return (false);
    }
};

// ==============================================================================

int main(int argc, char *argv[])
{
    //RUNTIME_VERBOSE_SWITCH = 0;
    
    std::string expt_f_name="exptfiles/experiments.json", expt_name="path_plan_1";
    if (argc > 1) {
        expt_f_name = argv[1];
        expt_name = argv[2];
    }
    
    searchProblem test_search_problem (expt_f_name, expt_name);
    test_search_problem.search();
    
    // -------------------------------------
    #if _STAT
    char statFname[1024];
    sprintf(statFname, "outfiles/%s_%s.txt", test_search_problem.imgPrefix.str().c_str(), test_search_problem.map_image_fName.c_str());
    FILE* pFile;
    pFile = fopen (statFname,"a+");
    if (pFile!=NULL)
        fprintf (pFile, "%s:\n", test_search_problem.expt_Name.c_str());
    #endif
    
    for (int i=0; i<test_search_problem.homotopyGoals.size(); ++i) {
        // get path
        auto path = test_search_problem.reconstructPath (test_search_problem.homotopyGoals[i]);
        double cost = 0.0;
        
        // Print and draw path
        #if _VIS
        //printf("\nPath: ");
        #endif
        myNode thisPt=test_search_problem.startNode, lastPt;
        std::vector<myNode> allPts;
        for (int a=path.size()-1; a>=0; --a) {
            lastPt = thisPt;
            thisPt = myNode(0.0, 0.0);
            for (auto it=path[a].begin(); it!=path[a].end(); ++it) {
                thisPt.x += it->second * it->first->x;
                thisPt.y += it->second * it->first->y;
            }
            allPts.push_back (thisPt);
            #if _VIS
            //printf ("[%f,%f]; ", (double)thisPt.x, (double)thisPt.y);
            cv::line (test_search_problem.image_to_display, 
                    test_search_problem.cv_plot_coord(thisPt.x,thisPt.y), test_search_problem.cv_plot_coord(lastPt.x,lastPt.y), 
                                cvScalar(200.0,100.0,100.0),
                                        test_search_problem.LINE_THICKNESS*test_search_problem.PLOT_SCALE);
            #endif
            cost += sqrt ((thisPt.x-lastPt.x)*(thisPt.x-lastPt.x) + (thisPt.y-lastPt.y)*(thisPt.y-lastPt.y));
        }
        
        #if _VIS
        //printf ("\n");
        for (int a=allPts.size()-1; a>=0; --a)
            cv::circle (test_search_problem.image_to_display, 
                    test_search_problem.cv_plot_coord(allPts[a].x,allPts[a].y),
                        test_search_problem.VERTEX_SIZE*test_search_problem.PLOT_SCALE, cvScalar(150.0,0.0,255.0), -1, 8);
        #endif
        
        printf ("\tclass: %d, cost: %f.\n", i, cost);
        #if _STAT
        if (pFile!=NULL)
            fprintf (pFile, "\tclass: %d, cost: %f.\n", i, cost);
        #endif
    }
    
    #if _VIS
    cv::imshow("Display window", test_search_problem.image_to_display);
    #endif
    
    #if _STAT
    if (pFile!=NULL)
        fclose (pFile);
    #endif
    
    test_search_problem.clear();
    
    #if _VIS
    #if SAVE_IMAGE!=0
    char imgFname[1024];
    sprintf(imgFname, "outfiles/%s_%s_%s_path.png", test_search_problem.imgPrefix.str().c_str(), 
                                                    test_search_problem.map_image_fName.c_str(), test_search_problem.expt_Name.c_str());
    cv::imwrite(imgFname, test_search_problem.image_to_display);
    #endif
    //printf ("\ncomputation time = %f\n", 0.0);
    cvWaitKey();
    #endif
}

