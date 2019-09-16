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
#include <fstream>
#include <vector>
// Open CV:
#include <opencv2/opencv.hpp> 
#include <opencv2/highgui.hpp> 

// DOSL library
#include <dosl/dosl>

// Local utility libraries/headers
#include <dosl/aux-utils/cvParseMap2d.hpp>
#include <dosl/aux-utils/double_utils.hpp>
#include <dosl/aux-utils/string_utils.hpp> // compute_program_path
#include "../../include-local/RSJparser.tcc"

// =======================

// Algorithm
#ifndef _DOSL_ALGORITHM // can pass at command line during compilation: -D_DOSL_ALGORITHM=AStar
    #define _DOSL_ALGORITHM  AStar // AStar, SStar or ThetaStar
#endif

// Graph/planning options
#define GRAPH_TYPE 8 // 6 or 8

// display options
#define _STAT 0
#define _VIS 1
#define VIS_INTERVAL 100
#define VERTEX_COLORS 0
#define SAVE_IMG_INTERVAL -1 // 0 to not save at all. -1 to save last frame only.


// ==============================================================================

#define COORD_TYPE double // can be 'int' only if (GRAPH_TYPE==8 && _DOSL_ALGORITHM!=SStar)

// A node of the graph
class myNode : public _DOSL_ALGORITHM::Node<myNode,double>
{
public:
    COORD_TYPE x, y;
    
    // Comparison operator:
    bool operator==(const myNode& n) const {
        #if GRAPH_TYPE == 8 // COORD_TYPE int
        return ((x==n.x) && (y==n.y));
        #else
        return (fabs(x-n.x)<INFINITESIMAL_DOUBLE  &&  fabs(y-n.y)<INFINITESIMAL_DOUBLE);
        #endif
    }
    
    // Inherited functions being overwritten
    int getHashBin (void) const {
        #if GRAPH_TYPE == 8 // COORD_TYPE int
        return (abs(x));
        #else
        return ( MAX(round(fabs(x)+INFINITESIMAL_DOUBLE), round(fabs(x)-INFINITESIMAL_DOUBLE)) );
        #endif
    }
    
    // --------------------------------------------
    
    // constructor
    myNode () { }
    myNode (COORD_TYPE xx, COORD_TYPE yy) : x(xx), y(yy) { put_in_grid(); }
    
    #if GRAPH_TYPE == 8
        void put_in_grid (void) { x = (COORD_TYPE)round(x); y = (COORD_TYPE)round(y); }
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
    
    // print
    void print (std::string head="", std::string tail="") const {
        _dosl_cout << _GREEN + head << " (" << this << ")" GREEN_ " x=" << x << ", y=" << y << ", lin=" << lineage_data.id << ", ";
        (g_score==std::numeric_limits<double>::max())? printf("INF") : printf("%0.8f", g_score);
        /*printf ("\n");
        _dosl_cout << "\tExpanded=%d" << expanded << _dosl_endl;
        _dosl_cout << "successors: ";
        for (auto it=successors.begin(); it!=successors.end(); ++it)
            printf ("%x (%f), ", it->first, it->second); */
        std::cout << tail << _dosl_endl;
    }

};


class UnorderedIndexPair {
public:
    int i, j;
    UnorderedIndexPair (int a, int b) : i(a), j(b) { }
    bool operator==(const UnorderedIndexPair& p) const {
        return ((i==p.i && j==p.j) || (i==p.j && j==p.i));
    }
};

// ==============================================================================

class searchProblem : public _DOSL_ALGORITHM::Algorithm<searchProblem,myNode,double>
{
public:
    // Fime names and JSON objects
    std::string   map_image_fName, expt_fName, expt_folderName, expt_Name, out_folderName;
    cvParseMap2d  my_map;
    RSJresource expt_container;
    
    // Image display variables / parameters
    cv::Mat image_to_display;
    double PLOT_SCALE;
    double VERTEX_SIZE, LINE_THICKNESS;
    
    // variables for saving image
    int frameno;
    std::ostringstream imgPrefix;
    
    // variables decsribing problem
    COORD_TYPE MAX_X, MIN_X, MAX_Y, MIN_Y;
    myNode lastExpanded;
    std::vector<myNode> startNodes;
    std::vector<CvScalar> colors;
    
    
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
    searchProblem (std::string expt_f_name, std::string expt_name, std::string out_folder_name)
    {
        expt_fName = expt_f_name; expt_Name = expt_name;
        out_folderName = out_folder_name;
        expt_folderName = expt_fName.substr(0, expt_fName.find_last_of("/\\")+1);
        
        // Read from file
        std::ifstream my_fstream (expt_fName);
        expt_container = RSJresource (my_fstream)[expt_Name];
        
        // obstacle map
        if (expt_container["environment"]["pixmap"].exists()) {
            map_image_fName = expt_folderName + expt_container["environment"]["pixmap"].as<std::string>();
            my_map = cvParseMap2d (map_image_fName, false);
        }
        else if (expt_container["environment"]["width"].exists() && expt_container["environment"]["height"].exists()) {
            my_map = cvParseMap2d (cv::Mat (expt_container["environment"]["height"].as<int>(), 
                                            expt_container["environment"]["width"].as<int>(),
                                                            CV_8UC3, cv::Scalar(255,255,255)), false);
        }
        
        // read data for planning
        MAX_X=my_map.width(); MIN_X=0; MAX_Y=my_map.height(); MIN_Y=0;
        
        if (expt_container["start"].type() == RSJ_ARRAY)
            for (int a=0; a<expt_container["start"].size(); ++a) {
                startNodes.push_back (myNode (expt_container["start"][a][0].as<COORD_TYPE>(), expt_container["start"][a][1].as<COORD_TYPE>()));
                startNodes[a].put_in_grid();
                startNodes[a].print("Start node: ");
                colors.push_back (cvScalar (rand()%256, rand()%256, rand()%256));
            }
        
        else if (expt_container["start"].type() == RSJ_OBJECT) {
            
            if (expt_container["start"]["list"].exists())
                for (int a=0; a<expt_container["start"]["list"].size(); ++a) {
                    startNodes.push_back (myNode (expt_container["start"]["list"][a][0].as<COORD_TYPE>(), 
                                                    expt_container["start"]["list"][a][1].as<COORD_TYPE>()));
                    startNodes[a].put_in_grid();
                    startNodes[a].print("Start node: ");
                    colors.push_back (cvScalar (rand()%256, rand()%256, rand()%256));
                }
            
            if (expt_container["start"]["random"].exists()) {
                int tx, ty;
                for (int a=0; a<expt_container["start"]["random"].as<int>(1); ++a) {
                    do {
                        tx = rand() % int(MAX_X);
                        ty = rand() % int(MAX_Y);
                    } while (my_map.isObstacle (tx, ty));
                    startNodes.push_back (myNode (COORD_TYPE(tx), COORD_TYPE(ty)) );
                    startNodes[a].put_in_grid();
                    startNodes[a].print("Start node: ");
                    colors.push_back (cvScalar (rand()%256, rand()%256, rand()%256));
                }
            }
            
            if (expt_container["start"]["random_triangulate"].exists()) {
                int u,v;
                COORD_TYPE tx, ty;
                int sep = expt_container["start"]["random_triangulate"]["minsep"].as<int>(1);
                myNode tn;
                int wiggle = expt_container["start"]["random_triangulate"]["randomize"].as<int>(0);
                for (int a=0; a<expt_container["start"]["random_triangulate"]["count"].as<int>(1); ++a) {
                    do {
                        u = rand()%int(4*MAX_X/sep) - 2*MAX_X/sep;
                        v = rand()%int(4*MAX_Y/sep) - 2*MAX_Y/sep;;
                        tn = myNode ( COORD_TYPE(u*sep) + COORD_TYPE(v*sep)*0.5 + COORD_TYPE(rand()%(wiggle+1)),
                                        COORD_TYPE(v*sep)*0.6427876 + COORD_TYPE(rand()%(wiggle+1)) );
                    } while (my_map.isOutsideFrame (round(tn.x), round(tn.y)) || my_map.isObstacle (round(tn.x), round(tn.y)) || 
                                find (startNodes.begin(),startNodes.end(),tn)!=startNodes.end() );
                    startNodes.push_back (tn);
                    startNodes[a].put_in_grid();
                    startNodes[a].print("Start node: ");
                    colors.push_back (cvScalar (rand()%256, rand()%256, rand()%256));
                }
            }
        }
        
        // display options
        PLOT_SCALE = expt_container["plot_options"]["plot_scale"].as<double>(1.0);
        VERTEX_SIZE = expt_container["plot_options"]["vertex_size"].as<double>(1.0);
        LINE_THICKNESS = expt_container["plot_options"]["line_thickness"].as<double>(2.0); // CV_FILLED
        
        // saving options
        frameno = 0;
        imgPrefix << MAKESTR(_DOSL_ALGORITHM) << GRAPH_TYPE << "_map2d_Delaunay_";
        
        #if _VIS
        image_to_display = my_map.getCvMat (COLOR_MAP);
        cv::resize (image_to_display, image_to_display, cv::Size(), PLOT_SCALE , PLOT_SCALE );
        cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE);
        cv::imshow("Display window", image_to_display);
        cv::waitKey(0);
        #endif
        
        // Set planner variables
        all_nodes_set_p->reserve ( ceil(MAX_X - MIN_X + 1) );
    }
    
    // -----------------------------------------------------------
    
    bool isNodeInWorkspace (const myNode& tn) {
        if ( tn.x<MIN_X || tn.x>=MAX_X || tn.y<MIN_Y || tn.y>=MAX_Y )  return (false);
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
    
    // ======================================================================================
    // All the following functions are virtual members of 'Algorithm' class being overwritten
    
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
            
            s->push_back(tn);
            c->push_back(1.0);
        }
        
        #endif
    }
    
    // -----------------------------------------------------------
    // 'isSegmentFree' is required by ThetaStar
    
    bool isSegmentFree (myNode &n1, myNode &n2, double* c)
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
    
    double getHeuristics (myNode& n)
    {
        return (0.0); // Dijkstra's
    }
    
    // -----------------------------------------------------------
    
    std::vector<myNode> getStartNodes (void) 
    {
        #if _VIS
        for (int a=0; a<startNodes.size(); ++a)
            cv::circle (image_to_display, cv_plot_coord(startNodes[a].x,startNodes[a].y), VERTEX_SIZE*PLOT_SCALE, 
                                                                cvScalar (200.0, 150.0, 150.0), -1, 8);
        #endif
        return (startNodes);
    }
    
    // -----------------------------------------------------------
    
    std::vector<UnorderedIndexPair> delaunay_edges;
    
    void nodeEvent (myNode &n, unsigned int e) 
    {
        #if _VIS
        
        CvScalar col = cvScalar(0.0, 0.0, 0.0);;
        int thickness = -1; //lineThickness;
        double radFactor = 1.0;
        
        bool pauseForVis=false, drawVertex=true;
        
        // --------------------------------------------
        if (e & EXPANDED) {
            
            lastExpanded = n;
            bool cameFromNull = false;
            #ifdef DOSL_ALGORITHM_SStar
            cameFromNull = (n.CameFromSimplex==NULL);
            #endif
            if (!cameFromNull) {
                #if VERTEX_COLORS
                col = colors[n.lineage_data.id]; //cvScalar(150.0, 255.0, 150.0);
                #else
                col = cvScalar(255.0, 255.0, 255.0);
                #endif
            }
            else {
                printf ("expanded, but came-from is NULL!!\n");
                // pauseForVis = true;
            }
            
            //n.print("expanding: ");
            for (auto it=n.successors.begin(); it!=n.successors.end(); ++it)
                if (it->first->expanded  &&  it->first->lineage_data.id != n.lineage_data.id  && 
                            std::find (delaunay_edges.begin(), delaunay_edges.end(), 
                                            UnorderedIndexPair(n.lineage_data.id,it->first->lineage_data.id) ) == delaunay_edges.end() )
                    delaunay_edges.push_back (UnorderedIndexPair(n.lineage_data.id,it->first->lineage_data.id));
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
            // printf ("Backtrack started!!\n");
        }
        
        else
            return;
                
        //-------------------------------------------
        if (drawVertex) {
            double nodeRad = radFactor*VERTEX_SIZE;
            if (n.x<MAX_X  &&  n.y<MAX_Y  &&  my_map.isFree(round(n.x), round(n.y)) )
                cvPlotPoint (cv_plot_coord(n.x,n.y), col, PLOT_SCALE);
        }
        
        if (expand_count%VIS_INTERVAL == 0  ||  node_heap_p->size() == 0) {
            cv::imshow("Display window", image_to_display);
            std::cout << std::flush;
            if (SAVE_IMG_INTERVAL>0 && expand_count%SAVE_IMG_INTERVAL == 0) {
                char imgFname[1024];
                sprintf(imgFname, "%s%s%05d.png", out_folderName.c_str(), imgPrefix.str().c_str(), expand_count);
                cv::imwrite(imgFname, image_to_display);
            }
            cvWaitKey(1); //(10);
        }
        
        if (pauseForVis) {
            cvWaitKey();
        }
        #endif
    }
    
    // ---------------------------------------
    
};

// ==============================================================================


int main(int argc, char *argv[])
{
    // Read command-line parameters:
    compute_program_path();
    
    std::string expt_f_name = program_path+"../files/expt/basic_experiments.json", expt_name="L457_partition_rand";
    if (argc == 2) {
        expt_name = argv[1];
    }
    else if (argc > 2) {
        expt_f_name = argv[1];
        expt_name = argv[2];
    }
    
    printf (_BOLD _YELLOW "Note: " YELLOW_ BOLD_ "Using algorithm " _YELLOW  MAKESTR(_DOSL_ALGORITHM)  YELLOW_ 
                ". Run 'make' to recompile with a different algorithm.\n");
    
    // -------------------------------------
    
    // search:
    searchProblem test_search_problem (expt_f_name, expt_name, program_path+"../files/out/");
    test_search_problem.search();
    
    // ---------------------------------------
    
    // print statistics
    #if _STAT
    char statFname[1024];
    sprintf(statFname, "%s%s_%s.txt", test_search_problem.out_folderName.c_str(), 
                    test_search_problem.imgPrefix.str().c_str(), test_search_problem.map_image_fName.c_str());
    FILE* pFile;
    pFile = fopen (statFname,"a+");
    if (pFile!=NULL) {
        fprintf (pFile, "%s, %f. %f;\n", test_search_problem.expt_Name.c_str(), 0.0, cost);
        fclose (pFile);
    }
    #endif
    
    double tc;
    
    char edgesFname[1024]; sprintf(edgesFname, "%s_DELAUNAY_edges_(id1,id2).csv", test_search_problem.out_folderName.c_str());
    FILE* eFile; eFile = fopen (edgesFname,"a+");
    for (auto it=test_search_problem.delaunay_edges.begin(); it!=test_search_problem.delaunay_edges.end(); ++it)
        if ( test_search_problem.isSegmentFree (test_search_problem.startNodes[it->i], test_search_problem.startNodes[it->j], &tc) ) {
            fprintf(eFile, "%d, %d\n", it->i, it->j);
            cv::line (test_search_problem.image_to_display, 
                        test_search_problem.cv_plot_coord (test_search_problem.startNodes[it->i].x, test_search_problem.startNodes[it->i].y), 
                            test_search_problem.cv_plot_coord (test_search_problem.startNodes[it->j].x, test_search_problem.startNodes[it->j].y), 
                                cvScalar(200.0,200.0,200.0), 2.0); //test_search_problem.LINE_THICKNESS*test_search_problem.PLOT_SCALE
        }
    fclose (eFile);
    
    char nodesFname[1024]; sprintf(nodesFname, "%s_DELAUNAY_vetices_(id,x,y).csv", test_search_problem.out_folderName.c_str());
    FILE* nFile; nFile = fopen (nodesFname,"a+");
    for (int a=0; a<test_search_problem.startNodes.size(); ++a) {
        fprintf(nFile, "%d, %f, %f\n", a, test_search_problem.startNodes[a].x, test_search_problem.startNodes[a].y);
        cv::circle (test_search_problem.image_to_display, 
                    test_search_problem.cv_plot_coord (test_search_problem.startNodes[a].x, test_search_problem.startNodes[a].y),
                        3.0, cvScalar(220.0,100.0,100.0), -1, 8);
    }
    fclose (nFile);
    
    cv::imshow("Display window", test_search_problem.image_to_display);
    
    // save image file
    #if _VIS
    if (SAVE_IMG_INTERVAL != 0) {
        char imgFname[1024];
        sprintf(imgFname, "%s%s_DELAUNAY.png", test_search_problem.out_folderName.c_str(), 
                                                test_search_problem.imgPrefix.str().c_str());
        cv::imwrite(imgFname, test_search_problem.image_to_display);
    }
    cvWaitKey();
    #endif
    
    
    // clear memory
    test_search_problem.clear();
}

