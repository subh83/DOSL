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

#define GRAPH_TYPE 8 // 6 or 8. 8 for regular bitmap image


// ==============================================================================

// A node of the graph
template <class doslAlgorithm>
class cvPPnode : public doslAlgorithm::template Node< cvPPnode<doslAlgorithm>, double>
{
public:
    int x, y;
    // std::vector<int> h;
    
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
    
    bool isCoordsEqual(const cvPPnode& n) const {
        #if false && GRAPH_TYPE == 8 // int
        return (isEqual_i(x,n.x) && isEqual_i(y,n.y));
        #else
        return (isEqual_d(x,n.x) && isEqual_d(y,n.y));
        #endif
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
        #if GRAPH_TYPE == 8 // int
        return (abs(x));
        #else
        return ( MAX(round(fabs(x)+INFINITESIMAL_DOUBLE), round(fabs(x)-INFINITESIMAL_DOUBLE)) );
        #endif
    }
    
    // print
    void print (std::string head="", std::string tail="\n") const {
        _dosl_cout << _GREEN + head << " (" << this << ")" GREEN_ " x=" << x << ", y=" << y << "; ";
        (this->G == std::numeric_limits<double>::max())? printf("INF") : printf("G = %0.8f \n", this->G);
        _dosl_cout << "Successors: ";
        for (auto it=this->Successors.begin(); it!=this->Successors.end(); ++it)
            printf ("%x (%f), ", it->first, it->second);
        std::cout << tail << _dosl_endl;
    }
};

// ==============================================================================


template <class doslAlgorithm>
class cvPathPlanner : public doslAlgorithm::template Algorithm< cvPPnode<doslAlgorithm>, double>
{
public:
    typedef cvPPnode<doslAlgorithm> cvPPnode_;
    
    cvParseMap2d  my_map;
    
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
    int MAX_X, MIN_X, MAX_Y, MIN_Y;
    cvPPnode_ startNode, goalNode, lastExpanded;
    
    // =========================================
    
    template<typename T>
    cv::Point cv_plot_coord(T x, T y) { return( cv::Point( round(PLOT_SCALE*(x-MIN_X)), round(PLOT_SCALE*(y-MIN_Y)) ) ); }
    
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
    
    
    void find_path (void) {
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
    cvPathPlanner (cv::Mat obs_map, cv::Point start, cv::Point goal,
                        bool vis=false)
    {
        // int obsSizeThresh = 0; // MIN(obs_map.cols,obs_map.rows) / 4;
        my_map = cvParseMap2d (obs_map, false);
        
        // read data for planning
        MAX_X=my_map.width(); MIN_X=0; MAX_Y=my_map.height(); MIN_Y=0;
        startNode = cvPPnode_ (start.x, start.y);
        goalNode = cvPPnode_ (goal.x, goal.y);
        startNode.put_in_grid(); goalNode.put_in_grid();
        
        // display options
        PLOT_SCALE = 1.0;
        VERTEX_SIZE = 1.0;
        LINE_THICKNESS = 2.0;
        VERTEX_COLORS = 1; // 1 or 0
        VIS_INTERVAL = 100;
        
        // saving options
        //frameno = 0;
        //imgPrefix << MAKESTR(_DOSL_ALGORITHM) << GRAPH_TYPE << "homotopy2d_";
        
        visualize = vis;
        
        // Set planner variables
        this->AllNodesSet.HashTableSize = ceil(MAX_X - MIN_X + 1);
        
        // compute path
        // find_path();
    }
    
    // -----------------------------------------------------------
    
    bool isNodeInWorkspace (const cvPPnode_& tn) {
        if ( tn.x<MIN_X || tn.x>=MAX_X || tn.y<MIN_Y || tn.y>=MAX_Y )  return (false);
        return (true);
    }
    
    bool isEdgeAccessible (const cvPPnode_& tn1, const cvPPnode_& tn2) {
        if ( (!isNodeInWorkspace(tn1)) || !(isNodeInWorkspace(tn2)) )  return (false);
        
        #if GRAPH_TYPE == 8 // int int
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
    
    void getSuccessors (cvPPnode_ &n, std::vector<cvPPnode_>* s, std::vector<double>* c) // *** This must be defined
    {
        // This function should account for obstacles and size of environment.
        cvPPnode_ tn;
        
        #if GRAPH_TYPE == 8
        for (int a=-1; a<=1; ++a)
            for (int b=-1; b<=1; ++b) {
                if (a==0 && b==0) continue;
                
                if (doslAlgorithm::AlgorithmName == "SStar") { // use triangulation for 'SStar'
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
        
        #elif GRAPH_TYPE == 6
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
        
        #endif
        
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
        /* double dx = goalNode.x - n.x;
        double dy = goalNode.y - n.y;
        return (sqrt(dx*dx + dy*dy)); */
        return (0.0);
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
            if (n.x<MAX_X  &&  n.y<MAX_Y  &&  my_map.isFree(round(n.x), round(n.y)) )
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
