#ifndef __CV_MULTI_CLASS_PATH_PLANNER_TCC
#define __CV_MULTI_CLASS_PATH_PLANNER_TCC

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


double wiggle = 1e-3;

// ==============================================================================

// A node of the graph
template <class doslAlgorithm>
class cvMPPnode : public doslAlgorithm::template Node< cvMPPnode<doslAlgorithm>, double>
{
public:
    int x, y;
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
    
    bool isCoordsEqual(const cvMPPnode& n) const {
        #if false && GRAPH_TYPE == 8 // int
        return (isEqual_i(x,n.x) && isEqual_i(y,n.y));
        #else
        return (isEqual_d(x,n.x) && isEqual_d(y,n.y));
        #endif
    }
    
    // *** This must be defined for the node
    bool operator==(const cvMPPnode& n) const {
        if (!isCoordsEqual(n)) return (false);
        if (h.size()!=n.h.size()) return (false);
	    for (int a=0; a<h.size(); a++)
	        if (h[a]!=n.h[a]) return (false);
        return (true);
    }
    
    // constructor
    cvMPPnode () { }
    cvMPPnode (int xx, int yy) : x(xx), y(yy) { put_in_grid(); }
    
    
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
        (this->G == std::numeric_limits<double>::max())? printf("INF") : printf("G = %0.8f; h = [", this->G);
        for (int a=0; a<h.size(); ++a){
            printf("%d", a);
            if (a!=h.size()-1) printf(", ");
        } 
        printf("]\n");
        _dosl_cout << "Successors: ";
        for (auto it=this->Successors.begin(); it!=this->Successors.end(); ++it)
            printf ("%x (%f), ", it->first, it->second);
        std::cout << tail << _dosl_endl;
    }
};

// ==============================================================================

enum H_CLASS_TYPE { HOMOLOGY_CLASS, HOMOTOPY_CLASS };

template <class doslAlgorithm, int HClassType=HOMOTOPY_CLASS>
class cvMulticlassPathPlanner : public doslAlgorithm::template Algorithm< cvMPPnode<doslAlgorithm>, double>
{
public:
    typedef cvMPPnode<doslAlgorithm> cvMPPnode_;
    
    cvParseMap2d  my_map;
    int nClassesToFind;
    
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
    cvMPPnode_ startNode, goalNode, lastExpanded;
    
    // homotopy classes
    int nClasses;
    std::vector<cvMPPnode_> homotopyGoals;
    
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
    
    std::vector< std::vector< cv::Point > > paths;
    std::vector<double> costs;
    
    
    void find_paths (void) {
        this->search();
        
        for (int i=0; i<homotopyGoals.size(); ++i) {
            // get path
            auto path = this->reconstructPath (homotopyGoals[i]);
            double cost = 0.0;
            CvScalar colr = cvScalar (50.0 + rand()%150, 50.0 + rand()%150, 50.0 + rand()%150);
            
            cv::Point2d lastPt, thisPt(startNode.x,startNode.y);
            std::vector<cv::Point> allCvPts;
            for (int a=path.size()-1; a>=0; --a) {
                lastPt = thisPt;
                thisPt = cv::Point2d(0.0, 0.0);
                for (auto it=path[a].begin(); it!=path[a].end(); ++it) {
                    thisPt.x += it->second * it->first->x;
                    thisPt.y += it->second * it->first->y;
                }
                allCvPts.push_back (cv::Point((int)round(thisPt.x),(int)round(thisPt.y)));
                if (visualize)
                    cv::line (image_to_display, cv_plot_coord(thisPt.x,thisPt.y), cv_plot_coord(lastPt.x,lastPt.y), 
                                            colr, LINE_THICKNESS*PLOT_SCALE);
                cost += sqrt ((thisPt.x-lastPt.x)*(thisPt.x-lastPt.x) + (thisPt.y-lastPt.y)*(thisPt.y-lastPt.y));
            }
            paths.push_back (allCvPts);
            costs.push_back (cost);
        }
        
        if (visualize) {
            cv::imshow("Display window", image_to_display);
            cv::waitKey(0);
        }
    }
    
    // -----------------------------------------------------------
    
    // Constructor
    cvMulticlassPathPlanner (cv::Mat obs_map, cv::Point start, cv::Point goal, int nPaths=1, 
                                bool vis=false, int obsSizeThresh=0)
    {
        // int obsSizeThresh = 0; // MIN(obs_map.cols,obs_map.rows) / 4;
        my_map = cvParseMap2d (obs_map, true, obsSizeThresh);
        
        // read data for planning
        MAX_X=my_map.width(); MIN_X=0; MAX_Y=my_map.height(); MIN_Y=0;
        startNode = cvMPPnode_ (start.x, start.y);
        goalNode = cvMPPnode_ (goal.x, goal.y);
        startNode.put_in_grid(); goalNode.put_in_grid(); 
        nClassesToFind = nPaths;
        nClasses = 0;
        
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
        if (visualize) {
            image_to_display = my_map.getCvMat (COLOR_MAP);
            cv::resize (image_to_display, image_to_display, cv::Size(), PLOT_SCALE , PLOT_SCALE );
            cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE);
            cv::imshow("Display window", image_to_display);
            //cv::waitKey(0);
        }
        
        // Set planner variables
        this->AllNodesSet.HashTableSize = ceil(MAX_X - MIN_X + 1);
        
        // compute path
        find_paths();
    }
    
    // -----------------------------------------------------------
    
    bool isNodeInWorkspace (const cvMPPnode_& tn) {
        if ( tn.x<MIN_X || tn.x>=MAX_X || tn.y<MIN_Y || tn.y>=MAX_Y )  return (false);
        return (true);
    }
    
    bool isEdgeAccessible (const cvMPPnode_& tn1, const cvMPPnode_& tn2) {
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
    
    void updateHSignature (cvMPPnode_ &n, cvMPPnode_ &tn) { // updates h/H-signature of tn
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
                    if (HClassType == HOMOTOPY_CLASS) { // homotopy
                        if (tn.h.size()>0 && tn.h[tn.h.size()-1]==-pm*(p+1))
                            tn.h.pop_back();
                        else
                            tn.h.push_back(pm*(p+1));
                    } else // homology
                        tn.h[p] += pm;
	            }
            }
        }
    }
    
    void getSuccessors (cvMPPnode_ &n, std::vector<cvMPPnode_>* s, std::vector<double>* c) // *** This must be defined
    {
        // This function should account for obstacles and size of environment.
        cvMPPnode_ tn;
        
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
    
    double getHeuristics (cvMPPnode_& n)
    {
        /* double dx = goalNode.x - n.x;
        double dy = goalNode.y - n.y;
        return (sqrt(dx*dx + dy*dy)); */
        return (0.0);
    }
    
    // -----------------------------------------------------------
    
    std::vector<cvMPPnode_> getStartNodes (void) 
    {
        std::vector<cvMPPnode_> startNodes;
        
        startNodes.push_back (startNode);
        if (visualize)
            cv::circle (image_to_display, cv_plot_coord(startNode.x,startNode.y), VERTEX_SIZE*PLOT_SCALE, 
                                                                cvScalar (200.0, 150.0, 150.0), -1, 8);
        
        return (startNodes);
    }
    
    // -----------------------------------------------------------
    
    void nodeEvent (cvMPPnode_ &n, unsigned int e) 
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
    
    bool stopSearch (cvMPPnode_ &n) {
        if (n.isCoordsEqual(goalNode)) {
            homotopyGoals.push_back (n);
            ++nClasses;
            // n.print ("Found a path to ");
            if (nClasses>=nClassesToFind)
                return (true);
        }
        return (false);
    }
    
    // =========================================
    
    void draw_paths (cv::Mat& in_map, cv::Mat& out_map, std::vector<CvScalar> colors=std::vector<CvScalar>(),
                                            int thickness=1, int lineType=8, int shift=0) {
        int nChannels = out_map.channels();
        
        if (&in_map != &out_map) {
            in_map.copyTo (out_map);
            if (out_map.channels()<3 && nChannels==3)
                cv::cvtColor (out_map, out_map, CV_GRAY2RGB);
        }
        
        if (colors.size()<paths.size()) {
            if (colors.size()==1) {
                colors.resize (paths.size());
                for (int i=0; i<paths.size(); ++i)
                    colors[i] = colors[0];
            }
            else {
                colors.resize (paths.size());
                for (int i=0; i<paths.size(); ++i)
                    colors[i] = cvScalar (50.0 + rand()%200, 50.0 + rand()%200, 50.0 + rand()%200);
            }
        }
        
        for (int i=0; i<paths.size(); ++i)
            for (int j=0; j<paths[i].size()-1; ++j)
                cv::line (out_map, paths[i][j], paths[i][j+1], colors[i], thickness, lineType, shift);
    }
    
    
    cv::Mat draw_paths (std::vector<CvScalar> colors=std::vector<CvScalar>(),
                                    int thickness=1, int lineType=8, int shift=0) {
        cv::Mat ret = my_map.getCvMat (COLOR_MAP, FREE_MAP, 0);
        draw_paths (ret, ret, colors, thickness, lineType, shift);
        return (ret);
    }
};

// ==============================================================================


#endif
