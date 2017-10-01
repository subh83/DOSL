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

// Other libraries:
// Open CV:
#include <opencv2/opencv.hpp> 
#include <opencv2/highgui.hpp> 

// DOSL library

#ifndef _DOSL_ALGORITHM // can pass at command line during compilation: -D_DOSL_ALGORITHM=AStar
    #define _DOSL_ALGORITHM  AStar
#endif
#include <dosl/dosl>

// Local libraries/headers
#include <dosl/aux-utils/cvParseMap2d.hpp>
#include <dosl/aux-utils/double_utils.hpp>
#include <dosl/aux-utils/string_utils.hpp> // compute_program_path
#include "../../include-local/RSJparser.tcc"
#include "../../include-local/mathevalAux.h"
#include "../../include-local/map2d.h"

// =======================
// Graph parameters
#define GRAPH_TYPE 8 // 6 or 8
#define COORD_TYPE double // only 'double' is supported

// output options
#define _STAT 0
#define _VIS 1
#define VIS_INTERVAL 100
#define VERTEX_COLORS 1

// ==============================================================================

// A node of the graph
class myNode : public _DOSL_ALGORITHM::Node<myNode,double>
{
public:
    COORD_TYPE x, y;
    
    #if GRAPH_TYPE == 8
        void put_in_grid (void) { x=(COORD_TYPE)round(x); y=(COORD_TYPE)round(y); }
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
    
    // *** This must be defined for the node
    bool operator==(const myNode& n) const {
        #if GRAPH_TYPE == 8 // COORD_TYPE int
        return ((x==n.x) && (y==n.y));
        #else
        return (fabs(x-n.x)<INFINITESIMAL_DOUBLE  &&  fabs(y-n.y)<INFINITESIMAL_DOUBLE);
        #endif
    }
    
    // constructor
    myNode () { }
    myNode (COORD_TYPE xx, COORD_TYPE yy) : x(xx), y(yy) { put_in_grid(); }
    template<class PtType> myNode (PtType p) : x(p.x), y(p.y) { put_in_grid(); }
    
    // Inherited functions being overwritten
    int getHashBin (void) const {
        #if GRAPH_TYPE == 8 // COORD_TYPE int
        return (abs(x));
        #else
        return ( MAX(round(fabs(x)+INFINITESIMAL_DOUBLE), round(fabs(x)-INFINITESIMAL_DOUBLE)) );
        #endif
    }
    
    // print
    void print (std::string head="", std::string tail="") const {
        _dosl_cout << _GREEN + head << " (" << this << ")" GREEN_ " x=" << x << ", y=" << y << ", ";
        (G==std::numeric_limits<double>::max())? printf("INF") : printf("%0.8f", G);
        printf ("\n");
        _dosl_cout << "\tExpanded=%d" << Expanded << _dosl_endl;
        _dosl_cout << "Successors: ";
        for (auto it=Successors.begin(); it!=Successors.end(); ++it)
            printf ("%x (%f), ", it->first, it->second);
        std::cout << tail << _dosl_endl;
    }
    
};

// ==============================================================================

class searchProblem : public _DOSL_ALGORITHM::Algorithm<myNode,double>
{
public:
    // Fime names and JSON objects
    std::string   map_image_fName, expt_fName, expt_folderName, expt_Name, out_folderName;
    RSJresource expt_container;
    cvParseMap2d  my_map;
    
    // Image display variables / parameters
    cv::Mat image_to_display;
    double plot_scale, vertex_size, line_thickness;
    int save_image_interval;
    bool do_animation, pause_for_visualization;
    
    // variables for saving image
    int frameno;
    std::ostringstream imgPrefix;
    
    // variables decsribing problem scales
    COORD_TYPE width, height;
    double x_min, x_max, y_min, y_max;
    
    // other parameters
    double intersection_check_steps_per_unit_length;
    
    // start & goal
    myNode startNode, goalNode, lastExpanded;
    
    // Computed:
    //----------
    std::unordered_map<std::string,double> sym_vals; // symbol values
    
    // metric map
    //std::vector< std::vector< std::vector<double> > > metricMap; // metricMap[x][y] == {g11, g12, g21, g22}
    std::vector<map2d> metric; // g11, g12, g21, g2
    map2d det;
    double detMax, detMin; // for plotting
    
    // tools
    PixelChartTransformation  pc_transform;
    
    // -----------------------------------------------------------
    
    // TODO: replace these with openCV's scaling function
    
    template<typename T>
    CvPoint cv_plot_coord(T x, T y) { return( cvPoint( round(plot_scale*(x)), round(plot_scale*(y)) ) ); }
    
    void cvPlotPoint (CvPoint pt, CvScalar colr, int size=1) { 
        for (int i=pt.x; i<=pt.x+(size-1); i++)
            for (int j=pt.y; j<=pt.y+(size-1); j++) {
                if (i<0 || i>=image_to_display.cols || j<0 || j>=image_to_display.rows) continue;
                image_to_display.at<cv::Vec3b>(j,i) = cv::Vec3b ((uchar)colr.val[0], (uchar)colr.val[1], (uchar)colr.val[2]);
            }
    }
    
    // -----------------------------------------------------------
    
    // Constructor
    searchProblem (std::string expt_f_name, std::string expt_name, std::string out_folder_name="outfiles/")
    {
        expt_fName = expt_f_name; expt_Name = expt_name;
        expt_folderName = expt_fName.substr(0, expt_fName.find_last_of("/\\")+1);
        out_folderName = out_folder_name;
        
        // Read environment details:
        // ------------------------
        std::ifstream my_fstream (expt_fName);
        expt_container = RSJresource (my_fstream)[expt_Name];
        
        // onstacle map
        if (expt_container["environment"]["pixmap"].exists()) {
            map_image_fName = expt_folderName + expt_container["environment"]["pixmap"].as<std::string>();
            my_map = cvParseMap2d (map_image_fName, false);
        }
        else if (expt_container["environment"]["width"].exists() && expt_container["environment"]["height"].exists()) {
            my_map = cvParseMap2d (cv::Mat (MathEvaluator(expt_container["environment"]["height"].as<std::string>()).eval(), 
                                            MathEvaluator(expt_container["environment"]["width"].as<std::string>()).eval(),
                                                            CV_8UC3, cv::Scalar(255,255,255)), false);
        }
        
        // read and compute scale data
        // ---------------------------
        width=my_map.width(); height=my_map.height();
        sym_vals = std::unordered_map<std::string,double> ({{"width",width}, {"height",height}});
        
        x_min = MathEvaluator(expt_container["coord_scale"]["x_min"].as<std::string>("-0.5"), sym_vals).eval();
        x_max = MathEvaluator(expt_container["coord_scale"]["x_max"].as<std::string>("width-0.5"), sym_vals).eval();
        y_min = MathEvaluator(expt_container["coord_scale"]["y_min"].as<std::string>("-0.5"), sym_vals).eval();
        y_max = MathEvaluator(expt_container["coord_scale"]["y_max"].as<std::string>("height-0.5"), sym_vals).eval();
        
        sym_vals = std::unordered_map<std::string,double> ( {
                                                {"width",width}, {"height",height}, 
                                                {"x_min",x_min}, {"x_max",x_max}, {"y_min",y_min}, {"y_max",y_max},
                                                {"x",0.0}, {"y",0.0} } );
        
        int chart_limits = (expt_container["coord_scale"]["chart_limits"].as<std::string>("pixel_boundary") 
                                                                        == "pixel_boundary")? PIXEL_BOUNDARY : PIXEL_CENTER;
        pc_transform = PixelChartTransformation (width, height, x_min, x_max, y_min, y_max, chart_limits);
        pc_transform.print ("Map dimensions:");
        
        
        
        // start and goal
        // --------------
        startNode = myNode ( ChartPoint (
                                        MathEvaluator (expt_container["start"][0].as<std::string>(), sym_vals).eval(),
                                        MathEvaluator (expt_container["start"][1].as<std::string>(), sym_vals).eval()  
                                        ).toMapPoint (pc_transform) );
        goalNode = myNode ( ChartPoint (
                                        MathEvaluator (expt_container["goal"][0].as<std::string>(), sym_vals).eval(),
                                        MathEvaluator (expt_container["goal"][1].as<std::string>(), sym_vals).eval()  
                                        ).toMapPoint (pc_transform) );
        startNode.put_in_grid();
        goalNode.put_in_grid(); 
        
        
        // cost/metric map
        // ---------------
        bool fun_to_pixmap = expt_container["cost"]["fun_to_pixmap"].as<bool>(false);
        metric.resize(4);
        
        if (expt_container["cost"]["pixmap"].exists()) {
            std::string cost_map_fName = expt_folderName + expt_container["cost"]["pixmap"].as<std::string>();
            metric[0] = map2d (cost_map_fName, pc_transform);
            metric[1] = map2d (0.0);
            metric[2] = map2d (0.0);
            metric[3] = map2d (cost_map_fName, pc_transform);
        }
        else if (expt_container["cost"]["fun"].exists()) {
            metric[0] = map2d (MathEvaluator (expt_container["cost"]["fun"].as<std::string>(), sym_vals), 
                                                                                    pc_transform, fun_to_pixmap);
            metric[1] = map2d (0.0);
            metric[2] = map2d (0.0);
            metric[3] = map2d (MathEvaluator (expt_container["cost"]["fun"].as<std::string>(), sym_vals), 
                                                                                    pc_transform, fun_to_pixmap);
        }
        else if (expt_container["cost"]["metric"].exists()) {
            metric[0] = map2d (MathEvaluator (expt_container["cost"]["metric"][0].as<std::string>(), sym_vals), 
                                                                                    pc_transform, fun_to_pixmap);
            metric[1] = map2d (MathEvaluator (expt_container["cost"]["metric"][1].as<std::string>(), sym_vals), 
                                                                                    pc_transform, fun_to_pixmap);
            metric[2] = map2d (MathEvaluator (expt_container["cost"]["metric"][2].as<std::string>(), sym_vals), 
                                                                                    pc_transform, fun_to_pixmap);
            metric[3] = map2d (MathEvaluator (expt_container["cost"]["metric"][3].as<std::string>(), sym_vals), 
                                                                                    pc_transform, fun_to_pixmap);
        }
        else {
            metric[0] = map2d (1.0);
            metric[1] = map2d (0.0);
            metric[2] = map2d (0.0);
            metric[3] = map2d (1.0);
        }
        
        
        detMax=std::numeric_limits<double>::min(), detMin=std::numeric_limits<double>::max();
        det = map2d (width, height);
        for (int x=0; x<width; ++x)
            for (int y=0; y<height; ++y) {
                det.at(x,y) = metric[0].evalPixel(Pixel(x,y))*metric[3].evalPixel(Pixel(x,y)) 
                                    - metric[1].evalPixel(Pixel(x,y))*metric[2].evalPixel(Pixel(x,y));
                if (det.at(x,y)>detMax) detMax = det.at(x,y);
                if (det.at(x,y)<detMin) detMin = det.at(x,y);
            }
        
        // other parameters
        // ----------------
        
        intersection_check_steps_per_unit_length = expt_container["search_options"]
                                                        ["intersection_check_steps_per_unit_length"].as<double>(10.0);
        
        
        // display options
        // -----------------
        
        plot_scale = MathEvaluator(expt_container["plot_options"]["plot_scale"].as<std::string>("1.0"), sym_vals).eval();
        vertex_size = expt_container["plot_options"]["vertex_size"].as<double>(1.0);
        line_thickness = expt_container["plot_options"]["line_thickness"].as<double>(2.0);
        do_animation = expt_container["plot_options"]["do_animation"].as<bool>(true);
        save_image_interval = expt_container["plot_options"]["save_image_interval"].as<int>(-1.0);
        pause_for_visualization = expt_container["plot_options"]["pause_for_visualization"].as<bool>(true);
        
        // Initial plotting
        #if _VIS
        image_to_display = my_map.getCvMat (COLOR_MAP);
        cv::resize (image_to_display, image_to_display, cv::Size(), plot_scale , plot_scale );
        cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE);
        #endif
        
        #if _VIS
        double scaledDet, baseIntensity=0.5;
        for (int x=0; x<width; ++x)
            for (int y=0; y<height; ++y) 
                if (my_map.isFree(x,y)) {
                    scaledDet = (det.at(x,y) - detMin) / (detMax - detMin);
                    cvPlotPoint (cv_plot_coord(x,y), 
                                    cvScalar(255.0*baseIntensity, 255.0*(1.0-(1-baseIntensity)*scaledDet), 
                                                255.0*(baseIntensity+(1-baseIntensity)*scaledDet)),
                                    plot_scale);
            }
        cv::imshow("Display window", image_to_display);
        if (pause_for_visualization)
            cv::waitKey(0);
        #endif
        
        // saving options
        frameno = 0;
        imgPrefix << MAKESTR(_DOSL_ALGORITHM) << GRAPH_TYPE << "_" << expt_Name << "_";
        
        // Set planner variables
        AllNodesSet.HashTableSize = ceil(width + 1);
    }
    
    // -----------------------------------------------------------
    
    double getSmallSegmentCost (MapPoint p1, MapPoint p2) { // pixel coordinates
        ChartPoint cp1 = p1.toChartPoint (pc_transform);
        ChartPoint cp2 = p2.toChartPoint (pc_transform);
        double dx=cp2.x-cp1.x, dy=cp2.y-cp1.y;
        double _g11 = (metric[0].evalMapPoint (p1) + metric[0].evalMapPoint (p2)) / 2.0,
               _g12 = (metric[1].evalMapPoint (p1) + metric[1].evalMapPoint (p2)) / 2.0,
               _g21 = (metric[2].evalMapPoint (p1) + metric[2].evalMapPoint (p2)) / 2.0,
               _g22 = (metric[3].evalMapPoint (p1) + metric[3].evalMapPoint (p2)) / 2.0;
       return (sqrt (dx*dx*_g11 + dx*dy*(_g12+_g21) + dy*dy*_g22) );
    }
   
    
    bool isNodeInWorkspace (MapPoint n) {
        return ( n.inDomain (pc_transform) );
    }
    
    double isNodeAccessible (MapPoint n) {
        Pixel pix = n.toPixel (pc_transform);
        return (
                isNodeInWorkspace (n)
                && my_map.isFree (pix.x, pix.y)
                );
    }
    
    // -----------------------------------------------------------
    
    // only for Theta*
    bool isSegmentFree (myNode &n1, myNode &n2, double* c) {
        *c = 0.0;
        double dx=(double)(n2.x-n1.x), dy=(double)(n2.y-n1.y);
        double step_count = intersection_check_steps_per_unit_length * ceil(sqrt(dx*dx+dy*dy));
        
        if (step_count==0.0)
            return (isNodeAccessible(n1));
        else if ( !isNodeAccessible(n1) || !isNodeAccessible(n2) )
            return (false);
        
        double xstep=dx/step_count, ystep=dy/step_count;
        //double tn1x=n1.x, tn1y=n1.y, tn2x, tn2y;
        MapPoint tn1(n1), tn2;
        for (int a=0; a<step_count; ++a) {
            tn2=tn1;
            tn2.x += xstep; tn2.y += ystep;
            if ( !isNodeAccessible(tn2) ) return (false);
            (*c) += getSmallSegmentCost (tn1, tn2);
            tn1=tn2;
        }
        return (true);
    }
    
    // -----------------------------------------------------------
    
    void getSuccessors (myNode &n, std::vector<myNode>* s, std::vector<double>* c) // *** This must be defined
    {
        // This function should account for obstacles and size of environment.
        myNode tn;
        double tc;
        
        #if GRAPH_TYPE == 8
        for (int a=-1; a<=1; ++a)
            for (int b=-1; b<=1; ++b) {
                if (a==0 && b==0) continue;
                
                #if defined(DOSL_ALGORITHM_SStar)
                int xParity = ((int)round(fabs(n.x))) % 2;
                if (xParity==0 && (a!=0 && b==-1)) continue;
                if (xParity==1 && (a!=0 && b==1)) continue;
                #endif
                
                tn.x = n.x + a;
                tn.y = n.y + b;
                
                //if (!isEdgeAccessible(tn,n)) continue;
                if  (!isSegmentFree (tn, n, &tc)) continue;
                
                s->push_back(tn);
                //c->push_back( getSmallSegmentCost(n,tn) ); 
                c->push_back (tc);
            }
        
        #elif GRAPH_TYPE == 6
        double th, dx, dy;
        for (int a=0; a<6; ++a) {
            th = a * PI_BY_3;
            tn.x = n.x + 1.0*cos(th);
            tn.y = n.y + 1.0*sin(th);
            
            //if (!isEdgeAccessible(tn,n)) continue;
            if  (!isSegmentFree (tn, n, &tc)) continue;
            
            s->push_back(tn);
            //c->push_back( getSmallSegmentCost(n,tn) );
            c->push_back (tc);
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
        cv::circle (image_to_display, cv_plot_coord(startNode.x,startNode.y), vertex_size*plot_scale, 
                                                                cvScalar (200.0, 150.0, 150.0), -1, 8);
        #endif
        
        return (startNodes);
    }
    
    // -----------------------------------------------------------
    
    void nodeEvent (myNode &n, unsigned int e) 
    {
        #if _VIS
        if (!do_animation) return;
        
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
                col = cvScalar(150.0, 255.0, 150.0);
                #else
                col = cvScalar(255.0, 255.0, 255.0);
                #endif
            }
            else {
                // printf ("Expanded, but came-from is NULL!!\n");
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
        }
        
        else
            return;
                
        //-------------------------------------------
        if (drawVertex) {
            double nodeRad = radFactor*vertex_size;
            if (n.x<width  &&  n.y<height  &&  my_map.isFree(round(n.x), round(n.y)) )
                cvPlotPoint (cv_plot_coord(n.x,n.y), col, plot_scale);
        }
        
        if (ExpandCount%VIS_INTERVAL == 0  ||  NodeHeap.size() == 0) {
            cv::imshow("Display window", image_to_display);
            std::cout << std::flush;
            if (save_image_interval>0 && ExpandCount%save_image_interval == 0) {
                char imgFname[1024];
                sprintf(imgFname, "%s%s%05d.png", out_folderName.c_str(), imgPrefix.str().c_str(), ExpandCount);
                cv::imwrite(imgFname, image_to_display);
            }
            cvWaitKey(1); //(10);
        }
        
        if (pause_for_visualization && pauseForVis) {
            cvWaitKey();
        }
        #endif
    }
    
    // ---------------------------------------
    
    bool stopSearch (myNode &n) {
        return (n==goalNode);
    }
};

// ==============================================================================

int main(int argc, char *argv[])
{
    //RUNTIME_VERBOSE_SWITCH = 0;
    compute_program_path();
    
    std::string expt_f_name=program_path+"../files/expt/nonuniform_experiments.json", expt_name="empty_nonuniform_small";
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
    
    // run search
    searchProblem test_search_problem (expt_f_name, expt_name, program_path+"../files/out/");
    test_search_problem.search();
    
    // get path (node, weight tuples)
    std::vector < std::unordered_map<myNode*,double> > path = test_search_problem.reconstructPath (test_search_problem.goalNode);
    
    // -------------------------------------
    
    // Print and draw path
    std::ostream* path_ostream_p = &(std::cout);
    std::ofstream path_ofstream;
    if (test_search_problem.expt_container["print_options"]["path_output_file"].exists()) {
        path_ofstream.open (test_search_problem.expt_container["print_options"]["path_output_file"].as<std::string>().c_str(),
                                            std::fstream::app);
        path_ostream_p = &(path_ofstream);
    }
    
    #if _VIS
    *(path_ostream_p) << "\nPaths.('"<<test_search_problem.imgPrefix.str()<<"') = [";
    #endif
    myNode thisPt=test_search_problem.startNode, lastPt;
    std::vector<myNode> allPts;
    double cost = 0.0, tc;
    
    for (int a=path.size()-1; a>=0; --a) {
        //compute convex combination (relevant to SStar)
        lastPt = thisPt;
        thisPt = myNode(0.0, 0.0);
        for (auto it=path[a].begin(); it!=path[a].end(); ++it) { // weighted addition
            thisPt.x += it->second * it->first->x;
            thisPt.y += it->second * it->first->y;
        }
        allPts.push_back (thisPt);
        
        test_search_problem.isSegmentFree (thisPt, lastPt, &tc); // non-uniform cost
        cost += tc; //sqrt ((thisPt.x-lastPt.x)*(thisPt.x-lastPt.x) + (thisPt.y-lastPt.y)*(thisPt.y-lastPt.y));
        #if _VIS
        ChartPoint cp = MapPoint(thisPt).toChartPoint(test_search_problem.pc_transform);
        //printf ("[%f,%f, %f]; ", cp.x, cp.y, cost);
        *(path_ostream_p) << "["<<cp.x<<","<<cp.y<<", "<<cost<<"]";
        if (a>0) *(path_ostream_p) << "; ";
        else *(path_ostream_p) << "];\n";
        cv::line (test_search_problem.image_to_display, 
                test_search_problem.cv_plot_coord(thisPt.x,thisPt.y), test_search_problem.cv_plot_coord(lastPt.x,lastPt.y), 
                            cvScalar(200.0,100.0,100.0),
                                    test_search_problem.line_thickness*test_search_problem.plot_scale);
        #endif
    }
    
    #if _VIS
    for (int a=allPts.size()-1; a>=0; --a)
        cv::circle (test_search_problem.image_to_display, 
                test_search_problem.cv_plot_coord(allPts[a].x,allPts[a].y),
                    test_search_problem.vertex_size*test_search_problem.plot_scale, cvScalar(150.0,0.0,255.0), -1, 8);
    
    cv::imshow("Display window", test_search_problem.image_to_display);
    #endif
    
    test_search_problem.clear();
    
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
    
    #if _VIS
    if (test_search_problem.save_image_interval != 0) {
        char imgFname[1024];
        sprintf(imgFname, "%s%s_path.png", test_search_problem.out_folderName.c_str(), 
                                                test_search_problem.imgPrefix.str().c_str());
        cv::imwrite(imgFname, test_search_problem.image_to_display);
    }
    printf ("\ncomputation time = %f, cost = %f\n", 0.0, cost);
    if (test_search_problem.pause_for_visualization)
        cvWaitKey();
    #endif
}

