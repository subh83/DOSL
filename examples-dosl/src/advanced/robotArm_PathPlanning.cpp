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
#include <vector>

// Open CV:
#include <opencv2/opencv.hpp> 
#include <opencv2/highgui.hpp> 

// DOSL library
#define _DOSL_VERBOSE_LEVEL 1        // 0 to 5
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

#define SMALL_EPS INFINITESIMAL_DOUBLE
#define _VIS 2 // bit 1: animation, bit 2: final path

// ==============================================================================

// Helper function

double vec_dist (std::vector<double> p1, std::vector<double> p2) {
    double distSq=0.0, delta;
    for (int a=0; a<p1.size(); ++a) {
        delta = p1[a] - p2[a];
        distSq += delta*delta;
    }
    return (sqrt (distSq));
}

// --------------------------------------------------------------

// A node of the graph
class ArmConfig : public _DOSL_ALGORITHM::Node<ArmConfig,double>
{
public:
    // main data
    std::vector<double> th, ef;
    
    void check_modulo (void) {
        for (int a=0; a<th.size(); ++a) {
            if (isEqual_d(th[a],2.0*PI)) th[a] = 0.0;
            else if (isEqual_d(th[a],0.0)) th[a] = 0.0;
            else th[a] -= 2.0*PI * floor (th[a] / (2.0*PI)); // brings it in [0.0, 2*PI)
        }
    }
    
    void put_in_grid (void) { check_modulo(); }
    
    // This must be defined for the node
    bool operator==(const ArmConfig& n) const {
        for (int a=0; a<th.size(); ++a)
            if ( !isEqual_d (th[a], n.th[a]) )
                return (false);
        return (true);
    }
    
    // constructor
    ArmConfig () { th.resize(3,0); }
    ArmConfig (std::vector<double> tth) : th(tth) {
        if (tth.size()!=3) {
            printf("WARNING: th.size()!=3");
            print();
        }
        put_in_grid();
    }
    
    // Inherited functions being overwritten
    int GetHashBin (void) const {
        unsigned int bin = 0;
        for (int a=0; a<th.size(); ++a)
            bin += ( ((unsigned int)round(th[a]*100.0)) << (((unsigned int)round(th[a]*100.0)) % 8) ) + 
                                ( ((unsigned int)round(th[a]*100.0)) >> (((unsigned int)round(th[a]*100.0)) % 5) );
        return (bin);
    }
    
    // print
    void print (std::string head="", std::string tail="") const {
        _dosl_cout << _GREEN + head << " (" << this << ")" GREEN_ "th = {" 
                    << th[0] <<","<< th[1] <<","<< th[2] << "}. hash bin = " << GetHashBin();
        if(ef.size()>0) { _dosl_cout << ". ef = {" << ef[0] <<","<< ef[1] << "}. "; }
        _dosl_cout << _dosl_endl;
    }
};

ArmConfig weighted_addition (const std::unordered_map<ArmConfig*,double>& configP_weight) {
    ArmConfig* refConfig_p = configP_weight.begin()->first;
    ArmConfig retConfig = ArmConfig( std::vector<double>(refConfig_p->th.size(), 0.0) );
    double thisTh;
    
    for (auto it=configP_weight.begin(); it!=configP_weight.end(); ++it) {
        for (int b=0; b<retConfig.th.size(); ++b) {
            thisTh = it->first->th[b];
            if ( it!=configP_weight.begin()  &&  
                    fabs(refConfig_p->th[b] - thisTh) > PI + SMALL_EPS ) { // local chart
                if (thisTh > refConfig_p->th[b])  thisTh -= 2*PI;
                else  thisTh += 2*PI;
            }
            retConfig.th[b] += it->second * thisTh;
        }
    }
    
    return (retConfig);
}


// ==============================================================================

//#ifdef S_STAR_ALGORITHM
double edges[][3] = { {-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}, {-1.0, -1.0, -1.0}, 
                        {1.0, 1.0, 1.0}, {1.0, -1.0, 0.0}, {-1.0, 1.0, 0.0}, {0.0, -1.0, -1.0}, {0.0, 1.0, 1.0}, 
                        {-1.0, 0.0, -1.0}, {1.0, 0.0, 1.0}, {0.0, 0.0, -1.0}, {0.0, 0.0, 1.0} };

class SearchProblem : public _DOSL_ALGORITHM::Algorithm<SearchProblem,ArmConfig,double>
{
public:
    RSJresource expt_container;
    
    //std::vector<double> MAXs, MINs;
    cvParseMap2d obs_map;
    std::vector<double> seg_lengths, normalized_seg_lengths; // pixel units
    std::vector<double> arm_base; // pixel units
    
    cv::Mat image_to_display;
    
    ArmConfig start_config, final_config;
    std::vector<double> goal_ef;
    
    double arm_collison_check_segment_dl; // 1.0 / 2 (fractions of 1.0)
    double graph_edge_check_interp_dl; // 1.0 / 2 (fractions of 1.0)
    std::vector<double> theta_step;
    double end_effector_target_thresh;
    
    // Constructor
    SearchProblem (std::string expt_fName="../files/expt/robot_arm_expt.json", std::string expt_name="path1")
    {
        // TODO: Read from file
        std::ifstream my_fstream (expt_fName);
        expt_container = RSJresource (my_fstream)[expt_name];
        std::unordered_map<std::string,double> sym_vals;
        
        obs_map = cvParseMap2d (get_path(expt_fName)+expt_container["obstacle_map"].as<std::string>("arm_workspace.png"), 0, 200);
        sym_vals["width"] = (double)(obs_map.width()); sym_vals["height"] = (double)(obs_map.height());
        
        seg_lengths = std::vector<double>({200.0, 50.0, 100.0}); // pixel units. Must be integers!
        //sym_vals["l0"] = seg_lengths[0]; sym_vals["l1"] = seg_lengths[1]; sym_vals["l2"] = seg_lengths[2];
        double lsum = 0.0;
        for (int a=0; a<seg_lengths.size(); ++a) {
            sym_vals["l"+std::to_string(a)] = seg_lengths[a]; // l0, l1, ...
            lsum += seg_lengths[a];
        }
        sym_vals["lsum"] = lsum;
        for (int a=0; a<seg_lengths.size(); ++a)
            normalized_seg_lengths.push_back (seg_lengths[a]/lsum);
        
        arm_base = std::vector<double>( {  
                        MathEvaluator(expt_container["arm_base"][0].as<std::string>("width/2"), sym_vals).eval(), 
                        MathEvaluator(expt_container["arm_base"][1].as<std::string>("height/2"), sym_vals).eval()    } );
        
        arm_collison_check_segment_dl = expt_container["search_options"]["arm_collison_check_segment_dl"].as<double>(1.0);
        graph_edge_check_interp_dl = 10.0;
        end_effector_target_thresh = expt_container["search_options"]["end_effector_target_thresh"].as<double>(5.0);
        
        for (int a=0; a<seg_lengths.size(); ++a) {
            theta_step.push_back ( MathEvaluator(expt_container["dth"][a].as<std::string>("2.0*pi/20.0"), sym_vals).eval() );
            sym_vals["dth"+std::to_string(a)] = theta_step[a]; // dth0, dth1, ...
        }
        
        std::vector<double> thsInitial;
        
        for (int a=0; a<seg_lengths.size(); ++a)
            thsInitial.push_back ( MathEvaluator(expt_container["start_config"][a].as<std::string>("0.0"), sym_vals).eval() );
        start_config = ArmConfig (thsInitial); // ef = {140.0, 0.0}
        
        goal_ef = std::vector<double> ({ // pixel units
                        MathEvaluator(expt_container["goal_ef"][0].as<std::string>("lsum"), sym_vals).eval(), 
                        MathEvaluator(expt_container["goal_ef"][1].as<std::string>("0.0"), sym_vals).eval()    }); 
        
        #if _VIS
        image_to_display = obs_map.getCvMat (COLOR_MAP, ORIGINAL_MAP);
        cv::circle (image_to_display, cvPoint(goal_ef[0]+arm_base[0],goal_ef[1]+arm_base[1]), 
                                                            end_effector_target_thresh, cvScalar(100.0,100.0,200.0), -1, 8);
        cv::namedWindow("Display_window", cv::WINDOW_AUTOSIZE);
        cv::imshow("Display_window", image_to_display);
        cv::waitKey(10);
        #endif
        
        // Set planner variables
        all_nodes_set_p->reserve (16384);
        progress_show_interval = 100;
    }
    
    // -----------------------------------------------------------
    
    #if _VIS
    void drawArm (ArmConfig& tn, CvScalar arm_col = cvScalar(200.0,100.0,100.0), CvScalar joint_col = cvScalar(50.0,50.0,50.0)) {
        std::vector<double> last_arm_point, arm_point = arm_base;
        std::vector< std::vector<double> > arm_points;
        double th_sum = 0.0;
        for (int a=0; a<tn.th.size(); ++a) {
            arm_points.push_back (arm_point);
            th_sum += tn.th[a];
            last_arm_point = arm_point;
            arm_point[0] += seg_lengths[a]*cos(th_sum);
            arm_point[1] += seg_lengths[a]*sin(th_sum);
            cv::line (image_to_display, cvPoint(last_arm_point[0], last_arm_point[1]), cvPoint(arm_point[0], arm_point[1]), 
                                        arm_col, 2.0);
        }
        for (int a=0; a<arm_points.size(); ++a)
            cv::circle (image_to_display, cvPoint(arm_points[a][0],arm_points[a][1]), 2.0, joint_col, -1, 8);
    }
    #endif
    
    bool isPointFree (std::vector<double> arm_point) {
        for (double a=-SMALL_EPS; a<=SMALL_EPS; a+=SMALL_EPS)
            for (double b=-SMALL_EPS; b<=SMALL_EPS; b+=SMALL_EPS) {
                if ( !obs_map.isFree ( (int)round(arm_point[0]+arm_base[0]+a), (int)round(arm_point[1]+arm_base[1]+b) ) )
                    return (false);
            }
        return (true);
    }
    
    // ------------------------------------------------------------------------
    
    
    bool isConfigurationValid (ArmConfig& tn, std::vector<double>* end_effector=NULL) {
        std::vector<double> arm_point = std::vector<double>({0.0, 0.0});
        double dl, dl_cos_th_sum, dl_sin_th_sum, th_sum=0.0;
        for (int a=0; a<tn.th.size(); ++a) {
            dl = arm_collison_check_segment_dl;
            th_sum += tn.th[a]; dl_cos_th_sum = dl*cos(th_sum); dl_sin_th_sum = dl*sin(th_sum);
            for (double l=0.0; !isGreater_d(l,seg_lengths[a]); l+=dl) {
                arm_point[0] += dl_cos_th_sum;
                arm_point[1] += dl_sin_th_sum;
                if (!isPointFree(arm_point)) return (false);
            }
        }
        tn.ef = arm_point;
        if (end_effector) *end_effector = arm_point;
        return (true);
    }
    
    // ------------------------------------------------------------------------
    
    bool isEdgeAccessible (ArmConfig& tn1, ArmConfig& tn2) { // works only for 8-connected grid!!
        std::vector<double> ef1, ef2, dths;
        
        if ( (!isConfigurationValid(tn1, &ef1)) || !(isConfigurationValid(tn2, &ef2)) )
            return (false);
        
        return (true);
    }
    
    // -----------------------------------------------------------
    
    
    
    void getSuccessors (ArmConfig &n, std::vector<ArmConfig>* s, std::vector<double>* c) // *** This must be defined
    {
        // This function should account for obstacles and size of environment.
        ArmConfig tn, ttn;
        double dth;
        
        for (int a=0; a<14; ++a) {
            tn = n;
            double energy = 0.0;
            double lengthSqSum = 0.0;
            for (int b=tn.th.size()-1; b>=0; --b) {
                dth = theta_step[b]*edges[a][b];
                tn.th[b] = tn.th[b] + dth;
                lengthSqSum += seg_lengths[b] * seg_lengths[b];
                energy += dth*dth*lengthSqSum; //dth*dth*seg_lengths[b]; // dth*dth
            }
            
            ttn = tn;
            tn.check_modulo();
            if (!isEdgeAccessible(tn,n)) continue;
            
            s->push_back (tn);
            c->push_back (sqrt(energy)); 
        }
        
    }
    
    // -----------------------------------------------------------
    
    double GetHeuristics (ArmConfig& n)
    {
        /* double dx = goal.x - n.x;
        double dy = goal.y - n.y;
        return (sqrt(dx*dx + dy*dy)); */
        return (0.0);
    }
    
    // -----------------------------------------------------------
    
    std::vector<ArmConfig> getStartNodes (void) 
    {
        std::vector<ArmConfig> starts;
        
        starts.push_back (start_config);
        
        return (starts);
    }
    
    // -----------------------------------------------------------
    
    void nodeEvent (ArmConfig &n, unsigned int e) 
    {
        
        if (e & EXPANDED) {
            //n.print("expanding: ");
            #if (_VIS & 1)
            drawArm (n);
            cv::imshow("Display_window", image_to_display);
            cv::waitKey(1);
            #endif
        }
        
        if (e & ERROR) {
            #if (_VIS & 1)
            drawArm (n, cvScalar(80.0,80.0,200.0));
            cv::imshow("Display_window", image_to_display);
            cv::waitKey();
            #endif
        }
    }
    
    // ---------------------------------------
    
    bool stopSearch (ArmConfig &n) {
        if (n.ef.size() == 0)  isConfigurationValid (n);
        
        if (vec_dist (goal_ef, n.ef) < end_effector_target_thresh) {
            final_config = n;
            return (true);
        }
        return (false);
    }
};

// ==============================================================================

int main(int argc, char *argv[])
{
    //DOSL_RUNTIME_VERBOSE_SWITCH = 0;
    compute_program_path();
    
    std::string expt_fName = program_path+"../files/expt/robot_arm_expt.json";
    std::string expt_name = "path1_fast";
    if (argc == 2) {
        expt_name = argv[1];
    }
    else if (argc > 2) {
        expt_fName = argv[1];
        expt_name = argv[2];
    }
    
    printf (_BOLD _YELLOW "Note: " YELLOW_ BOLD_ "Using algorithm " _YELLOW  MAKESTR(_DOSL_ALGORITHM)  YELLOW_ 
                ". Run 'make' to recompile with a different algorithm.\n" 
                "Note: This program currently does not support 'ThetaStar' algorithm.\n");
    
    SearchProblem test_search_problem (expt_fName, expt_name);
    std::cout << "before search..." << std::endl;
    test_search_problem.search();
    std::cout << "after search..." << std::endl;
    
    // Get path
    DOSL_RUNTIME_VERBOSE_SWITCH = 1;
    auto path = test_search_problem.reconstruct_weighted_pointer_path (test_search_problem.final_config);
    
    // -------------------------------------
    
    // Print and draw path
    printf("\nPath: ");
    ArmConfig thisPt = ArmConfig (test_search_problem.start_config);
    ArmConfig lastPt, interpPt;
    std::vector<ArmConfig> allPts;
    double cost = 0.0;
    
    for (int a=path.size()-1; a>=0; --a) {
        lastPt = thisPt;
        thisPt = weighted_addition (path[a]);
        allPts.push_back (thisPt);
        printf ("[%f,%f,%f]; ", (double)thisPt.th[0], (double)thisPt.th[1], (double)thisPt.th[2]);
        #if (_VIS & 2)
        int interp_steps = 1;
        for (int b=0; b<interp_steps; ++b) {
            // TODO
            double darkness = 50.0 + 200.0 * (1.0 - ((double)a)/((double)path.size()-1.0) );
            test_search_problem.drawArm (thisPt, cvScalar(255.0, 255.0-darkness, 255.0-0.5*darkness), 
                                                    cvScalar(255.0-0.4*darkness, 255.0-0.4*darkness, 255.0-0.4*darkness) );
        }
        #endif
    }
    printf ("\nLast point: [%d*dth0, %d*dth1, %d*dth2]\n", (int)round(thisPt.th[0]/test_search_problem.theta_step[0]), 
                                                    (int)round(thisPt.th[1]/test_search_problem.theta_step[1]),
                                                    (int)round(thisPt.th[2]/test_search_problem.theta_step[2]) );
    #if _VIS
    cv::imshow("Display_window", test_search_problem.image_to_display);
    if (true) {
        char imgFname[1024];
        sprintf(imgFname, (program_path+"../files/out/robot_arm_%s.png").c_str(), expt_name.c_str());
        cv::imwrite(imgFname, test_search_problem.image_to_display);
    }
    cv::waitKey();
    #endif
    
}
