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
#include <vector>

// Other libraries:
// Open CV:
/*#include <opencv/cv.h>
#include <opencv/cvaux.h>
#include <opencv/highgui.h> */
// Armadillo
/* #include <armadillo> */

// DOSL library

// Optional parameters -- to be declared before including dosl: 
// #define _DOSL_DEBUG 1                // 0,1, or 2
// #define _DOSL_AUTOCORRECT 1          // 0 or 1
#define _DOSL_VERBOSE_LEVEL 1        // 0 to 5
// #define _DOSL_VERBOSE_ITEMS "NONE"   // "fun1,fun2" or "NONE" or "ALL"
           // Common items: "reconstructPath,findCamefromPoint,getAllAttachedMaximalSimplices,getAllMaximalSimplicesFromSet"
// #define _DOSL_PRINT_COLORS 1
// #define _DOSL_EVENTHANDLER 1

#ifndef _DOSL_ALGORITHM // can pass at command line during compilation: -D_DOSL_ALGORITHM=AStar
    #define _DOSL_ALGORITHM  AStar
#endif
#include <dosl/dosl>


// ==============================================================================
// change these:
#define GRAPH_TYPE 6 // 6 (equilateral triangular grid)  or  8 (uniform square grid)

// ==============================================================================

#if defined(DOSL_ALGORITHM_SStar) || GRAPH_TYPE==6
    #define COORD_TYPE double
#else
    #define COORD_TYPE int
#endif

#define PIBY3                 1.0471975512
#define SQRT3BY2              0.86602540378
#define INFINITESIMAL_DOUBLE  1e-6


class myNode : public DOSL_CLASS(Node)<myNode,double> // AStarNode<myNode,double>
{
public:
    COORD_TYPE x, y;
    
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
    
    // constructors
    myNode () { }
    myNode (COORD_TYPE xx, COORD_TYPE yy) : x(xx), y(yy) { put_in_grid(); }
    
    // Comparison operator must be defined for the node
    bool operator==(const myNode& n) const {
        #if defined(DOSL_ALGORITHM_SStar) || GRAPH_TYPE==6
        return (fabs(x-n.x)<INFINITESIMAL_DOUBLE  &&  fabs(y-n.y)<INFINITESIMAL_DOUBLE);
        #else
        return ((x==n.x) && (y==n.y));
        #endif
    }
    
    void print (std::string head="", std::string tail="") {
        _dosl_cout << _GREEN << head << "[" << this << "]" GREEN_ << "["<< getHashBin() << "]" << " x=" << x << ", y=" << y << _dosl_endl;
        _dosl_cout << "Successors: " << Successors << _dosl_endl; // Successors is an 'unordered_map'
    }
    
    // Functions being overwritten
    int getHashBin (void) {
        #if GRAPH_TYPE==6
        /*int xi=((int)round(x/0.5)), yi=((int)round(y/SQRT3BY2));
        return (abs((xi>>4) + (yi<<3) + (xi<<4) + (yi>>3)));*/
        return ( MAX(round(fabs(x)+INFINITESIMAL_DOUBLE), round(fabs(x)-INFINITESIMAL_DOUBLE)) );
        #elif defined(DOSL_ALGORITHM_SStar)
        int xi=((int)round(x)), yi=((int)round(y));
        return (abs((xi>>4) + (yi<<3) + (xi<<4) + (yi>>3)));
        #else
        return (abs((x>>4) + (y<<3) + (x<<4) + (y>>3)));
        #endif
    }
    
    
    #ifdef DOSL_ALGORITHM_SStar
    // + and * operators for convex combination of nodes when using SStar
    
    myNode operator+(const myNode &n) const {
        myNode tn;
        tn.x=x+n.x; tn.y=y+n.y;
        return (tn);
    }
    
    myNode operator*(const double &c) const { // right scalar multiplication
        myNode tn;
        tn.x=x*c; tn.y=y*c;
        return (tn);
    }
    #endif
};

// ==============================================================================

class searchProblem : public DOSL_CLASS(Problem)<myNode,double> // AStarProblem<myNode,double>
{
public:
    // Constructors, if any
    //searchProblem () { }
    
    // user-defined variables
    myNode goal;
    
    // User-defined functions:
    bool isAccessible (myNode &n) { // used inside 'getSuccessor'
        // defines an obstacles of radius 20 centered at (0,0)
        return ( n.x*n.x + n.y*n.y > 400 );
    }
    
    // ==============================================================================
    // The following functions are use by 'AStarProblem<nodeType,costType>' class to define graph sructure and search parameters
    
    // -----------------------------------------------------------
    
    void getSuccessors (myNode &n, std::vector<myNode>* s, std::vector<double>* c) {
        // This function should account for obstacles and size of environment.
        myNode tn;
        // if (~isAccessible(n)) return;
        
        #if GRAPH_TYPE==8 // uniform square grid
        for (int a=-1; a<=1; a++)
            for (int b=-1; b<=1; b++) {
                if (a==0 && b==0) continue;
                
                #ifdef DOSL_ALGORITHM_SStar // to avoid degenerate 3-simplices when using SStar
                int xParity = ((int)fabs(n.x)) % 2;
                if (xParity==0 && (a!=0 && b==-1)) continue;
                if (xParity==1 && (a!=0 && b==1)) continue;
                #endif
                
                tn.x = n.x + a;
                tn.y = n.y + b;
                
                if (isAccessible(tn)) {
                    s->push_back(tn);
                    c->push_back(sqrt((double)(a*a+b*b)));
                }
            }
        
        #elif GRAPH_TYPE==6 // equilateral triangular grid
        double th;
        for (int a=0; a<6; ++a) {
            th = a * PIBY3;
            tn.x = n.x + 1.0*cos(th);
            tn.y = n.y + 1.0*sin(th);
            
            if (isAccessible(tn)) {
                s->push_back(tn);
                c->push_back(1.0);
            }
        }
        
        #endif
    }
    
    // -----------------------------------------------------------
    
    std::vector<myNode> getStartNodes (void) {
        std::vector<myNode> startNodes;
        for (int a=0; a<1; ++a) {
            myNode tn (-70, -40); // start node
            startNodes.push_back (tn);
        }
        return (startNodes);
    }
    
    // -----------------------------------------------------------
    
    bool stopsearch (myNode &n) {
        return (n==goal);
    }
};

// ==============================================================================

int main(int argc, char *argv[])
{
    searchProblem test_search_problem;
    test_search_problem.goal = myNode(20,30);
    test_search_problem.search();
    
    // get path
    std::vector<myNode*> path = test_search_problem.getPointerPathToNode (test_search_problem.goal);
    
    printf ("\nPath: [ ");
    for (int a=path.size()-1; a>=0; --a) {
        std::cout << "[" << path[a]->x << ", " << path[a]->y << "]";
        if (a>0) std::cout << "; ";
        else std::cout << " ]\n";
    }
    
    return (1);
}

