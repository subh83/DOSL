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

#ifdef DOSL_ALGORITHM_SStar
    #define COORD_TYPE double
#else
    #define COORD_TYPE int
#endif

class myNode : public DOSL_CLASS(Node)<myNode,double> // AStarNode<myNode,double>
{
public:
    COORD_TYPE x, y;
    
    myNode () { }
    myNode (COORD_TYPE xx, COORD_TYPE yy) : x(xx), y(yy) { }
    
    bool operator==(const myNode& n) const { return (x==n.x && y==n.y); } // This must be defined for the node
    void print (std::string head="", std::string tail="") const
            { _dosl_cout << _GREEN << head << "[" << this << "]" GREEN_ " x=" << x << ", y=" << y << _dosl_endl; }
    
    // Functions to be overwritten
    int getHashBin (void) { return (abs(((int)x>>4) + ((int)y<<3) + ((int)x<<4) + ((int)y>>3))); }
    
    
    #ifdef DOSL_ALGORITHM_SStar
    // + and * operators for convex combination of nodes when using SStar
    
    myNode operator+(const myNode &n) const {
        return (myNode(x+n.x,y+n.y));
    }
    
    myNode operator*(const double &c) const { // right scalar multiplication
        return (myNode(x*c,y*c));
    }
    #endif
};

// ==============================================================================

class searchProblem : public DOSL_CLASS(Problem)<myNode,double> // AStarProblem<myNode,double>
{
public:
    // Constructors, if any
    //searchProblem () { }
    
    // -----------------------------------------------------------
    
    void getSuccessors (myNode &n, std::vector<myNode>* s, std::vector<double>* c) {
        // This function should account for obstacles and size of environment.
        myNode tn;
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
                // TODO
                s->push_back(tn);
                c->push_back(sqrt((double)(a*a+b*b))); 
            }
    }
    
    // -----------------------------------------------------------
    
    std::vector<myNode> getStartNodes (void) {
        std::vector<myNode> startNodes;
        for (int a=0; a<1; ++a) {
            myNode tn (0, 0); // start node
            startNodes.push_back (tn);
        }
        return (startNodes);
    }
    
    // -----------------------------------------------------------
    
    bool stopsearch (myNode &n) {
        return (n.x==150 && n.y==100);
    }
};

// ==============================================================================

int main(int argc, char *argv[])
{
    searchProblem test_search_problem;
    test_search_problem.AllNodesSet.set_hash_table_size (8192);
    test_search_problem.search();
    
    // get path
    std::vector<myNode*> path = test_search_problem.getPointerPathToNode (myNode(150,100));
    
    printf ("\nPath: [ ");
    for (int a=path.size()-1; a>=0; --a) {
        std::cout << "[" << path[a]->x << ", " << path[a]->y << "]";
        if (a>0) std::cout << "; ";
        else std::cout << " ]\n";
    }
    
    return (1);
}

