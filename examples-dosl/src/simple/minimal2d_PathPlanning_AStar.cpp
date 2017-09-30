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

// DOSL library
#define _DOSL_VERBOSE_LEVEL 1   // 0 by default
#include <dosl/dosl>

// Lical includes
#include <dosl/aux-utils/double_utils.hpp>  // SQRT2

const double dists8connected[] = {SQRT2, 1, SQRT2, 1, 1, SQRT2, 1, SQRT2};

// ==============================================================================

// The following class defines the type of a vertex of the graph:

class myNode : public AStar::Node<myNode,double> // CRTP
{
public:
    int x, y;
    
    myNode () { }
    myNode (int xx, int yy) : x(xx), y(yy) { }
    
    // The comparison operator must be defined for the node type
    bool operator==(const myNode& n) const { return (x==n.x && y==n.y); }
    
    // An efficint hash function for node type is desired. optional.
    int getHashBin (void) { return (abs((x>>4)+(y<<3)+(x<<4)+(y>>3))); }
    
    // optional print:
    void print (std::string head="", std::string tail="") const
            { _dosl_cout << head << "x=" << x << ", y=" << y << _dosl_endl; }
    
};

// ==============================================================================

// The following class contains the main search and path reconstruction function.

class searchProblem : public AStar::Algorithm<myNode,double>
{
public:
    // user-defined problem parameters:
    myNode goal_node;
    // Constructors, if any
    searchProblem () { goal_node = myNode(200,150); }
    
    // -----------------------------------------------------------
    // The following functions are use by 'Algorithm' class to define graph sructure and search parameters
    
    /* Implementation of 'getSuccessors':
           virtual void getSuccessors (NodeType &n, std::vector<NodeType>* const s, std::vector<CostType>* const c);
       Takes in a vertex, n, and returns its neighbors/successors, s, and the costs/distances of the edges, c. */
    
    void getSuccessors (myNode &n, std::vector<myNode>* s, std::vector<double>* c) {
        // This function should account for obstacles, constraints and size of environment.
        myNode tn;
        int ct=0;
        for (int a=-1; a<=1; a++)
            for (int b=-1; b<=1; b++) {
                if (a==0 && b==0) continue;
                
                tn.x = n.x + a;
                tn.y = n.y + b;
                
                // A single disk-shaped obstacle of radius 40 centered at (110,80)
                int rx=tn.x-110, ry=tn.y-80;
                if (rx*rx + ry*ry < 40*40) continue;
                
                s->push_back (tn);
                c->push_back (dists8connected[ct]); 
                ++ct;
            }
    }
    
    /* Implementation of 'getStartNodes':
           virtual std::vector<NodeType> getStartNodes (void);
       Returns the list of vertices(s) to start the search with. */
    
    std::vector<myNode> getStartNodes (void) {
        std::vector<myNode> startNodes;
        for (int a=0; a<1; ++a) {
            myNode tn (0, 0); // start node
            startNodes.push_back (tn);
        }
        return (startNodes);
    }
    
    /* Implementation of 'stopSearch':
           virtual bool stopSearch (NodeType &n);
       Determines when to stop the search.
       Optional -- search will terminate only when heap is empty (all nodes in graph have been expanded). */
    
    bool stopSearch (myNode &n) {
        return (n==goal_node);
    }
    
};


// ==============================================================================

int main(int argc, char *argv[])
{
    searchProblem my_search_problem;
    my_search_problem.AllNodesSet.set_hash_table_size (8192); // optional.
    my_search_problem.search();
    
    // get path
    std::vector<myNode*>  path = my_search_problem.reconstructPointerPath (my_search_problem.goal_node);
    
    // print path
    printf ("\nPath (as Octave array): \npp = [");
    for (int a=path.size()-1; a>=0; --a) {
        std::cout << "[" << path[a]->x << ", " << path[a]->y << "]";
        if (a>0) std::cout << "; ";
        else std::cout << "]\n\n";
    }
    // Octave plot command: figure; axis equal; hold on; plot(pp(:,1),pp(:,2), 'r-', 'LineWidth',2.0); c=[110,80]; r=40; fill(c(1)+r*cos(0:pi/20:2*pi),c(2)+r*sin(0:pi/20:2*pi), 'k');
    
    return (1);
}

