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
#define _DOSL_ALGORITHM AStar  // directs DOSL to load the header for 'AStar' planner only
#include <dosl/dosl>

// ==============================================================================

// The following class defines the type of a vertex of the graph:

class myNode : public AStarNode<myNode,double> // DOSL provides 'AStarNode<node_type,cost_type>'
{
public:
    int x, y;
    
    myNode () { }
    myNode (int xx, int yy) : x(xx), y(yy) { }
    
    // The comparison operator must be defined for the node type
    bool operator==(const myNode& n) const { return (x==n.x && y==n.y); }
    
    // An efficint hash function for node type is desired. optional.
    int getHashBin (void) { return (abs(((int)x>>4) + ((int)y<<3) + ((int)x<<4) + ((int)y>>3))); }
    
    // optional:
    void print (std::string head="", std::string tail="") const
            { _dosl_cout << head << "x=" << x << ", y=" << y << _dosl_endl; }
};

// ==============================================================================

// The following class contains the main search and path reconstruction function.

class searchProblem : public AStarProblem<myNode,double> // DOSL provides 'AStarProblem<node_type,cost_type>'
{
public:
    // user-defined problem parameters:
    myNode goal_node;
    // Constructors, if any
    searchProblem () { goal_node = myNode(200,150); }
    
    // -----------------------------------------------------------
    // The following functions are use by 'AStarProblem<nodeType,costType>' class to define graph sructure and search parameters
    
    /* Implementation of 'getSuccessors':
        template <class nodeType, class costType> class AStarProblem 
            { virtual void getSuccessors (NodeType &n, std::vector<NodeType>* const s, std::vector<CostType>* const c); }
       Takes in a vertex, n, and returns its neighbors/successors, s, and the costs/distances of the edges, c. */
    
    void getSuccessors (myNode &n, std::vector<myNode>* s, std::vector<double>* c) {
        // This function should account for obstacles, constraints and size of environment.
        myNode tn;
        for (int a=-1; a<=1; a++)
            for (int b=-1; b<=1; b++) {
                if (a==0 && b==0) continue;
                
                tn.x = n.x + a;
                tn.y = n.y + b;
                
                s->push_back(tn);
                c->push_back(sqrt((double)(a*a+b*b))); 
            }
    }
    
    /* Implementation of 'getStartNodes':
        template <class nodeType, class costType> class AStarProblem 
            { virtual std::vector<NodeType> getStartNodes (void); }
       Returns the list of vertices(s) to start the search with. */
    
    std::vector<myNode> getStartNodes (void) {
        std::vector<myNode> startNodes;
        for (int a=0; a<1; ++a) {
            myNode tn (0, 0); // start node
            startNodes.push_back (tn);
        }
        return (startNodes);
    }
    
    /* Implementation of 'stopsearch':
        template <class nodeType, class costType> class AStarProblem 
            { virtual bool stopsearch (NodeType &n); }
       Determines when to stop the search.
       Optional -- search will terminate only when heap is empty (all nodes in graph have been expanded). */
    
    bool stopsearch (myNode &n) {
        return (n==goal_node);
    }
};

// ==============================================================================

int main(int argc, char *argv[])
{
    searchProblem test_search_problem;
    test_search_problem.AllNodesSet.set_hash_table_size (8192); // optional.
    test_search_problem.search();
    
    // get path
    std::vector<myNode*> path = test_search_problem.getPointerPathToNode (test_search_problem.goal_node);
    
    // print path
    printf ("\nPath (as Octave array): \n[");
    for (int a=path.size()-1; a>=0; --a) {
        std::cout << "[" << path[a]->x << ", " << path[a]->y << "]";
        if (a>0) std::cout << "; ";
        else std::cout << "]\n\n";
    }
    
    return (1);
}

