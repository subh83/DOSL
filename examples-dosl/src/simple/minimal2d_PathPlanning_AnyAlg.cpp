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

// DOSL library
#define _DOSL_VERBOSE_LEVEL 1   // 0 by default
#include <dosl/dosl>

// Lical includes
#include <dosl/aux-utils/double_utils.hpp>  // SQRT2, isEqual_d

const double dists8connected[] = { SQRT2,  1.0,  SQRT2,
                                     1.0,  0.0,  1.0,
                                   SQRT2,  1.0,  SQRT2 };

// ==============================================================================

// The following class defines the type of a vertex of the graph:

template <class AlgClass>
class myNode : public AlgClass::template Node< myNode<AlgClass>, double > // e.g., 'AStar::Node<node_type,cost_type>'
{
public:
    double x, y;
    
    myNode () { }
    myNode (int xx, int yy) : x(xx), y(yy) { }
    
    // The comparison operator must be defined for the node type
    bool operator==(const myNode& n) const { return (x==n.x && y==n.y); }
    
    // An efficint hash function for node type is desired. optional.
    int getHashBin (void) { return (abs( (iround(x)>>4) + (iround(y)<<3) + (iround(x)<<4) + (iround(y)>>3))); }
    
    // optional print:
    void print (std::string head="", std::string tail="") const
            { _dosl_cout << head << "x=" << x << ", y=" << y << _dosl_endl; }
    
    
    // ---------------------------------------------------------------
    // SStar only: Operator overloading for convex combination [for path recnstruction in SStar algorithm.]
    
    myNode operator+(const myNode &b) const { // n1 + n2
        myNode ret = (*this);
        ret.x += b.x;
        ret.y += b.y;
        return (ret);
    }
    
    myNode operator*(const double &c) const { // n1 * c (right scalar multiplication)
        myNode ret = (*this);
        ret.x *= c;
        ret.y *= c;
        return (ret);
    } 
};

// ==============================================================================

// The following class contains the main search and path reconstruction function.

template <class AlgClass>
class searchProblem : public AlgClass::template Algorithm< searchProblem<AlgClass>, myNode<AlgClass>, double>
{
public:
    typedef myNode<AlgClass> MyNode;
    // user-defined problem parameters:
    MyNode goal_node;
    // Constructors, if any
    searchProblem () { goal_node = MyNode(200,150); }
    
    // -----------------------------------------------------------
    // The following functions are use by 'Algorithm' class to define graph sructure and search parameters
    
    /* Implementation of 'getSuccessors':
           virtual void getSuccessors (NodeType &n, std::vector<NodeType>* const s, std::vector<CostType>* const c);
       Takes in a vertex, n, and returns its neighbors/successors, s, and the costs/distances of the edges, c. */
    
    void getSuccessors (MyNode &n, std::vector<MyNode>* s, std::vector<double>* c) {
        // This function should account for obstacles, constraints and size of environment.
        MyNode tn;
        int ct=-1;
        for (int a=-1; a<=1; a++)
            for (int b=-1; b<=1; b++) {
                ++ct;
                if (a==0 && b==0) continue;
                
                if (this->algorithm_name() == "SStar") { /// triangulation for SStar
                    int xParity = ((int)round(fabs(n.x))) % 2;
                    if (xParity==0 && (a!=0 && b==-1)) continue;
                    if (xParity==1 && (a!=0 && b==1)) continue;
                }
                
                tn.x = n.x + a;
                tn.y = n.y + b;
                
                // Allow nodes to be contained within rectangle with corners (-100,-100) and (300,300)
                if (tn.x<-100 || tn.y<-100 || tn.x>300 || tn.y>300) continue;
                
                s->push_back (tn);
                c->push_back (dists8connected[ct]); 
            }
    }
    
    /* Implementation of 'getStartNodes':
           virtual std::vector<NodeType> getStartNodes (void);
       Returns the list of vertices(s) to start the search with. */
    
    std::vector<MyNode> getStartNodes (void) {
        std::vector<MyNode> startNodes;
        for (int a=0; a<1; ++a) {
            MyNode tn (0, 0); // start node
            startNodes.push_back (tn);
        }
        return (startNodes);
    }
    
    /* Implementation of 'stopSearch':
           virtual bool stopSearch (NodeType &n);
       Determines when to stop the search.
       Optional -- search will terminate only when heap is empty (all nodes in graph have been expanded). */
    
    bool stopSearch (MyNode &n) {
        return (n==goal_node);
    }
    
    // ---------------------------------------------------------------
    // ThetaStar only: 'isSegmentFree' checks and computes cost of long segments
    
    bool isSegmentFree (MyNode &n1, MyNode &n2, double* c)
    {
        double dx=n2.x-n1.x, dy=n2.y-n1.y;
        *c = sqrt(dx*dx + dy*dy);
        return (true); // this should check for intersection with obstacles
    }
};

// ==============================================================================
// Main function that runs search and prints result

template <class AlgClass>
void run_search (void) {
    std::cout << "run_search: Algorithm name: " << AlgClass::AlgorithmName << std::endl;
    
    searchProblem<AlgClass> test_search_problem;
    test_search_problem.all_nodes_set_p->reserve (8192); // optional.
    test_search_problem.search();
    
    // get path
    auto path = test_search_problem.reconstruct_pointer_path (test_search_problem.goal_node);
    
    // print path
    printf ("\nPath (as Octave array): \n[");
    for (int a=path.size()-1; a>=0; --a) {
        std::cout << "[" << path[a]->x << ", " << path[a]->y << "]";
        if (a>0) std::cout << "; ";
        else std::cout << "]\n\n";
    }
}

// ==============================================================================

int main(int argc, char *argv[])
{
    std::string alg_name;
    if (argc < 2) {
        std::cout << "Please enter algorithm name [" 
                        _YELLOW "AStar" YELLOW_ "|" _YELLOW "SStar" YELLOW_ "|" _YELLOW "ThetaStar" YELLOW_ "]: ";
        std::cin >> alg_name;
    }
    else
        alg_name = argv[1];
    
    if (alg_name == "AStar")
        run_search<AStar>();
    else if (alg_name == "SStar")
        run_search<SStar>();
    else if (alg_name == "ThetaStar")
        run_search<ThetaStar>();
    else
        std::cout << "Unknown algorithm." << std::endl;
    
    return (1);
}

