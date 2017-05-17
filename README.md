```
** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal Search Library (DOSL)                                              *
*    A template-based C++ library for discrete search                                    *
*    Version 3.1                                                                         *
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
*************************************************************************************** **
```

Description: "Discrete Optimal Search Library (DOSL)" is a fast, efficient and easy-to-use library for construction of discrete representation (e.g., graph) and search (e.g., using algorithms like A-star, Dijkstra's, etc.) library written in C++, designed specifically for searching medium to large scale graphs for optimal paths. 

DOSL is designed to be:
* Fast (e.g., with integer coordinates for nodes but floating point cost as well as cost function needing to perform floating point operations online, and an average degree of the graph being 8, the library can expand about 150,000 nodes in the graph in just 1 second on a 1.8GHz processor machine with 8GB RAM.)
* Easy to use (being template-based, defining new arbitrary node-types, cost types, etc. is made easy. For graph connectivity, node accessibility tests, etc, user-defined classes can be used, which makes defining the graph structure very easy, yet highly flexible.) 

DOSL supports:
- Directed graphs with complex cost functions.
- In near future, other forms of discrete representations, such as simplicial complexes, will be supported.
- On the fly graph construction (i.e. no need to construct and store a complete graph before starting the search/planning process - makes it highly suitable for RRT-like graph construction).
- Arbitrary graph node type declaration.
- Arbitrary edge cost type declaration.
- Intermediate storage of paths during graph search.
- Multiple goals (or goal manifold) that determine when to stop search.
- Multiple start nodes for wave-front expansion type graph exploration.
- Event handling (i.e., call to user-defined functions upon generation, g-score updating or expansion of a node during the search process).
- Ability to write new planners with much ease. Comes with a weighted A-star (that includes Dijkstra's and normal A-star) planner by default.

Coming soon in near future:
* New planner with an implementation of the S* search algorithm (https://arxiv.org/abs/1607.07009)

NOTE: Discrete Optimal Search Library (DOSL) is a fork of
      the Yet Another Graph-Search Based Planning Library (YAGSBPL)
      hosted at https://github.com/subh83/YAGSBPL .
      YAGSBPL is now deprecated.

*******************************************************************************

Installation:
DOSL is (to a large extent) template-based.
There is nothing to build for the library itself.
Simply include the file "dosl/dosl" in your code.
You can (optionally) install the headers in the system folder by running
```
    sudo make install
```

Quick compilation of examples:
```
    make examples
```
Running the example programs:
```
    cd examples-dosl
    ./<program_name>
```
(note: you need to `cd` into the `examples-dosl` folder before running the executables)

*******************************************************************************

Output from example program `map2d_path_planning` showing the progress of A* search algorithm in finding shortest path in an 8-connected grid graph (requires OpenCV):

<p align="center"><table border=0 width=100%>
 <tr>
  <td><img src="http://subhrajit.net/files/externally-linked-files/images/github-DOSL/AStar8map2d_10000.png" width="200"/></td>
  <td><img src="http://subhrajit.net/files/externally-linked-files/images/github-DOSL/AStar8map2d_30000.png" width="200"/></td>
  <td><img src="http://subhrajit.net/files/externally-linked-files/images/github-DOSL/AStar8map2d__path.png" width="200"/></td>
 </tr>
 <tr>
  <td align="center">10000 vertices expanded</td>
  <td align="center">30000 vertices expanded</td>
  <td align="center">Final path</td>
 </tr>
</table></p>


******************************************************************************************

Basic Usage (with explanations):
-------------------------------
```C++
// standard headers
#include <stdio.h>
#include <vector>
#include <iostream>

// DOSL headers:
#define _DOSL_ALGORITHM AStar  // directs DOSL to load the header for 'AStar' planner only
#include <dosl/dosl>

// ==============================================================================
/* The following class defines the type of a vertex of the graph.
       Needs to be derived from DOSL-provided class template 'AStarNode<node_type,cost_type>' */

class myNode : public AStarNode<myNode,double>
{
public:
    int x, y; // (x,y) coordinates defining a vertex.
    
    myNode () { }
    myNode (int xx, int yy) : x(xx), y(yy) { }
    
    // The comparison operator must be defined for the node type.
    bool operator==(const myNode& n) const { return (x==n.x && y==n.y); }
    
    // An efficint hash function, 'getHashBin', for node type is desired, but is optional.
    int getHashBin (void) { return (abs(((int)x>>4) + ((int)y<<3) + ((int)x<<4) + ((int)y>>3))); }
    
    // optional printing function, 'print':
    void print (std::string head="", std::string tail="") const
            { _dosl_cout << head << "x=" << x << ", y=" << y << _dosl_endl; }
};

// ==============================================================================
/* The following class contains the description of the graph (graph connectivity)
       and search problem description (start and stop criteria).
       Needs to be derived from DOSL-provided class template 'AStarProblem<node_type,cost_type>' */

class searchProblem : public AStarProblem<myNode,double>
{
public:
    // user-defined problem parameters:
    myNode goal_node;
    // Constructors, if any
    searchProblem () { goal_node = myNode(150,100); }
    
    // -----------------------------------------------------------
    /* The following functions are use by the base class 'AStarProblem<myNode,double>' to determine 
           graph structure and search parameters. */
    
    /* Prototype for 'AStarProblem<>::getSuccessors':
         template <class nodeType, class costType> class AStarProblem {
             virtual void getSuccessors 
                 (NodeType& n, std::vector<NodeType>* const s, std::vector<CostType>* const s);
         }
       Description: Takes in a vertex, n, and returns its neighbors/successors, s, 
                    and the costs/distances of the edges, c. This defines graph connectivity. */
    
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
    
    /* Prototype for 'AStarProblem<>::getStartNodes':
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
    
    /* Prototype for 'AStarProblem<>::stopsearch':
        template <class nodeType, class costType> class AStarProblem 
            { virtual bool stopsearch (NodeType& n); }
       Description: Determines whether to stop the search when 'n' is being expanded.
       Optional -- If not provided, search will terminate only when heap is empty 
                   (all nodes in graph have been expanded). */
    
    bool stopsearch (myNode &n) {
        return (n==goal_node);
    }
};

// ==============================================================================

int main(int argc, char *argv[])
{
    searchProblem test_search_problem; // declare an instance of the search problem.
    test_search_problem.AllNodesSet.set_hash_table_size (8192); // optional parameters.
    test_search_problem.search(); // execute search.
    
    // get path from start to goal vertex.
    std::vector<myNode*> path = test_search_problem.getPointerPathToNode (test_search_problem.goal_node);
    
    // print path
    printf ("\nPath: \n[");
    for (int a=path.size()-1; a>=0; --a) {
        std::cout << "[" << path[a]->x << ", " << path[a]->y << "]";
        if (a>0) std::cout << "; ";
        else std::cout << "]\n\n";
    }
    
    return (1);
}
```


*******************************************************************************

Version history:
---------------

* May 2017: 3.1 released

* Discrete Optimal Search Library (DOSL) is a fork of
  the Yet Another Graph-Search Based Planning Library (YAGSBPL)
  hosted at https://github.com/subh83/YAGSBPL .
  YAGSBPL is now deprecated.

