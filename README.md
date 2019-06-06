<!--
/** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal Search Library (DOSL)                                              *
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
 -->

Discrete Optimal Search Library (DOSL)
--------------------------------------

### Description:
"Discrete Optimal Search Library (DOSL)" is a fast, efficient and easy-to-use library for construction of discrete representation (e.g., a graph or a simplicial complex) and search (e.g., using algorithms like A-star, Dijkstra's, etc.) library written in C++, designed specifically for searching medium to large scale graphs or simplicial complexes for optimal paths. 

*******************************************************************************

Output from example program `map2d_PathPlanning` showing the progress of A* search algorithm in finding shortest path in an 8-connected grid graph (requires OpenCV):

<p align="center"><table border=0 width=100%>
 <tr>
  <td><img src="http://subhrajit.net/files/externally-linked-files/images/github-DOSL/AStar8map2d_10000.png" width="180"/></td>
  <td><img src="http://subhrajit.net/files/externally-linked-files/images/github-DOSL/AStar8map2d_20000.png" width="180"/></td>
  <td><img src="http://subhrajit.net/files/externally-linked-files/images/github-DOSL/AStar8map2d_30000.png" width="180"/></td>
  <td><img src="http://subhrajit.net/files/externally-linked-files/images/github-DOSL/AStar8map2d__path.png" width="180"/></td>
 </tr>
 <tr>
  <td align="center">10000 vertices expanded</td>
  <td align="center">20000 vertices expanded</td>
  <td align="center">30000 vertices expanded</td>
  <td align="center">Final path</td>
 </tr>
</table></p>


See the [Overview of Selected Example Programs](https://github.com/subh83/DOSL/wiki/Overview-of-Selected-Example-Programs) for more sample outputs from the provided example programs.

*******************************************************************************

### DOSL is designed to be:
* Fast (e.g., with integer coordinates for nodes but floating point cost as well as cost function needing to perform floating point operations online, and an average degree of the graph being 8, the A* search algorithm of the library can expand about 150,000 nodes in the graph in about 0.53 second on a 1.8GHz quad-core processor machine with 16GB RAM.)
* Easy to use / versatile (being template-based, defining new arbitrary node-types, cost types, etc. is made easy. For graph connectivity, node accessibility tests, etc, user-defined classes can be used, which makes defining the graph structure very easy, yet highly flexible.) 

### DOSL supports:
- Directed graphs with complex cost functions.
- On the fly graph construction (i.e. no need to construct and store a complete graph before starting the search/planning process - makes it highly suitable for RRT-like graph construction).
- Arbitrary graph node type declaration.
- Arbitrary edge cost type declaration.
- Intermediate storage of paths during graph search.
- Multiple goals (or goal manifold) that determine when to stop search.
- Multiple start nodes for wave-front expansion type graph exploration.
- Event handling (i.e., call to user-defined functions upon generation, g-score updating or expansion of a node during the search process).
- Other forms of discrete representations, such as simplicial complexes (S-star algorithm), is supported by specific planners.
- Ability to write new planners with ease. Comes with a weighted A-star (that includes Dijkstra's and normal A-star), Theta-star and S-star planner by default.

### New:
* Backward compatibility issue (user-end change since v3.3): The `[AlgName]::Algorithm` class template now takes in the derived class as its first template parameter (a [CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern)). Thus, definition of a `searchProblem` class should now be made as follows:
```C++
    class searchProblem : public AStar::Algorithm <searchProblem, myNode, double>
    { /* ... */};                                  // ^^^^ new template parameter
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This change ([replacing virtual functions with CRTP](https://stackoverflow.com/a/16988933)) has helped some example programs to run twice as fast. See the example programs for illustration and details.

* New planner with an implementation of the S-star search algorithm for finding optimal path through simplicial complexes (https://arxiv.org/abs/1607.07009)

*NOTE:* Discrete Optimal Search Library (DOSL) is a fork of the Yet Another Graph-Search Based Planning Library (YAGSBPL)       hosted at https://github.com/subh83/YAGSBPL . YAGSBPL is now deprecated.

*******************************************************************************

Installation and Compilation of Examples:
----------------------------------------

**Installation** (_optional_)**:** 
DOSL is (to a large extent) template-based. There is nothing to build for the library itself.
Simply include the header "dosl/dosl" in your C++ code to use DOSL.
You can (optionally) install the headers in the system folder by running
```
    sudo make install
```

**Compilation:**
Quick compilation of all the examples in the `examples-dosl` folder (this will ask you to choose an algorithm):
```
    cd examples-dosl
    make all
```
Alternatively, run  `make simple`  or  `make advanced` to compile a selection of the example programs only (see `examples-dosl/makefile` for other make rules).
All executables are created in the `examples-dosl/bin` folder. 


**Running the Examples:**
After compilation, to run an example program:
```
    ./bin/<program_name>
```
Or to run all compiled example programs under the `./bin` folder:
```
    make run
```

******************************************************************************************

Basic Usage:
------------

DOSL uses graphs as the most basic form of discrete representation. In order to describe a graph to DOSL, the user needs to define two classes:

1. A class defining the type of node/vertex in the graph. This class needs to be derived from DOSL's `[AlgName]::Node` class (where`[AlgName]` is the name of the algorithm being used). The user must overload `operator==` for this class so that DOSL knows how to tell whether two nodes are the same or different.

2. A class defining the other aspects of the graph structure and the search problem. This class needs to be derived from DOSL's `[AlgName]::Algorithm` class. Typical members of this class (virtual members of `[AlgName]::Algorithm`) that the user needs to define are `getSuccessors`, `getStartNodes` and `stopSearch`.


Below is a bare-bones example, with explanations in comments, illustrating the use of DOSL with the AStar algorithm:


```C++
// standard headers
#include <stdio.h>
#include <vector>
#include <iostream>

// DOSL headers:
#include <dosl/dosl>

// ==============================================================================
/* The following class defines the type of a vertex of the graph.
       Needs to be derived from DOSL-provided class template 'AStar::Node<node_type,cost_type>' */

class myNode : public AStar::Node<myNode,double>
{
public:
    int x, y; // (x,y) coordinates defining a point on plane.
    
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
       Needs to be derived from DOSL-provided class template 'AStar::Algorithm<search_problem_class,node_type,cost_type>' */

class searchProblem : public AStar::Algorithm<searchProblem,myNode,double>
{
public:
    // user-defined problem parameters:
    myNode goal_node;
    // Constructors, if any
    searchProblem () { goal_node = myNode(150,100); }
    
    // -----------------------------------------------------------
    /* The following functions are use by the base class 'AStar::Algorithm' to determine 
           graph structure and search parameters. */
    
    /* Prototype for 'AStar::Algorithm<>::getSuccessors':
           void getSuccessors 
               (NodeType& n, std::vector<NodeType>* const s, std::vector<CostType>* const s);
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
    
    /* Prototype for 'AStar::Algorithm<>::getStartNodes':
           std::vector<NodeType> getStartNodes (void);
       Description: Should return the list of vertices(s) to start the search with. */
    
    std::vector<myNode> getStartNodes (void) {
        std::vector<myNode> startNodes;
        for (int a=0; a<1; ++a) {
            myNode tn (0, 0); // start node
            startNodes.push_back (tn);
        }
        return (startNodes);
    }
    
    /* Prototype for 'AStar::Algorithm<>::stopSearch':
           bool stopSearch (NodeType& n);
       Description: Determines whether to stop the search when 'n' is being expanded.
       Optional -- If not provided, search will terminate only when heap is empty. */
    
    bool stopSearch (myNode &n) {
        return (n==goal_node);
    }
};

// ==============================================================================

int main(int argc, char *argv[])
{
    searchProblem test_search_problem; // declare an instance of the search problem.
    test_search_problem.search(); // execute search.
    
    // get path from start to goal vertex.
    std::vector<myNode*> path = test_search_problem.reconstruct_pointer_path (test_search_problem.goal_node);
    
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

Encapsulations:
---------------

For common types of graph, encaptulations can be used to hide some of the details involved in defining the graph type. Currently there are two encaptsulations that are provided with DOSL:

__Single shortest path planning in an OpenCV image:__

Include `dosl/encapsulations/cvPathPlanner.tcc` in your code. Then you can compute the optimal path from a `start` pixel to a `goal` pixel in a OpenCV image `obs_map` by creating a variable of type `cvPathPlanner<[AlgName]>` (where `[AlgName]` is the name of the search algorithm) using the constructor
```C++
cvPathPlanner::cvPathPlanner (cv::Mat obs_map, ...);
```
and then calling the member
```C++
void cvPathPlanner::find_path (cv::Point start, cv::Point goal, bool vis=false);
```
Set `vis` to `true` to visualize the search.
The shortest path is stored in the member `std::vector<cv::Point> path`. You can directly draw the path in a matrix using either of the following functions:
```C++
void cvPathPlanner::draw_path (cv::Mat& in_map, cv::Mat& out_map, ...);
cv::Mat cvPathPlanner::draw_path (...);
```
See the example in `examples-dosl/src/simple/map2d_encapsulated_PathPlanning.cpp` for more details.

__Multiple shortest paths in different topological classes in an OpenCV image:__

Include `dosl/encapsulations/cvMulticlassPathPlanner.tcc` in your code. Then create an instance of `cvMulticlassPathPlanner<[AlgName]>` (where `[AlgName]` is the name of the search algorithm) using the constructor
```C++
cvMulticlassPathPlanner::cvMulticlassPathPlanner (cv::Mat obs_map, ...)
```
and compute the paths using member
```C++
void cvMulticlassPathPlanner::find_paths (cv::Point start, cv::Point goal, int nPaths=1, bool vis=false, int obsSizeThresh=0);
```
`obsSizeThresh` is the minimum size of obstacle that would create multiple classes of paths.
The shortest path is stored in the member `std::vector< std::vector< cv::Point > > paths`. You can directly draw the path in a matrix using either of the following functions:
```C++
void cvMulticlassPathPlanner::draw_paths (cv::Mat& in_map, cv::Mat& out_map, ...);
cv::Mat cvMulticlassPathPlanner::draw_paths (...)
```
See the example in `examples-dosl/src/simple/homotopy2d_encapsulated_PathPlanning.cpp` for more details.
******************************************************************************************

Documentation:
--------------

DOSL wiki is under construction: https://github.com/subh83/DOSL/wiki


*******************************************************************************

Citation:
--------

If you found this library useful in your research, please cite it in your paper as follows:

Suggested citation format:
_Subhrajit Bhattacharya, "Discrete Optimal Search Library (DOSL): A template-based C++ library for discrete optimal search", 2017. Available at https://github.com/subh83/DOSL ._

Bibtex entry:
```
 @misc{dosl,
     title = {Discrete Optimal Search Library (DOSL): A template-based C++ library for discrete optimal search},
     author = {Subhrajit Bhattacharya},
     note = {Available at https://github.com/subh83/DOSL },
     url = { https://github.com/subh83/DOSL },
     year = {2017}
 }
```

If you, in particular, use the S-star search algorithm, please use the following citation:

_Subhrajit Bhattacharya, "Towards optimal path computation in a simplicial complex", The International Journal of Robotics Research (IJRR), online-first, June, 2019. DOI: 10.1177/0278364919855422._

Bibtex entry:
```
@ARTICLE { Simplicial:star:19,
    AUTHOR = { Subhrajit Bhattacharya },
    TITLE = { Towards optimal path computation in a simplicial complex },
    JOURNAL = { The International Journal of Robotics Research (IJRR) },
    MONTH = { June },
    YEAR = { 2019 },
    NOTE = { DOI: 10.1177/0278364919855422 },
    URL = { https://doi.org/10.1177/0278364919855422 }
}
```

*******************************************************************************

Version history:
---------------

* May 2019: version 3.31: 'multiple definition' bug fix.

* May 2019: version 3.3 released: virtual functions replaced by CRTP in `[AlgName]::Algorithm` classes, thus making some programs run twice as fast; replaced node and simplex containers/heaps with pointers, making way for multi-thread programming in future; renaming of some variables and functions (back compatibility included); bug fixes.

* Nov 2018: version 3.26 released: Some bug fixes. Changes made to encapsulation member functions.

* Sep 2017: version 3.25 released: Added ThetaStar and SStar search algorithms; Organized components of algorithm under nested classes; More extensive and organized examples grouped into "simple" and "advanced"; Simple encapsulations for path planning in OpenCV matrices added.

* May 2017: version 3.1 released

* Nov 2016: version 3.0a released

* Discrete Optimal Search Library (DOSL) is a fork of
  the Yet Another Graph-Search Based Planning Library (YAGSBPL)
  hosted at https://github.com/subh83/YAGSBPL .
  YAGSBPL is now deprecated.



*******************************************************************************

License:
--------

```
/** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal Search Library (DOSL)                                              *
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
```
