#include <stdio.h>
#include <stdlib.h>
#include <string>
// Open CV:
#include <opencv/cv.h>
#include <opencv/cvaux.h>
#include <opencv/highgui.h> 

//#include <armadillo>

// ==============================================================================
// DOSL library

// Optional parameters: 
// #define _DOSL_DEBUG 0
// #define _DOSL_AUTOCORRECT 1
// #define _DOSL_VERBOSE 1  // TODO: to be replaced competely by _DOSL_VERBOSE_ITEMS in future.
// #define _DOSL_EVENTHANDLER 1

#define _DOSL_ALGORITHM  AStar
#include "dosl.h"

class myNode : public DOSL_ALGORITHM(Node)<myNode,double>
{
public:
    int x, y;
    myNode () { }
    myNode (int xx, int yy) : x(xx), y(yy) { }
    
    bool operator==(const myNode& n) const { return (x==n.x && y==n.y); } // This must be defined for the node
    void print (std::string head="", std::string tail="") const
            { std::cout << head << "x=" << x << ", y=" << y << std::endl; }
    
    // Functions to be overwritten
    int GetHashBin (void) const { return (abs((x>>4) + (y<<3) + (x<<4) + (y>>3))); }
};

// ==============================================================================

class SearchProblem : public DOSL_ALGORITHM(Problem)<myNode,double>
{
public:
    // Constructors, if any
    //SearchProblem () { }
    
    // -----------------------------------------------------------
    
    void GetSuccessors (myNode &n, std::vector<myNode>* s, std::vector<double>* c) {
        // This function should account for obstacles and size of environment.
        myNode tn;
        for (int a=-1; a<=1; a++)
            for (int b=-1; b<=1; b++) {
                if (a==0 && b==0) continue;
                
                tn.x = n.x + a;
                tn.y = n.y + b;
                // TODO
                s->push_back(tn);
                c->push_back(sqrt((double)(a*a+b*b))); 
            }
    }
    
    // -----------------------------------------------------------
    
    std::vector<myNode> GetStartNodes (void) {
        std::vector<myNode> startNodes;
        for (int a=0; a<1; ++a) {
            myNode tn (-150, -100);
            startNodes.push_back (tn);
        }
        return (startNodes);
    }
    
    // -----------------------------------------------------------
    
    bool stopSearch (myNode &n) {
        return (n.x==150 && n.y==100);
    }
};

// ==============================================================================

int main(int argc, char *argv[])
{
    SearchProblem test_search_problem;
    test_search_problem.AllNodesSet.set_hash_table_size (8192);
    test_search_problem.Search();
    
    // Get path
    std::vector<myNode*> path = test_search_problem.GetPointerPathToNode (myNode(120,80));
    
    printf ("\nPath: ");
    for (int a=path.size()-1; a>=0; --a)
        printf ("(%d,%d), ", path[a]->x, path[a]->y);
    
    return (1);
}

