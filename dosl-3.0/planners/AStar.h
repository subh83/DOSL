/** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal Search Library (DOSL)                                              *
*    A template-based C++ library for discrete (graph) search                            *
*    Version 3.0                                                                         *
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
*    Contact:  subhrajit@gmail.com ,  http://subhrajit.net/                              *
*                                                                                        *
*                                                                                        *
*************************************************************************************** **/

#ifndef __DOSL_AStar_H
#define __DOSL_AStar_H

// includes

#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <limits>
#include <type_traits>

#include "utils/macros_constants.h"

// ====================================================================

#if _DOSL_EVENTHANDLER
enum AStarNodeEventType {
    // Expanded: bit 0
        EXPANDED = 1,
    // Heap event: bits 1, 2
        HEAP = 2+4,
        PUSHED = 2, UPDATED = 4, POPPED = 2+4,
    // for compatibility with S-star
        UNEXPANDED = 0
};
#endif

// --------------------------------------------------------------------

class LineageDataType
{
public:
    bool defined;
    int id, generation;
    
    LineageDataType () : defined(false) { }
    LineageDataType (int i, int g=0) : defined(true), id(i), generation(g) { }
    
    bool is_set (void) { return (defined); }
    LineageDataType next_generation(void) { LineageDataType ret(*this); ++(ret.generation); return (ret); }
};


// To be derived by user node type
template <class nodeType, class costType=double> // CRTP
class AStarNode // : public HeapItem<costType>, public HashItem
{
public:
    typedef costType  CostType;
    typedef nodeType  NodeType;
    
    // For keeping track of hash insertion (new node creation)
    bool PostHashInsertInitiated;
    
    // Node specific variables
    CostType F, G;
    bool Expanded; // Whether in closed list or not
    int nodeHeapPos;
    LineageDataType lineageData; // stores which start node the node came from
    
    // successors
    _DOSL_SMALL_MAP<NodeType*,CostType> Successors;
    bool SuccessorsCreated;
    
    // for reconstruction
    NodeType* CameFrom;
    
    // -------------------------------------
    // constructors
    AStarNode() : PostHashInsertInitiated(false), nodeHeapPos(-1),
                  Expanded(false), SuccessorsCreated(false), CameFrom(NULL), lineageData(LineageDataType()) { }
    
    void soft_reset (void) { // reset everything other than successor list
        PostHashInsertInitiated = false; nodeHeapPos = -1; Expanded = false;
        CameFrom = NULL; lineageData = LineageDataType();
    }
    
    // TODO: Copy constructor to prevent copying these members
    
    // Define virtual functions of derived class
    inline CostType HeapKey() { return (F); }
    
    // ----------------------------------------------------------------------
    // Functions to be overwritten by user node type
    //      (need not be virtual since use is of only derived class members).
    // Need to have virtual members with same name in the problem class.
    inline int GetHashBin (void) { return (0); }
    inline void GetSuccessors (std::vector<nodeType>* s, std::vector<CostType>* c) { }
    inline CostType GetHeuristics (void) { return ((CostType)0.0); }
    inline bool bookmarkNode (void) { return (false); }
    inline bool stopSearch (void) { return (false); }
    #if _DOSL_EVENTHANDLER
    /* inline void eventUpdated (void) { }
    inline void eventExpanded (void) { } */
    void Event (unsigned int e) { }
    #endif
    #if _DOSL_VERBOSE
    inline virtual void print (std::string head="", int nTab=0, std::string tail="\n") { 
        std::cout << head << " (";
        printf("%x", this);
        std::cout << ")" << tail;
    }
    #endif
    // Derived functions. Can also be directly overwritten.
    inline CostType GetHeapKey (double subopEps) { return (G + (CostType)(subopEps*GetHeuristics())); }
};


template <class nodeType, class costType=double>
class AStarProblem
{
public:
    typedef costType  CostType;
    typedef nodeType  NodeType;
    
    // parameters for problem
    double SubopEps;
    int ExpandCount;
    #if _DOSL_VERBOSE
    int ProgressShowInterval;
    clock_t  StartClock;
    time_t StartSecond;
    #endif
    // Member variables
    std::vector<NodeType> StartNodes;
    std::vector<NodeType*> BookmarkNodePointers;
    
    // Heaps and Hashes
    _DOSL_LARGE_UNORDERED_SET <NodeType, AStarProblem, AStarProblem>  AllNodesSet;
    bool operator()(const NodeType& n1, const NodeType& n2) { return (n1==n2); } // equal_to
    // Heaps
    _DOSL_HEAP <NodeType*, AStarProblem, AStarProblem> NodeHeap;
    bool operator()(NodeType* const & np1, NodeType* const & np2) 
            { return (GetHeapKey(*np1) < GetHeapKey(*np2)); } // less_than
    
    // Constructors and initiators
    AStarProblem () {
        // default parameters
        SubopEps = 1.0;
        
        // node hash table
        AllNodesSet.set_equal_to_function (this); // will call this->operator()(const NodeType& n1, const NodeType& n2)
        AllNodesSet.set_hash_function (this, &AStarProblem::GetHashBin); // will call this->GetHashBin(const NodeType& n)
        // node heap
        NodeHeap.set_compare_function (this);
        NodeHeap.set_heappos_function (this, &AStarProblem::GetNodeHeapPos);
        
        #if _DOSL_VERBOSE
        ProgressShowInterval = 10000;
        #endif
    }
    
    // Main search functions
    void Search (void);
    void Clear (bool clearHeap=true, bool resetNodesInHash=true, unsigned int clearHash=0u);
    
    // Functions for reading paths to arbitrary nodes
    inline NodeType* GetNodePointer (NodeType n) { return (AllNodesSet.get(n)); }
    inline CostType GetCostsToNodes (NodeType n) { return (AllNodesSet.get(n)->G); }
    std::vector<nodeType*> GetPointerPathToNode (nodeType const & n);
    // bookmark nodes
    inline std::vector<NodeType*> GetBookmarkNodePointers (void) { return (BookmarkNodePointers); }
    std::vector<CostType> GetCostsToBookmarkNodes (void);
    std::vector< std::vector<NodeType*> > GetPointerPathsToBookmarkNodes (void);
    // for compatibility with s-star:
    std::vector < std::unordered_map <nodeType*, costType> > ReconstructPath (nodeType const & n);
    
    // -----------------------------------------------------
    // functions to be overwritten by user problem instance.
    // Also in node class
    inline virtual unsigned int GetHashBin (const NodeType &n) { return (n.GetHashBin()); }
    inline virtual void GetSuccessors (NodeType &n, std::vector<NodeType>* const s, std::vector<CostType>* const c) 
            { return (n.GetSuccessors(s,c)); }
    inline virtual CostType GetHeuristics (NodeType& n) { return (n.GetHeuristics()); }
    inline virtual bool bookmarkNode (NodeType &n) { return (n.bookmarkNode()); }
    inline virtual bool stopSearch (NodeType &n) { return (n.stopSearch()); }
    #if _DOSL_EVENTHANDLER
    inline virtual void NodeEvent (NodeType &n, unsigned int e) { n.Event(e); }
    #endif
    #if _DOSL_VERBOSE
    inline virtual void print (NodeType &n, std::string head) { n.print(head); }
    #endif
    // Derived functions. Can also be directly overwritten.
    inline virtual CostType GetHeapKey (NodeType& n) {
        return (n.G + (CostType)(SubopEps*GetHeuristics(n)));
        // TODO: Make more efficient by reading 'F'. But be careful about updated G values.
    }
    // Specific to searcher class
    inline virtual std::vector<NodeType> GetStartNodes (void) { return (std::vector<NodeType>()); }
    
    // -------------------------------------------------
    // Derived and helper functions
    int& GetNodeHeapPos (NodeType* const & np) { return (np->nodeHeapPos); }
    
    void GenerateSuccessors (NodeType* nodeInHash_p);
    
    // Temporary variables
    std::vector<NodeType> thisSuccessors;
    std::vector<CostType> thisTransitionCosts;
    int a;
};

// =====================================================================================

template <class nodeType, class costType>
void AStarProblem<nodeType,costType>::GenerateSuccessors (NodeType* nodeInHash_p)
{
    if ( !(nodeInHash_p->SuccessorsCreated) ) // Successors were not generated previously
    {
        thisSuccessors.clear();
        thisTransitionCosts.clear();
        GetSuccessors (*nodeInHash_p, &thisSuccessors, &thisTransitionCosts);
        
        #if _DOSL_DEBUG
        if (thisSuccessors.size()!=thisTransitionCosts.size())
            _dosl_err("Number of successors (%d) is different from numer of transition costs (%d) as returned by 'GetSuccessors'.", thisSuccessors.size(), thisTransitionCosts.size());
        #endif
        
        nodeInHash_p->Successors.reserve (thisSuccessors.size()); // reserve space for fast pushing
        for (a=0; a<thisSuccessors.size(); a++) {
            thisSuccessors[a].soft_reset(); // in case the successors were created from copies of this
            nodeInHash_p->Successors.insert (
                _DOSL_SMALL_MAP_pairfun (AllNodesSet.get(thisSuccessors[a]), thisTransitionCosts[a]) );
        }
        nodeInHash_p->SuccessorsCreated = true;
    }
}

// -------------------------------------------------------------------------------------

template <class nodeType, class costType>
void AStarProblem<nodeType,costType>::Search (void)
{
    StartNodes = GetStartNodes();
    
    #if _DOSL_DEBUG
    if (StartNodes.size()==0)
        _dosl_err("No start node given! Please define the 'GetStartNodes' function in the problem class.");
    #endif
    
    ExpandCount = 0;
    #if _DOSL_VERBOSE
    float timediff = 0.0;
    StartClock = clock();
    StartSecond = time(NULL);
    #endif
    
    // Temporary variables
    NodeType* thisNodeInHash_p;
    NodeType* thisNeighbourNodeInHash_p;
    CostType thisTransitionCost, test_g_val;
    
    // -----------------------------------
    // Insert start nodes into open list
    for (int a=0; a<StartNodes.size(); a++) {
        // put start node in hash
        thisNodeInHash_p = AllNodesSet.get (StartNodes[a]);
        // set member variables for start node
        thisNodeInHash_p->CameFrom = NULL;
        if (!(thisNodeInHash_p->lineageData.is_set()))
            thisNodeInHash_p->lineageData = LineageDataType(a); // lineage a, generation 0
        thisNodeInHash_p->Expanded = false;
        thisNodeInHash_p->G = (CostType)0.0;
        thisNodeInHash_p->F = GetHeapKey(*thisNodeInHash_p); // TODO: this is not being used!
        if (thisNodeInHash_p->nodeHeapPos < 0) { // should not push the same node twice into heap
            NodeHeap.insert (thisNodeInHash_p);
            #if _DOSL_EVENTHANDLER
            NodeEvent (*thisNodeInHash_p, PUSHED);
            #endif
        }
        thisNodeInHash_p->PostHashInsertInitiated = true;
    }
    
    // -----------------------------------
    // Main loop
    while (!(NodeHeap.empty()))
    {
        #if _DOSL_VERBOSE
        if(ProgressShowInterval>0) {
            if (ExpandCount % ProgressShowInterval == 0) {
                if (timediff>=0.0)  timediff = ((float)(clock()-StartClock)) / ((float)CLOCKS_PER_SEC);
                printf("Number of states Expanded: %d. Heap size: %d. Time elapsed: %F s.\n", 
                        ExpandCount, NodeHeap.size(), ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
            }
        }
        #endif
        ExpandCount++;
        
        // Get the node with least F-value
        thisNodeInHash_p = NodeHeap.pop();
        #if _DOSL_VERBOSE > 1
        thisNodeInHash_p->print("Now expanding: ");
        printf("G-value at expanding node: %f. F-value at expanding node: %f. ", thisNodeInHash_p->G, thisNodeInHash_p->F);
        #endif
        
        // Generate the neighbours if they are already not generated
        GenerateSuccessors (thisNodeInHash_p);
        #if _DOSL_VERBOSE > 1
        printf ("Number of successors = %d\n", thisNodeInHash_p->Successors.size());
        #endif
        
        // -----------------------------------------------------
        // Expand
        
        thisNodeInHash_p->Expanded = true; // Put in closed list
        
        #if _DOSL_EVENTHANDLER
        NodeEvent (*thisNodeInHash_p, EXPANDED|POPPED);
        #endif
        
        // Check if we need to bookmark the node being Expanded
        if ( bookmarkNode (*thisNodeInHash_p) )
        {
            BookmarkNodePointers.push_back (thisNodeInHash_p);
            #if _DOSL_VERBOSE
            if (timediff>=0.0)  timediff = ((float)(clock()-StartClock)) / ((float)CLOCKS_PER_SEC);
            thisNodeInHash_p->print ("Bookmarked a node: ");
            printf("... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.\n", 
                    ExpandCount, NodeHeap.size(), ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
            #endif
        }
        // Check if we need to stop furthur expansion
        if ( stopSearch (*thisNodeInHash_p) ) {
            #if _DOSL_VERBOSE
            if (ProgressShowInterval>0) {
                if (timediff>=0.0)  timediff = ((float)(clock()-StartClock)) / ((float)CLOCKS_PER_SEC);
                thisNodeInHash_p->print ("Stopping search ('stopSearch' returned true) at: ");
                printf("... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.\n", 
                        ExpandCount, NodeHeap.size(), ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
            }
            #endif
            return;
        }
        
        // Initiate the neighbours (if required) and update their G & F values
        for (auto it = thisNodeInHash_p->Successors.begin(); it!=thisNodeInHash_p->Successors.end(); ++it)
        {
            thisNeighbourNodeInHash_p = it->first;
            thisTransitionCost = it->second;
            
            // An uninitiated neighbour node - definitely G & F values not set either.
            if (!(thisNeighbourNodeInHash_p->PostHashInsertInitiated)) {
                thisNeighbourNodeInHash_p->CameFrom = thisNodeInHash_p;
                thisNeighbourNodeInHash_p->lineageData = thisNodeInHash_p->lineageData.next_generation();
                thisNeighbourNodeInHash_p->G = thisNodeInHash_p->G + thisTransitionCost;
                thisNeighbourNodeInHash_p->F = GetHeapKey (*thisNeighbourNodeInHash_p);
                thisNeighbourNodeInHash_p->Expanded = false;
                // Put in open list and continue to next neighbour
                thisNeighbourNodeInHash_p->PostHashInsertInitiated = true;
                                        // ^^^ Always set this when other variables have already been set
                
                #if _DOSL_VERBOSE > 1
                thisNeighbourNodeInHash_p->print("\tChild generated for the first time: ");
                printf("\t\tTransition cost to child: %f \n", thisTransitionCost);
                #endif
                
                // push
                NodeHeap.push (thisNeighbourNodeInHash_p);
                #if _DOSL_EVENTHANDLER
                NodeEvent (*thisNeighbourNodeInHash_p, PUSHED);
                #endif
                
                continue;
            }
            
            // Neighbour that is in closed list is to be skipped
            if (thisNeighbourNodeInHash_p->Expanded)
                continue;
            
            // Update CameFrom, G and F values if better
            test_g_val = thisNodeInHash_p->G + thisTransitionCost;
            if (test_g_val < thisNeighbourNodeInHash_p->G) {
                thisNeighbourNodeInHash_p->G = test_g_val;
                thisNeighbourNodeInHash_p->F = GetHeapKey (*thisNeighbourNodeInHash_p);
                thisNeighbourNodeInHash_p->CameFrom = thisNodeInHash_p;
                thisNeighbourNodeInHash_p->lineageData = thisNodeInHash_p->lineageData.next_generation();
                
                // Since thisNeighbourGraphNode->F is changed, re-arrange it in heap
                NodeHeap.update (thisNeighbourNodeInHash_p);
                #if _DOSL_EVENTHANDLER
                NodeEvent (*thisNeighbourNodeInHash_p, UPDATED);
                #endif
            }
        }
    }
    #if _DOSL_VERBOSE
    if (ProgressShowInterval>0 && NodeHeap.empty())
        printf("Stopping search!! Heap is empty... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.\n", 
                   ExpandCount, NodeHeap.size(), ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
    #endif
}

// -------------------------------------------------------------------------------------

template <class nodeType, class costType>
void AStarProblem<nodeType,costType>::Clear (bool clearHeap, bool resetNodesInHash, unsigned int clearHash)
{
    NodeType* thisNodeInHash_p;
    if (clearHeap) NodeHeap.clear();
        /* while (!(NodeHeap.empty()))
            NodeHeap.pop(); */
    
    if (clearHash>0)
        AllNodesSet.clear ( ((clearHash==1)?false:true) );
    
    else if (resetNodesInHash)
        for (int a=0; a<AllNodesSet.HashTableSize; ++a)
            for (int b=0; b<AllNodesSet.HashTable[a].size(); ++b) {
                AllNodesSet.HashTable[a][b]->soft_reset();
        }
}

// -------------------------------------------------------------------------------------
// For reading results:

template <class nodeType, class costType>
std::vector< std::vector<nodeType*> > AStarProblem<nodeType,costType>::GetPointerPathsToBookmarkNodes(void)
{
    std::vector< std::vector<NodeType*> > paths;
    std::vector<NodeType*> thisPath;
    for (auto it = BookmarkNodePointers.begin(); it != BookmarkNodePointers.end(); ++it) {
        thisPath.clear();
        // Reconstruct path
        NodeType* thisNodeInHash_p = *it;
        while (thisNodeInHash_p) {
            thisPath.push_back (thisNodeInHash_p);
            thisNodeInHash_p = thisNodeInHash_p->CameFrom;
        }
        paths.push_back(thisPath);
    }
    return (paths);
}

template <class nodeType, class costType>
std::vector<costType> AStarProblem<nodeType,costType>::GetCostsToBookmarkNodes(void)
{
    std::vector<CostType> ret;
    for (auto it = BookmarkNodePointers.begin(); it != BookmarkNodePointers.end(); ++it)
        ret.push_back (it->G);
    return (ret);
}

template <class nodeType, class costType>
std::vector<nodeType*> AStarProblem<nodeType,costType>::GetPointerPathToNode (nodeType const & n)
{
    std::vector<NodeType*> thisPath;
    // Reconstruct path
    NodeType* thisNodeInHash_p = AllNodesSet.get(n);
    #if _DOSL_VERBOSE > 1
    n.print("Reconstruct path called");
    printf("Node in hash: %x\n", thisNodeInHash_p);
    #endif
    while (thisNodeInHash_p) {
        #if _DOSL_VERBOSE > 1
        thisNodeInHash_p->print();
        #endif
        thisPath.push_back (thisNodeInHash_p);
        thisNodeInHash_p = thisNodeInHash_p->CameFrom;
    }
    return (thisPath);
}

// ===========================
// For compatibity with s-star

template <class nodeType, class costType>
    std::vector < std::unordered_map <nodeType*, costType> >
        AStarProblem<nodeType,costType>::ReconstructPath (nodeType const & n)
{
    std::vector<nodeType*> npVec = GetPointerPathToNode (n);
    
    std::vector < std::unordered_map <nodeType*, costType> >  ret;
    for (auto it=npVec.begin(); it!=npVec.end(); ++it) {
        std::unordered_map <nodeType*, costType> thisPt;
        thisPt[(*it)] = 1.0;
        ret.push_back (thisPt);
    }
    return (ret);
}

#endif
