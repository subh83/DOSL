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

#ifndef __DOSL_ThetaStar_TCC
#define __DOSL_ThetaStar_TCC
// user-readable
#define DOSL_ALGORITHM_ThetaStar

// includes

#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <limits>
#include <type_traits>

#include "../utils/macros_constants.tcc"
#include "../utils/stl_utils.tcc"
#include "planner_bits.hpp"

// ====================================================================

class ThetaStar {
public:
    declare_alg_name("ThetaStar"); // macro from '_planner_bits'

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
    class Node
    {
    public:
        typedef costType  CostType;
        typedef nodeType  NodeType;
        
        // Utility for looking up algorithm name:
        declare_alg_name("ThetaStar"); // macro from '_planner_bits'
        
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
        NodeType* CameFrom; // *** Theta-star: same as 'grand-parent'
        
        // -------------------------------------
        // constructors
        Node() : PostHashInsertInitiated(false), nodeHeapPos(-1),
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
        inline int getHashBin (void) { return (0); }
        inline void getSuccessors (std::vector<nodeType>* s, std::vector<CostType>* c) 
            { _dosl_default_fun_warn("ThetaStar::[Algorithm|Node]::getSuccessors"); }
        inline CostType getHeuristics (void) { return ((CostType)0.0); }
        inline bool bookmarkNode (void) { return (false); }
        inline bool stopSearch (void) { return (false); }
        #if _DOSL_EVENTHANDLER
        void nodeEvent (unsigned int e) { }
        #endif
        inline void print (std::string head="", std::string tail="") {
            _dosl_printf(_GREEN "%s (%x) %s" GREEN_, head.c_str(), this, tail.c_str());
        }
        // Derived functions. Can also be directly overwritten.
        inline CostType getHeapKey (double subopEps) { return (G + (CostType)(subopEps*getHeuristics())); }
    };


    template <class nodeType, class costType=double>
    class Algorithm
    {
    public:
        typedef costType  CostType;
        typedef nodeType  NodeType;
        NodeType dummy_node; // forces compiler to generate code for the possibly template class 'nodeType'
        
        // Utility for looking up algorithm name:
        declare_alg_name("ThetaStar"); // macro from '_planner_bits'
        
        // Constants
        #if _DOSL_EVENTHANDLER
        enum nodeEventType {
            // Expanded: bit 0
                EXPANDED = 1,
            // Heap event: bits 1, 2
                HEAP = 2+4,
                PUSHED = 2, UPDATED = 4, POPPED = 2+4,
            // for compatibility with S-star
                UNEXPANDED = 0,
            // errors
                ERROR = 16 // error if bit 4 is on. bit 0-3 will then contain information about the error.
        };
        #endif
        
        // parameters for problem
        double SubopEps;
        int ExpandCount;
        // only for verbode:
        int ProgressShowInterval;
        float timediff;
        clock_t  StartClock;
        time_t StartSecond;
        // Member variables
        std::vector<NodeType> StartNodes;
        std::vector<NodeType*> BookmarkNodePointers;
        
        // Heaps and Hashes
        _DOSL_LARGE_UNORDERED_SET <NodeType, Algorithm, Algorithm>  AllNodesSet;
        bool operator()(const NodeType& n1, const NodeType& n2) { return (n1==n2); } // equal_to
        // Heaps
        _DOSL_HEAP <NodeType*, Algorithm, Algorithm> NodeHeap;
        bool operator()(NodeType* const & np1, NodeType* const & np2) 
                { return (getHeapKey(*np1) < getHeapKey(*np2)); } // less_than
        
        // Constructors and initiators
        Algorithm () {
            // default parameters
            SubopEps = 1.0;
            
            // node hash table
            AllNodesSet.set_equal_to_function (this); // will call this->operator()(const NodeType& n1, const NodeType& n2)
            AllNodesSet.set_hash_function (this, &Algorithm::getHashBin); // will call this->getHashBin(const NodeType& n)
            // node heap
            NodeHeap.set_compare_function (this);
            NodeHeap.set_heappos_function (this, &Algorithm::getNodeHeapPos);
            
            ProgressShowInterval = 10000;
        }
        
        // Main search functions
        void search (void);
        void clear (bool clearHeap=true, bool resetNodesInHash=true, unsigned int clearHash=0u);
        
        // Functions for reading paths to arbitrary nodes
        std::vector < std::unordered_map <nodeType*, costType> > reconstructPath (nodeType n); // for compatibility with SStar
        std::vector<nodeType*> reconstructPointerPath (nodeType n);
        // access other node data
        inline NodeType* getNodePointer (NodeType n) { return (AllNodesSet.get(n)); }
        inline CostType getCostsToNodes (NodeType n) { return (AllNodesSet.get(n)->G); }
        // bookmark nodes
        inline std::vector<NodeType*> getBookmarkNodePointers (void) { return (BookmarkNodePointers); }
        
        
        // -----------------------------------------------------
        // functions to be overwritten by user problem instance.
        // Also in node class
        inline virtual unsigned int getHashBin (NodeType &n) { return (n.getHashBin()); }
        inline virtual void getSuccessors (NodeType &n, std::vector<NodeType>* const s, std::vector<CostType>* const c) 
                                            { return (n.getSuccessors(s,c)); }
        inline virtual bool isSegmentFree (NodeType &n1, NodeType &n2, CostType* c) // *** Theta-star specific
            { _dosl_default_fun_warn("ThetaStar::Algorithm::isSegmentFree"); *c=(CostType)0.0; return (true); } 
        inline virtual CostType getHeuristics (NodeType& n) { return (n.getHeuristics()); }
        inline virtual std::vector<NodeType> getStartNodes (void)
            { _dosl_default_fun_warn("ThetaStar::Algorithm::getStartNodes"); return (std::vector<NodeType>()); }
        inline virtual bool bookmarkNode (NodeType &n) { return (n.bookmarkNode()); }
        inline virtual bool stopSearch (NodeType &n) { return (n.stopSearch()); }
        #if _DOSL_EVENTHANDLER
        inline virtual void nodeEvent (NodeType &n, unsigned int e) { n.nodeEvent(e); }
        #endif
        inline virtual void print (NodeType &n, std::string head) { n.print(head); }
        // Derived functions. Can also be directly overwritten.
        inline virtual CostType getHeapKey (NodeType& n) { return (n.G + (CostType)(SubopEps*getHeuristics(n))); }
        
        // -------------------------------------------------
    private:
        // Derived and helper functions
        int& getNodeHeapPos (NodeType* const & np) { return (np->nodeHeapPos); }
        void GenerateSuccessors (NodeType* nodeInHash_p);
        
        // Temporary variables
        std::vector<NodeType> thisSuccessors;
        std::vector<CostType> thisTransitionCosts;
        int a;
    };

};

// =====================================================================================

template <class nodeType, class costType>
void ThetaStar::Algorithm<nodeType,costType>::GenerateSuccessors (NodeType* nodeInHash_p)
{
    if ( !(nodeInHash_p->SuccessorsCreated) ) // Successors were not generated previously
    {
        thisSuccessors.clear();
        thisTransitionCosts.clear();
        getSuccessors (*nodeInHash_p, &thisSuccessors, &thisTransitionCosts);
        
        #if _DOSL_DEBUG > 0
        if (thisSuccessors.size()!=thisTransitionCosts.size())
            _dosl_err("Number of successors (%d) is different from numer of transition costs (%d) as returned by 'getSuccessors'.", thisSuccessors.size(), thisTransitionCosts.size());
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
void ThetaStar::Algorithm<nodeType,costType>::search (void)
{
    _dosl_verbose_head(1);
    
    StartNodes = getStartNodes();
    
    #if _DOSL_DEBUG > 0
    if (StartNodes.size()==0)
        _dosl_err("No start node given! Please define the 'getStartNodes' function in the problem class.");
    #endif
    
    ExpandCount = 0;
    if (_dosl_verbose_on(0)) {
        timediff = 0.0;
        StartClock = clock();
        StartSecond = time(NULL);
    }
    
    // Temporary variables
    NodeType* thisNodeInHash_p;
    NodeType* thisNeighbourNodeInHash_p;
    NodeType* test_cameFrom;
    CostType thisTransitionCost, grandparentSegmentCost, test_g_val, temp_g_val;
    
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
        thisNodeInHash_p->F = getHeapKey(*thisNodeInHash_p); // TODO: this is not being used!
        if (thisNodeInHash_p->nodeHeapPos < 0) { // should not push the same node twice into heap
            NodeHeap.insert (thisNodeInHash_p);
            #if _DOSL_EVENTHANDLER
            nodeEvent (*thisNodeInHash_p, PUSHED);
            #endif
        }
        thisNodeInHash_p->PostHashInsertInitiated = true;
    }
    
    // -----------------------------------
    // Main loop
    while (!(NodeHeap.empty()))
    {
        if (_dosl_verbose_on(0))
            if (ProgressShowInterval>0) {
                if (ExpandCount % ProgressShowInterval == 0) {
                    if (timediff>=0.0)  timediff = ((float)(clock()-StartClock)) / ((float)CLOCKS_PER_SEC);
                    _dosl_printf("Number of states Expanded: %d. Heap size: %d. Time elapsed: %F s.", 
                            ExpandCount, NodeHeap.size(), ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
                }
            }
        ExpandCount++;
        
        // get the node with least F-value
        thisNodeInHash_p = NodeHeap.pop();
        if (_dosl_verbose_on(1)) {
            thisNodeInHash_p->print("Now expanding: ");
            _dosl_printf("G-value at expanding node: %f. F-value at expanding node: %f.", thisNodeInHash_p->G, thisNodeInHash_p->F);
        }
        
        // Generate the neighbours if they are already not generated
        GenerateSuccessors (thisNodeInHash_p);
        if (_dosl_verbose_on(1)) {
            _dosl_printf("Number of successors = %d.", thisNodeInHash_p->Successors.size());
        }
        
        // -----------------------------------------------------
        // Expand
        
        thisNodeInHash_p->Expanded = true; // Put in closed list
        
        #if _DOSL_EVENTHANDLER
        nodeEvent (*thisNodeInHash_p, EXPANDED|POPPED);
        #endif
        
        // Check if we need to bookmark the node being Expanded
        if ( bookmarkNode (*thisNodeInHash_p) )
        {
            BookmarkNodePointers.push_back (thisNodeInHash_p);
            if (_dosl_verbose_on(0)) {
                if (timediff>=0.0)  timediff = ((float)(clock()-StartClock)) / ((float)CLOCKS_PER_SEC);
                thisNodeInHash_p->print ("Bookmarked a node: ");
                _dosl_printf("... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                        ExpandCount, NodeHeap.size(), ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
            }
        }
        // Check if we need to stop furthur expansion
        if ( stopSearch (*thisNodeInHash_p) ) {
            if (_dosl_verbose_on(0))
                if (ProgressShowInterval>0) {
                    if (timediff>=0.0)  timediff = ((float)(clock()-StartClock)) / ((float)CLOCKS_PER_SEC);
                    thisNodeInHash_p->print ("Stopping search ('stopSearch' returned true) at: ");
                    _dosl_printf("... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                            ExpandCount, NodeHeap.size(), ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
                }
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
                thisNeighbourNodeInHash_p->F = getHeapKey (*thisNeighbourNodeInHash_p);
                thisNeighbourNodeInHash_p->Expanded = false;
                // Put in open list and continue to next neighbour
                thisNeighbourNodeInHash_p->PostHashInsertInitiated = true;
                                        // ^^^ Always set this when other variables have already been set
                
                if (_dosl_verbose_on(2)) {
                    thisNeighbourNodeInHash_p->print("Child (generated for the first time): ");
                    _dosl_printf("Transition cost to child: %f.", thisTransitionCost);
                }
                
                // push
                NodeHeap.push (thisNeighbourNodeInHash_p);
                #if _DOSL_EVENTHANDLER
                nodeEvent (*thisNeighbourNodeInHash_p, PUSHED);
                #endif
                
                continue;
            }
            
            if (_dosl_verbose_on(2)) {
                thisNeighbourNodeInHash_p->print("Child (generated previously): ");
                _dosl_printf("Transition cost to child: %f.", thisTransitionCost);
            }
            
            // Neighbour that is in closed list is to be skipped
            if (thisNeighbourNodeInHash_p->Expanded)
                continue;
            
            // Update CameFrom, G and F values if better
            // path 1
            test_cameFrom = thisNodeInHash_p;
            test_g_val = thisNodeInHash_p->G + thisTransitionCost;
            // path 2: *** Theta-star specific
            if (isSegmentFree (*(thisNodeInHash_p->CameFrom), *thisNeighbourNodeInHash_p, &grandparentSegmentCost)) { // path 2
                //printf("Theta-star: grandparentSegmentCost = %f\n", grandparentSegmentCost);
                temp_g_val = thisNodeInHash_p->CameFrom->G + grandparentSegmentCost;
                if (temp_g_val < test_g_val) {
                    test_cameFrom = thisNodeInHash_p->CameFrom;
                    test_g_val = temp_g_val;
                }
            }
            
            if (test_g_val < thisNeighbourNodeInHash_p->G) {
                thisNeighbourNodeInHash_p->G = test_g_val;
                thisNeighbourNodeInHash_p->F = getHeapKey (*thisNeighbourNodeInHash_p);
                thisNeighbourNodeInHash_p->CameFrom = test_cameFrom;
                thisNeighbourNodeInHash_p->lineageData = test_cameFrom->lineageData.next_generation();
                
                // Since thisNeighbourGraphNode->F is changed, re-arrange it in heap
                NodeHeap.update (thisNeighbourNodeInHash_p);
                #if _DOSL_EVENTHANDLER
                nodeEvent (*thisNeighbourNodeInHash_p, UPDATED);
                #endif
            }
            
            
            
        }
    }
    
    if (_dosl_verbose_on(0)) {
        if (ProgressShowInterval>0 && NodeHeap.empty()) {
            _dosl_printf("Stopping search!! Heap is empty... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                       ExpandCount, NodeHeap.size(), ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
        }
    }
}

// -------------------------------------------------------------------------------------

template <class nodeType, class costType>
void ThetaStar::Algorithm<nodeType,costType>::clear (bool clearHeap, bool resetNodesInHash, unsigned int clearHash)
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
std::vector<nodeType*> ThetaStar::Algorithm<nodeType,costType>::reconstructPointerPath (nodeType n)
{
    _dosl_verbose_head(1);
        
    std::vector<NodeType*> thisPath;
    // Reconstruct path
    NodeType* thisNodeInHash_p = AllNodesSet.get(n);
    if (_dosl_verbose_on(0))  n.print("Reconstruct path called on ");
    if (_dosl_verbose_on(1))  _dosl_printf("Node in hash: %x.", thisNodeInHash_p);
    while (thisNodeInHash_p) {
        if (_dosl_verbose_on(1)) thisNodeInHash_p->print();
        thisPath.push_back (thisNodeInHash_p);
        thisNodeInHash_p = thisNodeInHash_p->CameFrom;
    }
    return (thisPath);
}

// ===========================
// For compatibity with s-star

template <class nodeType, class costType>
    std::vector < std::unordered_map <nodeType*, costType> >
        ThetaStar::Algorithm<nodeType,costType>::reconstructPath (nodeType n)
{
    std::vector<nodeType*> npVec = reconstructPointerPath (n);
    
    std::vector < std::unordered_map <nodeType*, costType> >  ret;
    for (auto it=npVec.begin(); it!=npVec.end(); ++it) {
        std::unordered_map <nodeType*, costType> thisPt;
        thisPt[(*it)] = 1.0;
        ret.push_back (thisPt);
    }
    return (ret);
}

#endif
