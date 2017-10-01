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

#ifndef __DOSL_S_STAR_TCC
#define __DOSL_S_STAR_TCC
// user-readable
#define DOSL_ALGORITHM_SStar

// includes

#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <limits>
#include <type_traits>

#include "../utils/macros_constants.tcc"
#include "../utils/metric_simplex.tcc"
#include "../utils/stl_utils.tcc"
#include "_planner_bits.hpp"

// ====================================================================

#ifndef _S_STAR_DOSL_CASCADED_BACKTRACE  // everything seems to work just fine even without this! why???
#define _S_STAR_DOSL_CASCADED_BACKTRACE false
#endif

// 2: proctively fix directionality; 1: report when error due to directionality; 0: no check
#ifndef _DOSL_CHECK_GRAPH_DIRECTIONALITY
#define _DOSL_CHECK_GRAPH_DIRECTIONALITY 1
#endif

// --------------------------------------------------------------------
class SStar {
public:
    declare_alg_name("SStar"); // macro from '_planner_bits'
    
    // To be derived by user node type
    template <class nodeType, class costType=double> // CRTP. // costType needs to be a floating point type
    class Node : public MetricSimplexVertex <nodeType*, costType> // MetricSimplexVertex provides 'G' & 'Successors'
    {
    public:
        using MetricSimplexVertex<nodeType*,costType>::G;
        using MetricSimplexVertex<nodeType*,costType>::Expanded;
        using MetricSimplexVertex<nodeType*,costType>::Successors;

        // Basic vertex related types
        typedef costType  CostType;
        typedef nodeType  NodeType;
        
        // Simplex related types
        typedef costType DoubleType;
        typedef AMetricSimplex<NodeType*,DoubleType>  MetricSimplexType;
        typedef MetricSimplexCollection <NodeType*,DoubleType>  MetricSimplexContainerType;
        
        // Utility for looking up algorithm name
        declare_alg_name("SStar"); // macro from '_planner_bits'
        
        // For keeping track of hash insertion (new node creation)
        bool PostHashInsertInitiated;
        
        // Node specific variables
        int nodeHeapPos;
        int backTrackVertex;
        bool SuccessorsCreated;
        // _DOSL_SMALL_MAP<NodeType*,CostType> Successors; // inherited from 'MetricSimplexVertex'
        
        MetricSimplexType* CameFromSimplex;
        
        // -------------------------------------
        // constructors
        Node() : PostHashInsertInitiated(false), nodeHeapPos(-1),
                      backTrackVertex(0), SuccessorsCreated(false), CameFromSimplex(NULL) { }
        
        // TODO: Copy constructor to prevent copying these members
        
        void soft_reset (void) { // reset everything other than successor list
            PostHashInsertInitiated = false; nodeHeapPos = -1; Expanded = false;
            backTrackVertex = 0; SuccessorsCreated = false; CameFromSimplex = NULL;
        }
        
        // Define virtual functions of derived class
        inline CostType HeapKey() { return (G+getHeuristics()); }
        
        // ----------------------------------------------------------------------
        // Functions to be overwritten by user node type
        //      (need not be virtual since use is of only derived class members).
        // Need to have virtual members with same name in the problem class.
        inline int getHashBin (void) { return (0); }
        inline void getSuccessors (std::vector<nodeType>* s, std::vector<CostType>* c)
            { _dosl_default_fun_warn("SStar::[Algorithm|Node]::getSuccessors"); }
        inline CostType getHeuristics (void) { return ((CostType)0.0); }
        inline bool bookmarkNode (void) { return (false); }
        inline bool stopSearch (void) { return (false); }
        #if _DOSL_EVENTHANDLER
        void nodeEvent (unsigned int e) { }
        #endif
        inline void print (std::string head="", std::string tail="") { 
            _dosl_cout << _GREEN << head << " (" << this << ")" << GREEN_ << _dosl_endl;
            _dosl_cout << "Successors: " << Successors << _dosl_endl;
        }
        // Derived functions. Can also be directly overwritten.
        inline CostType getHeapKey (double subopEps) { return (G + (CostType)(subopEps*getHeuristics())); }
        
        // ----------------------------------------------------------------------
        // For path reconstruction using getPointerPathToNode. These operators need to be replaced in the derived class.
        virtual NodeType operator+(const NodeType &b) const { // n1 + n2 == n1
            _dosl_default_fun_warn ("SStar::Node::operator+");
            return (*((NodeType*)this));
        }
        virtual NodeType operator*(const DoubleType &c) const { // n1 * c  == n1 (right scalar multiplication)
            _dosl_default_fun_warn ("SStar::Node::operator*");
            return (*((NodeType*)this));
        } 
    };

    // ---------------------------------------------

    template <class nodeType, class costType=double>
    class Algorithm
    {
    // 'nodeType' should be derived fom 'SStar::Node<CostType>'
    public:
        typedef costType  CostType;
        typedef nodeType  NodeType;
        NodeType dummy_node; // forces compiler to generate code for the possibly template class 'nodeType'
        
        typedef costType DoubleType;
        typedef AMetricSimplex<NodeType*,/*MetricSimplexPointersUnorderedSetType,*/DoubleType>  MetricSimplexType;
        typedef MetricSimplexCollection <NodeType*,DoubleType>  MetricSimplexContainerType;
        typedef _DOSL_SMALL_MAP <NodeType*,DoubleType>  SimplexPointType;
        
        // Utility for looking up algorithm name
        declare_alg_name("SStar"); // macro from '_planner_bits'
        
        // Constants
        #if _DOSL_EVENTHANDLER
        enum nodeEventType {
            // Expanded: bit 0
                EXPANDED = 1,
            // Heap event: bits 1, 2
                HEAP = 2+4,
                PUSHED = 2, UPDATED = 4, POPPED = 2+4,
            // Unexpanding events: bit 3
                UNEXPANDED = 8,
            // errors
                ERROR = 16 // error if bit 4 is on. bit 0-3 will then contain information about the error.
        };
        #endif
        
        // parameters for problem
        double SubopEps;
        int ExpandCount;
        // verbose only
        int ProgressShowInterval;
        float timediff;
        clock_t  StartClock;
        time_t StartSecond;
        // Member variables
        std::vector<NodeType> StartNodes;
        std::vector<NodeType*> BookmarkNodePointers;
        
        // Node Set (Hash table)
        _DOSL_LARGE_UNORDERED_SET <NodeType, Algorithm, Algorithm>  AllNodesSet;
        bool operator()(const NodeType& n1, const NodeType& n2) { return (n1==n2); } // equal_to // TODO: NodeType const &
        // Heaps
        _DOSL_HEAP <NodeType*, Algorithm, Algorithm> NodeHeap;
        bool operator()(NodeType* const & np1, NodeType* const & np2) { return (getHeapKey(*np1) < getHeapKey(*np2)); } // less_than
        // store pointers to simplies, just for keeping track
        MetricSimplexContainerType  AllSimplices;
        
        // Constructors and initiators
        Algorithm ()
        { 
            // default parameters
            SubopEps = 1.0;
            
            // node hash table
            AllNodesSet.set_equal_to_function (this);
            AllNodesSet.set_hash_function (this, &Algorithm::getHashBin); // will call this->getHashBin(const NodeType& n)
            // node heap
            NodeHeap.set_compare_function (this);
            NodeHeap.set_heappos_function (this, &Algorithm::getNodeHeapPos);
            
            ProgressShowInterval = 10000;
        }
        
        // Main search functions
        void search (void);
        void clear (bool clearHeap=true, bool resetNodesInHash=true, unsigned int clearHash=0u);
        
        // Functions for reading results
        std::vector < std::unordered_map <NodeType*, DoubleType> >  reconstructPath (NodeType n);
        std::vector<nodeType*> reconstructPointerPath (nodeType n, bool convertWeightsToInt=true); // for compatibility with AStar
        // access other node data
        inline NodeType* getNodePointer (NodeType n) { return (AllNodesSet.get(n)); }
        inline CostType getCostsToNodes (NodeType n) { return (AllNodesSet.get(n)->G); }
        // bookmark nodes
        inline std::vector<NodeType*> getBookmarkNodePointers (void) { return (BookmarkNodePointers); }
        
        
        // -----------------------------------------------------
        // functions to be overwritten by user problem instance.
        // Also in node class
        inline virtual unsigned int getHashBin (NodeType& n) { return (n.getHashBin()); }
        inline virtual void getSuccessors (NodeType &n, std::vector<NodeType>* const s, std::vector<CostType>* const c)
                                            { return (n.getSuccessors(s,c)); }
        inline virtual CostType getHeuristics (NodeType& n) { return (n.getHeuristics()); }
        inline virtual std::vector<NodeType> getStartNodes (void)
            { _dosl_default_fun_warn("SStar::Algorithm::getStartNodes"); return (std::vector<NodeType>()); }
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
        // S-star specific:
        void CheckAndPushNodeIntoHeap (NodeType* nodeInHash_p, bool heapKeyUpdated=true, bool cascadeToChildren=false);
        
        // Temporary variables
        std::vector<NodeType> thisSuccessors;
        std::vector<CostType> thisTransitionCosts;
        int a;
    };

};

// =====================================================================================

// Basic member functions (not for end-user use)

template <class nodeType, class costType>
void SStar::Algorithm<nodeType,costType>::GenerateSuccessors (NodeType* nodeInHash_p)
{
    if ( !(nodeInHash_p->SuccessorsCreated) ) // Successors were not generated previously
    {
        thisSuccessors.clear();
        thisTransitionCosts.clear();
        getSuccessors (*nodeInHash_p, &thisSuccessors, &thisTransitionCosts);
        
        #if _DOSL_DEBUG
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
void SStar::Algorithm<nodeType,costType>::CheckAndPushNodeIntoHeap (NodeType* nodeInHash_p, bool heapKeyUpdated,
                                                                            bool cascadeToChildren)
{
    nodeInHash_p->Expanded = false;
    
    if (nodeInHash_p->nodeHeapPos < 0) { // may already in heap due to earlier initiated backtracking
        NodeHeap.insert (nodeInHash_p);
        #if _DOSL_EVENTHANDLER
        nodeEvent (*nodeInHash_p, PUSHED);
        #endif
    } else if (heapKeyUpdated) {
        NodeHeap.update (nodeInHash_p);
        #if _DOSL_EVENTHANDLER
        nodeEvent (*nodeInHash_p, UPDATED);
        #endif
    }
    
    // TODO: Implement for  cascadeToChildren = true
    if (cascadeToChildren)
        for (auto it=nodeInHash_p->ChildrenInfluenced.begin(); it!=nodeInHash_p->ChildrenInfluenced.end(); ++it)
            CheckAndPushNodeIntoHeap (*it, heapKeyUpdated, cascadeToChildren);
}

// -------------------------------------------------------------------------------------

template <class nodeType, class costType>
void SStar::Algorithm<nodeType,costType>::search (void)
{
    _dosl_verbose_head(1);
    
    StartNodes = getStartNodes();
    
    #if _DOSL_DEBUG
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
    CostType thisTransitionCost, test_g_val;
    
    // -----------------------------------
    // Insert start nodes into open list
    for (a=0; a<StartNodes.size(); a++) {
        // put start node in hash
        thisNodeInHash_p = AllNodesSet.get (StartNodes[a]);
        // set member variables for start node
        thisNodeInHash_p->Expanded = false;
        thisNodeInHash_p->backTrackVertex = 0;
        thisNodeInHash_p->G = (CostType)0.0;
        // push in heap
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
            if(ProgressShowInterval>0) {
                if (ExpandCount % ProgressShowInterval == 0) {
                    if (timediff>=0.0)  timediff = ((float)(clock()-StartClock)) / ((float)CLOCKS_PER_SEC);
                    _dosl_printf("Number of states Expanded: %d. Number of simplices ceated: %d. Node heap size: %d." 
                                        "Time elapsed: %F s.", ExpandCount, AllSimplices.size(), NodeHeap.size(), 
                                        ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
                }
            }
        ExpandCount++;
        
        // get the node with least F-value
        thisNodeInHash_p = NodeHeap.pop();
        
        // Generate the neighbours if they are already not generated
        GenerateSuccessors (thisNodeInHash_p);
        
        // Generate successors of unexpanded neighbors
        for (auto it=thisNodeInHash_p->Successors.begin(); it!=thisNodeInHash_p->Successors.end(); ++it) {
            thisNeighbourNodeInHash_p = it->first;
            if (!(thisNeighbourNodeInHash_p->Expanded)) 
                GenerateSuccessors (thisNeighbourNodeInHash_p);
        }
        
        #if _DOSL_CHECK_GRAPH_DIRECTIONALITY >= 2
        for (auto it=thisNodeInHash_p->Successors.begin(); it!=thisNodeInHash_p->Successors.end(); ++it)
            if ( it->first->SuccessorsCreated  &&  
                    it->first->Successors.find(thisNodeInHash_p)==it->first->Successors.end() ) {
                thisNodeInHash_p->print("node being expanded: ");
                it->first->print("neighbor not containg expanding node in its successor: ");
                _dosl_warn("Asymmetric successors (not an undirected graph as is required by SStar). Autofixing.\n");
                it->first->Successors [thisNodeInHash_p] = thisNodeInHash_p->Successors[it->first];
            }
        for (auto it=thisNeighbourNodeInHash_p->Successors.begin(); it!=thisNeighbourNodeInHash_p->Successors.end(); ++it)
            if ( it->first->SuccessorsCreated  &&  
                    it->first->Successors.find(thisNeighbourNodeInHash_p)==it->first->Successors.end() ) {
                thisNeighbourNodeInHash_p->print("node being updated: ");
                it->first->print("neighbor not containg expanding node in its successor: ");
                _dosl_warn("Asymmetric successors (not an undirected graph as is required by SStar). Autofixing.\n");
                it->first->Successors [thisNeighbourNodeInHash_p] = thisNeighbourNodeInHash_p->Successors[it->first];
            }
        #endif
        
        // -----------------------------------------------------
        // Expand
        
        thisNodeInHash_p->Expanded = true; // Put in closed list
        
        if (_dosl_verbose_on(1)) {
            _dosl_printf (_YELLOW "=========================================================" YELLOW_);
            thisNodeInHash_p->print("Now expanding: ");
            _dosl_printf ("G-value at expanding node: %f.", thisNodeInHash_p->G);
        }
        
        #if _DOSL_EVENTHANDLER
        //eventExpanded (*thisNodeInHash_p);
        nodeEvent (*thisNodeInHash_p, EXPANDED|POPPED);
        #endif
        
        // Check if we need to bookmark the node being Expanded
        if ( bookmarkNode (*thisNodeInHash_p) )
        {
            BookmarkNodePointers.push_back (thisNodeInHash_p);
            if (_dosl_verbose_on(0)) {
                if (timediff>=0.0)  timediff = ((float)(clock()-StartClock)) / ((float)CLOCKS_PER_SEC);
                thisNodeInHash_p->print ("Bookmarked a node: ");
                _dosl_printf ("... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                        ExpandCount, NodeHeap.size(), ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
            }
        }
        // Check if we need to stop furthur expansion
        if ( stopSearch (*thisNodeInHash_p) ) {
            if (_dosl_verbose_on(0)) 
                if (ProgressShowInterval>0) {
                    if (timediff>=0.0)  timediff = ((float)(clock()-StartClock)) / ((float)CLOCKS_PER_SEC);
                    thisNodeInHash_p->print ("Stopping search ('stopSearch' returned true) at: ");
                    _dosl_printf ("... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                            ExpandCount, NodeHeap.size(), ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
                }
            return;
        }
        
        // -----------------------------------------------------
        
        // Initiate the neighbours (if required) and update their G & F values
        for (auto it = thisNodeInHash_p->Successors.begin(); it!=thisNodeInHash_p->Successors.end(); ++it)
        {
            thisNeighbourNodeInHash_p = it->first;
            thisTransitionCost = it->second;
            
            if (_dosl_verbose_on(1)) {
                DOSL_INDENT; printf (_YELLOW "--------------------------------------" YELLOW_);
                printf("\n");
            }
            
            // ----------------------------------------
            bool neighborJustCreated = false;
            
            // An uninitiated neighbour node - definitely G & F values not set either.
            if (!(thisNeighbourNodeInHash_p->PostHashInsertInitiated)) {
                thisNeighbourNodeInHash_p->Expanded = false;
                thisNeighbourNodeInHash_p->backTrackVertex = 0;
                thisNeighbourNodeInHash_p->PostHashInsertInitiated = true; // Always set this when other variables have already been set
                
                neighborJustCreated = true;
                if (_dosl_verbose_on(1)) {
                    thisNeighbourNodeInHash_p->print("Child generated for the first time: ");
                    _dosl_printf ("Transition cost to child: %f.", thisTransitionCost);
                }
                
            }
            
            bool isTestingForUnexpanding = false;
            if (thisNeighbourNodeInHash_p->Expanded) {
                isTestingForUnexpanding = true;
                if (_dosl_verbose_on(1)) {
                    thisNeighbourNodeInHash_p->print("Child ALREADY expanded (will check improvement): ");
                }
            }
            else {
                if (_dosl_verbose_on(1)) {
                    thisNeighbourNodeInHash_p->print("Child not expanded (will check improvement): ");
                }
            }
            
            // --------------------------------------------
            
            // temporary variables
            MetricSimplexType* new_simplex = NULL;
            MetricSimplexType* best_simplex = NULL;
            MetricSimplexType* one_simplex = NULL;
            MetricSimplexType* tmp_new_simplex = NULL;
            MetricSimplexType* tmptmp_new_simplex = NULL;
            MetricSimplexType* test_came_from_simplex = NULL;
            NodeType* thisCommonNeighbor_p = NULL;
            
            // ----------------------------------------------------
            // Attached simplices
            
            one_simplex  = AllSimplices.createNewOneSimplex (thisNeighbourNodeInHash_p, thisNodeInHash_p, thisTransitionCost);
            if (_dosl_verbose_on(1)) {
                one_simplex->print ("Edge (1-simplex) created:");
            }
            
            
            std::unordered_set<MetricSimplexType*> attachedMaximalSimplices = 
                            //AllSimplices.constructAllAttachedMaximalSimplices (one_simplex, &allCommonNeighbors);
                            AllSimplices.getAllAttachedMaximalSimplices (one_simplex,
                                                                            COMPUTE_ALL, true, false, // <-- defauts
                                                                            true // <-- use ony expanded vertices
                                                                                ); // also computed G-score
            
            if (_dosl_verbose_on(2)) {
                _dosl_printf (_YELLOW "Number of attached maximal simplices of %x (possibly including self) = %d" YELLOW_, 
                                                    one_simplex, attachedMaximalSimplices.size());
            }
            
            // choose the one with lowest g-val
            test_g_val = std::numeric_limits<CostType>::max();
            for (auto it2=attachedMaximalSimplices.begin(); it2!=attachedMaximalSimplices.end(); ++it2) {
                if (_dosl_verbose_on(2)) {
                    (*it2)->print ("Attached maximal simplex:");
                }
                if ((*it2)->G < test_g_val) {
                    test_came_from_simplex = *it2;
                    test_g_val = (*it2)->G; // G-score through the different simplices.
                }
            }
            
            
            // ----------------------------------------------------
            
            if (test_came_from_simplex) {
                if (isTestingForUnexpanding)
                    ++(test_came_from_simplex->backTrackSimplex);
            }
            #if _DOSL_CHECK_GRAPH_DIRECTIONALITY >= 1
            else { // test_came_from_simplex==NULL
                for (auto it=thisNodeInHash_p->Successors.begin(); it!=thisNodeInHash_p->Successors.end(); ++it)
                    if ( it->first->SuccessorsCreated  &&  
                            it->first->Successors.find(thisNodeInHash_p)==it->first->Successors.end() ) {
                        thisNodeInHash_p->print("node being expanded: ");
                        it->first->print("neighbor not containg expanding node in its successor: ");
                        /*#if _DOSL_AUTO_FIX_GRAPH_DIRECTIONALITY
                        _dosl_warn("Asymmetric successors (not an undirected graph as is required by SStar). Autofixing.\n");
                        it->first->Successors [thisNodeInHash_p] = thisNodeInHash_p->Successors[it->first];
                        #else*/
                        nodeEvent (*thisNodeInHash_p, ERROR);
                        nodeEvent (*(it->first), ERROR);
                        _dosl_err("Asymmetric successors (not an undirected graph as is required by SStar).\n");
                        //#endif
                    }
                for (auto it=thisNeighbourNodeInHash_p->Successors.begin(); it!=thisNeighbourNodeInHash_p->Successors.end(); ++it)
                    if ( it->first->SuccessorsCreated  &&  
                            it->first->Successors.find(thisNeighbourNodeInHash_p)==it->first->Successors.end() ) {
                        thisNeighbourNodeInHash_p->print("node being updated: ");
                        it->first->print("neighbor not containg expanding node in its successor: ");
                        /*#if _DOSL_AUTO_FIX_GRAPH_DIRECTIONALITY
                        _dosl_warn("Asymmetric successors (not an undirected graph as is required by SStar). Autofixing.\n");
                        it->first->Successors [thisNeighbourNodeInHash_p] = thisNeighbourNodeInHash_p->Successors[it->first];
                        #else */
                        nodeEvent (*thisNeighbourNodeInHash_p, ERROR);
                        nodeEvent (*(it->first), ERROR);
                        _dosl_err("Asymmetric successors (not an undirected graph as is required by SStar).\n");
                        // #endif
                    }
            }
            #endif
            
            // ---------------------------------------------------------
            if (_dosl_verbose_on(1)) {
                _dosl_printf_nobreak ("Test g-score of child (%x) = %f (through simplex %x). Current g-score = ", 
                                    thisNeighbourNodeInHash_p, test_g_val, test_came_from_simplex);
                (thisNeighbourNodeInHash_p->G==std::numeric_limits<double>::max())?
                                        printf("INF\n") : printf("%f\n", thisNeighbourNodeInHash_p->G);
            }
            
            // compare & update
            if (neighborJustCreated) { // it's a 1-simplex for sure
                thisNeighbourNodeInHash_p->G = test_g_val;
                thisNeighbourNodeInHash_p->CameFromSimplex = test_came_from_simplex;
                thisNeighbourNodeInHash_p->CameFromSimplex->SetChildInfluence();
                if (_dosl_verbose_on(1)) {
                    _dosl_printf (_YELLOW "Set G of child (%x) to %f" YELLOW_, thisNeighbourNodeInHash_p, thisNeighbourNodeInHash_p->G);
                }
                
                NodeHeap.push (thisNeighbourNodeInHash_p);
                #if _DOSL_EVENTHANDLER
                nodeEvent (*thisNeighbourNodeInHash_p, PUSHED);
                #endif
            }
            else if (thisNeighbourNodeInHash_p->G - test_g_val > _MS_DOUBLE_EPS) {
                double previous_g_val =  thisNeighbourNodeInHash_p->G;
                
                thisNeighbourNodeInHash_p->G = test_g_val;
                thisNeighbourNodeInHash_p->CameFromSimplex->UnsetChildInfluence(); // unset for previous came-from-simplex
                thisNeighbourNodeInHash_p->CameFromSimplex = test_came_from_simplex;
                thisNeighbourNodeInHash_p->CameFromSimplex->SetChildInfluence();
                if (_dosl_verbose_on(1)) {
                    _dosl_printf (_YELLOW "Updated: G of child (%x) to %f." YELLOW_, 
                                                    thisNeighbourNodeInHash_p, thisNeighbourNodeInHash_p->G);
                    test_came_from_simplex->print ("New came-from simplex of the child:");
                }
                
                if (isTestingForUnexpanding) { // Previously expanded. Need to backtack since g-score is updated
                    // Unexpand
                    CheckAndPushNodeIntoHeap (thisNeighbourNodeInHash_p, true, _S_STAR_DOSL_CASCADED_BACKTRACE);
                    ++(thisNeighbourNodeInHash_p->backTrackVertex);
                    
                    // stop expanding this current vertex if the backtrack has put it back into heap.
                    if (_S_STAR_DOSL_CASCADED_BACKTRACE  &&  thisNodeInHash_p->nodeHeapPos>=0)
                        break; // will break the for loop over successors
                    
                    if (_dosl_verbose_on(1)) {
                    _dosl_printf ("Unexpanding started at %x: g-value changed from %f to %f. ", 
                                    thisNeighbourNodeInHash_p, previous_g_val, thisNeighbourNodeInHash_p->G);
                    }
                    #if _DOSL_EVENTHANDLER
                    nodeEvent (*thisNeighbourNodeInHash_p, UNEXPANDED);
                    #endif
                }
                else {
                    // Since thisNeighbourGraphNode->F is changed, re-arrange it in heap
                    NodeHeap.update (thisNeighbourNodeInHash_p);
                    
                    #if _DOSL_EVENTHANDLER
                    nodeEvent (*thisNeighbourNodeInHash_p, UPDATED);
                    #endif
                }
            }
        }
    }
    if (_dosl_verbose_on(0)) 
        if (ProgressShowInterval>0 && NodeHeap.empty()) {
            _dosl_printf ("Stopping search!! Heap is empty... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                       ExpandCount, NodeHeap.size(), ((timediff>=0.0) ? timediff : difftime(time(NULL),StartSecond)) );
        }
}

// -------------------------------------------------------------------------------------

template <class nodeType, class costType>
void SStar::Algorithm<nodeType,costType>::clear (bool clearHeap, bool resetNodesInHash, unsigned int clearHash)
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
                AllNodesSet.HashTable[a][b]->Expanded = false;
                AllNodesSet.HashTable[a][b]->G = std::numeric_limits<CostType>::infinity();
                AllNodesSet.HashTable[a][b]->PostHashInsertInitiated = false;
        }
}

// -------------------------------------------------------------------------------------
// path reconstruction

int DOUBLE_TO_INT_FACTOR = 1000000;

/* template <class nodeType, class DoubleType, bool convertWeightToInt=true>
nodeType computeWeightedSum (std::unordered_map <nodeType*, DoubleType> node_weight);

template <class nodeType, class DoubleType>
nodeType computeWeightedSum<nodeType,DoubleType,true> (std::unordered_map <nodeType*, DoubleType> node_weight) { }

template <class nodeType, class DoubleType>
nodeType computeWeightedSum<nodeType,DoubleType,false> (std::unordered_map <nodeType*, DoubleType> node_weight) { } */

// -----

template <class nodeType, class costType>
    std::vector < std::unordered_map <nodeType*, costType> > // node-weight pairs
        SStar::Algorithm<nodeType,costType>::reconstructPath (nodeType n) // n should be a node that is in the S-star hash
{
    _dosl_verbose_head(1);
    
    // create a 0-simplex out of np
    NodeType* np = AllNodesSet.get(n);
    
    std::unordered_map <nodeType*, DoubleType>  thisPt;
    std::vector < std::unordered_map <nodeType*, DoubleType> >  ret;
    
    // isert the first point
    thisPt.insert (std::make_pair(np, 1.0));
    ret.push_back (thisPt);
    
    // next point
    typename MetricSimplexContainerType::PathPointType  path_pt =
            typename MetricSimplexContainerType::PathPointType 
                        (np->CameFromSimplex->w, np->CameFromSimplex->G_cameFromPoint, thisPt);
    
    if (_dosl_verbose_on(0))  n.print("Reconstruct path called");
    if (_dosl_verbose_on(1))  _dosl_printf("Node in hash: %x.", np);
    while (path_pt.p.size() > 0) {
        if (_dosl_verbose_on(1)) {
            path_pt.print("path point: ");
        }
        ret.push_back (path_pt.p);
        path_pt = AllSimplices.findCamefromPoint (path_pt);
    }
    
    return (ret);
}

template <class nodeType, class costType>
    std::vector<nodeType*>  // convex combined using overloaded operators + and * of 'nodeType'
        SStar::Algorithm<nodeType,costType>::reconstructPointerPath (nodeType n, bool convertWeightsToInt)
{
    std::vector < std::unordered_map <nodeType*, DoubleType> >  path = reconstructPath (n);
    std::vector<nodeType*>  ret;
    nodeType tn;
    
    // weighted combination
    for (a=0; a<path.size(); ++a) {
        auto it = path[a].begin();
        
        if (convertWeightsToInt)
            tn = *(it->first) * ((int)(it->second * DOUBLE_TO_INT_FACTOR));
        else
            tn = *(it->first) * it->second;
        ++it;
        
        while (it!=path[a].end()) {
            if (convertWeightsToInt)
                tn = tn  +  *(it->first) * ((int)(it->second * DOUBLE_TO_INT_FACTOR));
            else
                tn = tn  +  *(it->first) * it->second;
        
            // tn = tn + (*(it->first) * it->second);
            ++it;
        }
        
        if (convertWeightsToInt)
            tn = tn * (1.0/DOUBLE_TO_INT_FACTOR);
        ret.push_back (AllNodesSet.get(tn));
    }
    
    return (ret);
}

// =====================================================================================

#endif
