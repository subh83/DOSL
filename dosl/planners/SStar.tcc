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
#include <vector>
#include <limits>
#include <memory>
#include <type_traits>

#include "../utils/macros_constants.tcc"
#include "../utils/metric_simplex.tcc"
#include "../utils/stl_utils.tcc"
#include "planner_bits.hpp"

// ====================================================================

#ifndef _S_STAR_DOSL_CASCADED_BACKTRACE  // everything seems to work just fine even without this! why???
#define _S_STAR_DOSL_CASCADED_BACKTRACE false
#endif

#ifndef _S_STAR_DOSL_GENERATE_GRAND_CHILDREN
#define _S_STAR_DOSL_GENERATE_GRAND_CHILDREN true
#endif

// 2: proctively fix directionality; 1: report when error due to directionality; 0: no check
#ifndef _DOSL_CHECK_GRAPH_DIRECTIONALITY
#define _DOSL_CHECK_GRAPH_DIRECTIONALITY 2
#endif

// --------------------------------------------------------------------
class SStar {
public:
    declare_alg_name("SStar"); // macro from '_planner_bits'
    
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
    
    
    enum DataClearingMode {
        // for nodes:
        CLEAR_NODE_SUCCESSORS = 1, // bit 0
        CLEAR_NODE_LINEAGE = 2, // bit 1
        CLEAR_NODE_DATA = 3,
        // search data
        CLEAR_NODES_HEAP = 4, // bit 2
        CLEAR_NODES_HASHTABLE = 8 // bit 3
    };
    
    // To be derived by user node type
    template <class NodeType, class _CostType=double> // CRTP. // CostType needs to be a floating point type
    class Node : public BinaryHeapElement, public MetricSimplexVertex <NodeType*, _CostType> // MetricSimplexVertex provides 'g_score' & 'successors'
    {
    public:
        // Basic vertex related types
        typedef _CostType  CostType;
        
        using MetricSimplexVertex<NodeType*,CostType>::g_score;
        using MetricSimplexVertex<NodeType*,CostType>::expanded;
        using MetricSimplexVertex<NodeType*,CostType>::successors;
        
        
        // Simplex related types
        typedef CostType DoubleType;
        typedef AMetricSimplex<NodeType*,DoubleType>  MetricSimplexType;
        typedef MetricSimplexCollection <NodeType*,DoubleType>  MetricSimplexContainerType;
        
        // Utility for looking up algorithm name
        declare_alg_name("SStar"); // macro from '_planner_bits'
        
        // For keeping track of hash insertion (new node creation)
        bool post_hash_insert_initiated;
        
        LineageDataType lineage_data;
        
        // Node specific variables
        int backTrackVertex;
        bool successors_created;
        // _DOSL_SMALL_MAP<NodeType*,CostType> successors; // inherited from 'MetricSimplexVertex'
        
        MetricSimplexType* CameFromSimplex;
        
        // -------------------------------------
        // constructors
        Node() : post_hash_insert_initiated(false),
                      backTrackVertex(0), successors_created(false), CameFromSimplex(NULL), lineage_data(LineageDataType()) { }
        // TODO: Copy constructor to prevent copying these members
        
        // pseudo-destructor
        void clear_search_data (unsigned int mode = CLEAR_NODE_SUCCESSORS);
        
        // ----------------------------------------------------------------------
        // Functions to be overwritten by user node type
        inline int getHashBin (void) { return (0); }
        bool operator==(const NodeType& n)
            { _dosl_default_fun_warn("'SStar::Node::operator==' OR 'SStar::Algorithm::equalTo'"); }
        inline void getSuccessors (std::vector<NodeType>* s, std::vector<CostType>* c)
            { _dosl_default_fun_warn("SStar::[Algorithm|Node]::getSuccessors"); }
        inline CostType getHeuristics (void) { return ((CostType)0.0); }
        inline bool bookmarkNode (void) { return (false); }
        inline bool stopSearch (void) { return (false); }
        #if _DOSL_EVENTHANDLER
        void nodeEvent (unsigned int e) { }
        #endif
        inline void print (std::string head="", std::string tail="") { 
            _dosl_cout << _GREEN << head << " (" << this << ")" << GREEN_ << _dosl_endl;
            _dosl_cout << "successors: " << successors << _dosl_endl;
        }
        // Derived functions. Can also be directly overwritten.
        inline CostType getHeapKey (double subopEps) { return (g_score + (CostType)(subopEps*getHeuristics())); }
        
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

    template <class AlgDerived, class NodeType, class _CostType=double>
    class Algorithm
    {
    // 'NodeType' should be derived fom 'SStar::Node<CostType>'
    public:
        typedef _CostType  CostType;
        NodeType dummy_node; // forces compiler to generate code for the possibly template class 'NodeType'
        
        typedef CostType DoubleType;
        typedef AMetricSimplex<NodeType*,/*MetricSimplexPointersUnorderedSetType,*/DoubleType>  MetricSimplexType;
        typedef MetricSimplexCollection <NodeType*,DoubleType>  MetricSimplexContainerType;
        typedef _DOSL_SMALL_MAP <NodeType*,DoubleType>  SimplexPointType;
        
        // Utility for looking up algorithm name
        declare_alg_name("SStar"); // macro from '_planner_bits'
        
        // Constants
        #if _DOSL_EVENTHANDLER
        enum NodeEventType {
            // expanded: bit 0
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
        
        // user definable parameters for problem
        double subopt_eps;
        
        // counters
        int expand_count;
        // verbose only
        int progress_show_interval;
        ChronoTimer timer;
        // Member variables
        std::vector<NodeType> start_nodes;
        std::vector<NodeType*> bookmarked_node_pointers;
        
        // =====================================================
        // Main search data containers.
        
        // ---------------------------
        // Node Hash
        
        FUNCTOR_BEGIN (NodeHasherFunc, AlgDerived) // declares 'data_p' of type AlgDerived*
            unsigned int operator()(NodeType& n) { return data_p->getHashBin(n); }
        FUNCTOR_END (node_hasher_instance) // instanciates Node_hasher
        
        FUNCTOR_BEGIN (NodeEqualToFunc, AlgDerived) // declares 'data_p' of type AlgDerived*
            bool operator()(const NodeType& n1, const NodeType& n2) { return data_p->equalTo(n1,n2); }
        FUNCTOR_END (node_equal_to_instance) // instanciates Node_equal_to
        
        typedef _DOSL_LARGE_UNORDERED_SET<NodeType, NodeHasherFunc, NodeEqualToFunc>  NodesSetType;
        std::shared_ptr<NodesSetType>  all_nodes_set_p;
        
        // ---------------------------
        // Heaps
        
        FUNCTOR_BEGIN (NodeKeyLessThanFunc, AlgDerived) // declares 'data_p' of type AlgDerived*
            bool operator()(NodeType*& np1, NodeType*& np2) { return (data_p->getHeapKey(*np1) < data_p->getHeapKey(*np2)); }
        FUNCTOR_END (node_key_less_than_instance) // instanciates Node_hasher
        
        typedef _DOSL_HEAP <NodeType*, NodeKeyLessThanFunc>  NodeHeapType;
        std::shared_ptr<NodeHeapType>  node_heap_p;
        
        // ---------------------------
        // store pointers to simplies, just for keeping track
        MetricSimplexContainerType  AllSimplices;
        
        // =====================================================
        // ---------------------------
        // Constructors and initiators
        Algorithm (AlgDerived* shared_instance_p = NULL);
        
        // Main search functions
        void search (void);
        void clear (unsigned int mode = CLEAR_NODE_DATA | CLEAR_NODES_HEAP);
        
        // Functions for reading results
        std::vector < std::unordered_map <NodeType*, DoubleType> >  reconstruct_weighted_pointer_path (NodeType n);
        std::vector<NodeType*> reconstruct_pointer_path (NodeType n, bool convertWeightsToInt=true); // for compatibility with AStar
        // access other node data
        inline NodeType* get_node_pointer (NodeType n) { return (all_nodes_set_p->get(n)); }
        inline CostType get_costs_to_nodes (NodeType n) { return (all_nodes_set_p->get(n)->g_score); }
        // bookmark nodes
        inline std::vector<NodeType*> get_bookmark_node_pointers (void) { return (bookmarked_node_pointers); }
        
        
        // -----------------------------------------------------
        // functions to be overwritten by user problem instance.
        // Also in node class
        inline unsigned int getHashBin (NodeType& n) { return (n.getHashBin()); }
        bool equalTo (const NodeType &n1, const NodeType &n2) { return (n1==n2); }
        inline void getSuccessors (NodeType &n, std::vector<NodeType>* const s, std::vector<CostType>* const c)
                                            { return (n.getSuccessors(s,c)); }
        inline CostType getHeuristics (NodeType& n) { return (n.getHeuristics()); }
        inline std::vector<NodeType> getStartNodes (void)
            { _dosl_default_fun_warn("SStar::Algorithm::getStartNodes"); return (std::vector<NodeType>()); }
        inline bool bookmarkNode (NodeType &n) { return (n.bookmarkNode()); }
        inline bool stopSearch (NodeType &n) { return (n.stopSearch()); }
        #if _DOSL_EVENTHANDLER
        inline void nodeEvent (NodeType &n, unsigned int e) { n.nodeEvent(e); }
        #endif
        inline virtual void print (NodeType &n, std::string head) { n.print(head); }
        // Derived functions. Can also be directly overwritten.
        inline CostType getHeapKey (NodeType& n) { return (n.g_score + (CostType)(subopt_eps * _this->getHeuristics(n))); }
        
        // -------------------------------------------------
    private:
        AlgDerived* _this;
        
        // Derived and helper functions
        void generate_successors (NodeType* node_in_hash_p);
        // S-star specific:
        void CheckAndPushNodeIntoHeap (NodeType* node_in_hash_p, bool heapKeyUpdated=true, bool cascadeToChildren=false);
        
        // Temporary variables
        std::vector<NodeType> this_successors;
        std::vector<CostType> this_transition_costs;
        int a;
        
        // Preventing making copies of an instance of the search class
        Algorithm (const Algorithm& alg);
        Algorithm& operator=(const Algorithm& alg);
    };

};

// =====================================================================================

template <class AlgDerived, class NodeType, class CostType>
SStar::Algorithm<AlgDerived,NodeType,CostType>::Algorithm (AlgDerived* shared_instance_p) : _this (static_cast<AlgDerived*>(this))
{ 
    // Copy pointers to containers
    if (shared_instance_p)
        all_nodes_set_p = shared_instance_p->all_nodes_set_p;
    else {
        node_hasher_instance = NodeHasherFunc (_this);
        node_equal_to_instance = NodeEqualToFunc (_this);
        all_nodes_set_p = std::shared_ptr<NodesSetType>(new NodesSetType(4096, node_hasher_instance, node_equal_to_instance) );
    }
    
    if (shared_instance_p)
        node_heap_p = shared_instance_p->node_heap_p;
    else {
        node_key_less_than_instance = NodeKeyLessThanFunc (_this);
        node_heap_p = std::shared_ptr<NodeHeapType>(new NodeHeapType(node_key_less_than_instance) );
    }
    
    // default parameters
    subopt_eps = 1.0;
    progress_show_interval = 10000;
}

// Basic member functions (not for end-user use)

template <class AlgDerived, class NodeType, class CostType>
void SStar::Algorithm<AlgDerived,NodeType,CostType>::generate_successors (NodeType* node_in_hash_p)
{
    if ( !(node_in_hash_p->successors_created) ) // successors were not generated previously
    {
        this_successors.clear();
        this_transition_costs.clear();
        _this->getSuccessors (*node_in_hash_p, &this_successors, &this_transition_costs);
        
        #if _DOSL_DEBUG
        if (this_successors.size()!=this_transition_costs.size())
            _dosl_err("Number of successors (%d) is different from numer of transition costs (%d) as returned by 'getSuccessors'.", this_successors.size(), this_transition_costs.size());
        #endif
        
        node_in_hash_p->successors.reserve (this_successors.size()); // reserve space for fast pushing
        for (a=0; a<this_successors.size(); a++) {
            this_successors[a].clear_search_data (CLEAR_NODE_SUCCESSORS); // in case the successors were created from copies of this
            node_in_hash_p->successors.insert (
                _DOSL_SMALL_MAP_pairfun (all_nodes_set_p->get(this_successors[a]), this_transition_costs[a]) );
        }
        node_in_hash_p->successors_created = true;
    }
}

// -------------------------------------------------------------------------------------

template <class AlgDerived, class NodeType, class CostType>
void SStar::Algorithm<AlgDerived,NodeType,CostType>::CheckAndPushNodeIntoHeap (NodeType* node_in_hash_p, bool heapKeyUpdated,
                                                                            bool cascadeToChildren)
{
    node_in_hash_p->expanded = false;
    
    if (node_in_hash_p->heap_pos < 0) { // may already in heap due to earlier initiated backtracking
        node_heap_p->insert (node_in_hash_p);
        #if _DOSL_EVENTHANDLER
        _this->nodeEvent (*node_in_hash_p, PUSHED);
        #endif
    } else if (heapKeyUpdated) {
        node_heap_p->update (node_in_hash_p);
        #if _DOSL_EVENTHANDLER
        _this->nodeEvent (*node_in_hash_p, UPDATED);
        #endif
    }
    
    // TODO: Implement for  cascadeToChildren = true
    if (cascadeToChildren)
        for (auto it=node_in_hash_p->children_influenced.begin(); it!=node_in_hash_p->children_influenced.end(); ++it)
            CheckAndPushNodeIntoHeap (*it, heapKeyUpdated, cascadeToChildren);
}

// -------------------------------------------------------------------------------------

template <class AlgDerived, class NodeType, class CostType>
void SStar::Algorithm<AlgDerived,NodeType,CostType>::search (void)
{
    _dosl_verbose_head(1);
    
    start_nodes = _this->getStartNodes();
    
    #if _DOSL_DEBUG
    if (start_nodes.size()==0)
        _dosl_err("No start node given! Please define the 'getStartNodes' function in the problem class.");
    #endif
    
    expand_count = 0;
    if (_dosl_verbose_on(0)) {
        timer.start();
    }
    
    // Temporary variables
    NodeType* thisNodeInHash_p;
    NodeType* this_neighbour_node_in_hash_p;
    CostType thisTransitionCost, test_g_val;
    
    // -----------------------------------
    // Insert start nodes into open list
    for (a=0; a<start_nodes.size(); a++) {
        start_nodes[a].clear_search_data (CLEAR_NODE_SUCCESSORS); // in case this node is output of a previous search
        // put start node in hash
        thisNodeInHash_p = all_nodes_set_p->get (start_nodes[a]);
        // set member variables for start node
        if (!(thisNodeInHash_p->lineage_data.is_set() && thisNodeInHash_p->lineage_data.generation==0))
            thisNodeInHash_p->lineage_data = LineageDataType(a);
        thisNodeInHash_p->expanded = false;
        thisNodeInHash_p->backTrackVertex = 0;
        thisNodeInHash_p->g_score = (CostType)0.0;
        // push in heap
        if (thisNodeInHash_p->heap_pos < 0) { // should not push the same node twice into heap
            node_heap_p->insert (thisNodeInHash_p);
            #if _DOSL_EVENTHANDLER
            _this->nodeEvent (*thisNodeInHash_p, PUSHED);
            #endif
        }
        thisNodeInHash_p->post_hash_insert_initiated = true;
    }
    
    // -----------------------------------
    // Main loop
    while (!(node_heap_p->empty()))
    {
        if (_dosl_verbose_on(0))
            if(progress_show_interval>0) {
                if (expand_count % progress_show_interval == 0) {
                    _dosl_printf("Number of states expanded: %d. Number of simplices ceated: %d. Node heap size: %d." 
                                        "Time elapsed: %f_score s.", expand_count, AllSimplices.size(), node_heap_p->size(), timer.read());
                }
            }
        expand_count++;
        
        // get the node with least f_score-value
        thisNodeInHash_p = node_heap_p->pop();
        
        // Generate the neighbours if they are already not generated
        generate_successors (thisNodeInHash_p);
        
        // Generate successors of unexpanded neighbors
        #if _S_STAR_DOSL_GENERATE_GRAND_CHILDREN
        for (auto it=thisNodeInHash_p->successors.begin(); it!=thisNodeInHash_p->successors.end(); ++it) {
            this_neighbour_node_in_hash_p = it->first;
            if (!(this_neighbour_node_in_hash_p->expanded)) 
                generate_successors (this_neighbour_node_in_hash_p);
        }
        #endif
        
        #if _DOSL_CHECK_GRAPH_DIRECTIONALITY >= 2
        for (auto it=thisNodeInHash_p->successors.begin(); it!=thisNodeInHash_p->successors.end(); ++it) {
            if ( it->first->successors_created  &&  
                    it->first->successors.find(thisNodeInHash_p)==it->first->successors.end() ) {
                thisNodeInHash_p->print("node being expanded: ");
                it->first->print("neighbor not containg expanding node in its successor: ");
                _dosl_warn("Asymmetric successors (not an undirected graph as is required by SStar). Autofixing.\n");
                it->first->successors [thisNodeInHash_p] = thisNodeInHash_p->successors[it->first];
            }
            // fixing second generation (grand-children) connections:
            this_neighbour_node_in_hash_p = it->first;
            for (auto it2=this_neighbour_node_in_hash_p->successors.begin(); it2!=this_neighbour_node_in_hash_p->successors.end(); ++it2)
                if ( it2->first->successors_created  &&  
                        it2->first->successors.find(this_neighbour_node_in_hash_p)==it2->first->successors.end() ) {
                    this_neighbour_node_in_hash_p->print("node being updated: ");
                    it2->first->print("neighbor not containg expanding node in its successor: ");
                    _dosl_warn("Asymmetric successors (not an undirected graph as is required by SStar). Autofixing.\n");
                    it2->first->successors [this_neighbour_node_in_hash_p] = this_neighbour_node_in_hash_p->successors[it2->first];
                }
        }
        #endif
        
        // -----------------------------------------------------
        // Expand
        
        thisNodeInHash_p->expanded = true; // Put in closed list
        
        if (_dosl_verbose_on(1)) {
            _dosl_printf (_YELLOW "=========================================================" YELLOW_);
            thisNodeInHash_p->print("Now expanding: ");
            _dosl_printf ("g_score-value at expanding node: %f.", thisNodeInHash_p->g_score);
        }
        
        #if _DOSL_EVENTHANDLER
        //eventExpanded (*thisNodeInHash_p);
        _this->nodeEvent (*thisNodeInHash_p, EXPANDED|POPPED);
        #endif
        
        // Check if we need to bookmark the node being Expanded
        if ( _this->bookmarkNode (*thisNodeInHash_p) )
        {
            bookmarked_node_pointers.push_back (thisNodeInHash_p);
            if (_dosl_verbose_on(0)) {
                thisNodeInHash_p->print ("Bookmarked a node: ");
                _dosl_printf ("... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                        expand_count, node_heap_p->size(), timer.read());
            }
        }
        // Check if we need to stop furthur expansion
        if ( _this->stopSearch (*thisNodeInHash_p) ) {
            if (_dosl_verbose_on(0)) 
                if (progress_show_interval>0) {
                    thisNodeInHash_p->print ("Stopping search ('stopSearch' returned true) at: ");
                    _dosl_printf ("... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                            expand_count, node_heap_p->size(), timer.read());
                }
            return;
        }
        
        // -----------------------------------------------------
        
        // Initiate the neighbours (if required) and update their g_score & f_score values
        for (auto it = thisNodeInHash_p->successors.begin(); it!=thisNodeInHash_p->successors.end(); ++it)
        {
            this_neighbour_node_in_hash_p = it->first;
            thisTransitionCost = it->second;
            
            if (_dosl_verbose_on(1)) {
                DOSL_INDENT; printf (_YELLOW "--------------------------------------" YELLOW_);
                printf("\n");
            }
            
            // ----------------------------------------
            bool neighborJustCreated = false;
            
            // An uninitiated neighbour node - definitely g_score & f_score values not set either.
            if (!(this_neighbour_node_in_hash_p->post_hash_insert_initiated)) {
                this_neighbour_node_in_hash_p->expanded = false;
                this_neighbour_node_in_hash_p->lineage_data = thisNodeInHash_p->lineage_data.next_generation();
                this_neighbour_node_in_hash_p->backTrackVertex = 0;
                this_neighbour_node_in_hash_p->post_hash_insert_initiated = true; // Always set this when other variables have already been set
                
                neighborJustCreated = true;
                if (_dosl_verbose_on(1)) {
                    this_neighbour_node_in_hash_p->print("Child generated for the first time: ");
                    _dosl_printf ("Transition cost to child: %f.", thisTransitionCost);
                }
                
            }
            
            bool isTestingForUnexpanding = false;
            if (this_neighbour_node_in_hash_p->expanded) {
                isTestingForUnexpanding = true;
                if (_dosl_verbose_on(1)) {
                    this_neighbour_node_in_hash_p->print("Child ALREADY expanded (will check improvement): ");
                }
            }
            else {
                if (_dosl_verbose_on(1)) {
                    this_neighbour_node_in_hash_p->print("Child not expanded (will check improvement): ");
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
            
            one_simplex  = AllSimplices.create_new_one_simplex (this_neighbour_node_in_hash_p, thisNodeInHash_p, thisTransitionCost);
            if (_dosl_verbose_on(1)) {
                one_simplex->print ("Edge (1-simplex) created:");
            }
            
            
            std::unordered_set<MetricSimplexType*> attachedMaximalSimplices = 
                            //AllSimplices.constructAllAttachedMaximalSimplices (one_simplex, &allCommonNeighbors);
                            AllSimplices.get_all_attached_maximal_simplices (one_simplex,
                                                                            COMPUTE_ALL, // things_to_compute: default
                                                                            true, // return_base_if_maximal: default
                                                                            false, // force_compute: default
                                                                            true // use_only_expnded_nodes: use ony expanded vertices
                                                                                ); // <-- this also computed g_score
            
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
                if ((*it2)->g_score < test_g_val) {
                    test_came_from_simplex = *it2;
                    test_g_val = (*it2)->g_score; // g_score-score through the different simplices.
                }
            }
            
            
            // ----------------------------------------------------
            
            if (test_came_from_simplex) {
                if (isTestingForUnexpanding)
                    ++(test_came_from_simplex->back_track_simplex);
            }
            #if _DOSL_CHECK_GRAPH_DIRECTIONALITY >= 1
            else { // test_came_from_simplex==NULL
                for (auto it=thisNodeInHash_p->successors.begin(); it!=thisNodeInHash_p->successors.end(); ++it) {
                    if ( it->first->successors_created  &&  
                            it->first->successors.find(thisNodeInHash_p)==it->first->successors.end() ) {
                        thisNodeInHash_p->print("node being expanded: ");
                        it->first->print("neighbor not containg expanding node in its successor: ");
                        /*#if _DOSL_AUTO_FIX_GRAPH_DIRECTIONALITY
                        _dosl_warn("Asymmetric successors (not an undirected graph as is required by SStar). Autofixing.\n");
                        it->first->successors [thisNodeInHash_p] = thisNodeInHash_p->successors[it->first];
                        #else*/
                        _this->nodeEvent (*thisNodeInHash_p, ERROR);
                        _this->nodeEvent (*(it->first), ERROR);
                        _dosl_err("Asymmetric successors (not an undirected graph as is required by SStar).\n");
                        //#endif
                    }
                    // checking second generation (grand-children) connections:
                    this_neighbour_node_in_hash_p = it->first;
                    for (auto it2=this_neighbour_node_in_hash_p->successors.begin(); it2!=this_neighbour_node_in_hash_p->successors.end(); ++it2)
                        if ( it2->first->successors_created  &&  
                                it2->first->successors.find(this_neighbour_node_in_hash_p)==it2->first->successors.end() ) {
                            this_neighbour_node_in_hash_p->print("node being updated: ");
                            it2->first->print("neighbor not containg expanding node in its successor: ");
                            /*#if _DOSL_AUTO_FIX_GRAPH_DIRECTIONALITY
                            _dosl_warn("Asymmetric successors (not an undirected graph as is required by SStar). Autofixing.\n");
                            it2->first->successors [this_neighbour_node_in_hash_p] = this_neighbour_node_in_hash_p->successors[it2->first];
                            #else */
                            _this->nodeEvent (*this_neighbour_node_in_hash_p, ERROR);
                            _this->nodeEvent (*(it2->first), ERROR);
                            _dosl_err("Asymmetric successors (not an undirected graph as is required by SStar).\n");
                            // #endif
                        }
                }
            }
            #endif
            
            // ---------------------------------------------------------
            if (_dosl_verbose_on(1)) {
                _dosl_printf_nobreak ("Test g-score of child (%x) = %f (through simplex %x). Current g-score = ", 
                                    this_neighbour_node_in_hash_p, test_g_val, test_came_from_simplex);
                (this_neighbour_node_in_hash_p->g_score==std::numeric_limits<double>::max())?
                                        printf("INF\n") : printf("%f\n", this_neighbour_node_in_hash_p->g_score);
            }
            
            // compare & update
            if (neighborJustCreated) { // it's a 1-simplex for sure
                this_neighbour_node_in_hash_p->g_score = test_g_val;
                this_neighbour_node_in_hash_p->lineage_data = test_came_from_simplex->p[0]->lineage_data.next_generation();
                this_neighbour_node_in_hash_p->CameFromSimplex = test_came_from_simplex;
                this_neighbour_node_in_hash_p->CameFromSimplex->SetChildInfluence();
                if (_dosl_verbose_on(1)) {
                    _dosl_printf (_YELLOW "Set g_score of child (%x) to %f" YELLOW_, this_neighbour_node_in_hash_p, this_neighbour_node_in_hash_p->g_score);
                }
                
                node_heap_p->push (this_neighbour_node_in_hash_p);
                #if _DOSL_EVENTHANDLER
                _this->nodeEvent (*this_neighbour_node_in_hash_p, PUSHED);
                #endif
            }
            else if (this_neighbour_node_in_hash_p->g_score - test_g_val > _MS_DOUBLE_EPS) {
                double previous_g_val =  this_neighbour_node_in_hash_p->g_score;
                
                this_neighbour_node_in_hash_p->g_score = test_g_val;
                this_neighbour_node_in_hash_p->lineage_data = test_came_from_simplex->p[0]->lineage_data.next_generation();
                this_neighbour_node_in_hash_p->CameFromSimplex->UnsetChildInfluence(); // unset for previous came-from-simplex
                this_neighbour_node_in_hash_p->CameFromSimplex = test_came_from_simplex;
                this_neighbour_node_in_hash_p->CameFromSimplex->SetChildInfluence();
                if (_dosl_verbose_on(1)) {
                    _dosl_printf (_YELLOW "Updated: g_score of child (%x) to %f." YELLOW_, 
                                                    this_neighbour_node_in_hash_p, this_neighbour_node_in_hash_p->g_score);
                    test_came_from_simplex->print ("New came-from simplex of the child:");
                }
                
                if (isTestingForUnexpanding) { // Previously expanded. Need to backtack since g-score is updated
                    // Unexpand
                    CheckAndPushNodeIntoHeap (this_neighbour_node_in_hash_p, true, _S_STAR_DOSL_CASCADED_BACKTRACE);
                    ++(this_neighbour_node_in_hash_p->backTrackVertex);
                    
                    // stop expanding this current vertex if the backtrack has put it back into heap.
                    if (_S_STAR_DOSL_CASCADED_BACKTRACE  &&  thisNodeInHash_p->heap_pos>=0)
                        break; // will break the for loop over successors
                    
                    if (_dosl_verbose_on(1)) {
                    _dosl_printf ("Unexpanding started at %x: g-value changed from %f to %f. ", 
                                    this_neighbour_node_in_hash_p, previous_g_val, this_neighbour_node_in_hash_p->g_score);
                    }
                    #if _DOSL_EVENTHANDLER
                    _this->nodeEvent (*this_neighbour_node_in_hash_p, UNEXPANDED);
                    #endif
                }
                else {
                    // Since thisNeighbourGraphNode->f_score is changed, re-arrange it in heap
                    node_heap_p->update (this_neighbour_node_in_hash_p);
                    
                    #if _DOSL_EVENTHANDLER
                    _this->nodeEvent (*this_neighbour_node_in_hash_p, UPDATED);
                    #endif
                }
            }
        }
    }
    if (_dosl_verbose_on(0)) 
        if (progress_show_interval>0 && node_heap_p->empty()) {
            _dosl_printf ("Stopping search!! Heap is empty... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                       expand_count, node_heap_p->size(), timer.read());
        }
}

// -------------------------------------------------------------------------------------

template <class NodeType, class CostType>
void SStar::Node<NodeType,CostType>::clear_search_data (unsigned int mode) { // reset everything other than successor list
    post_hash_insert_initiated = false; expanded = false; //came_from = NULL;
    heap_pos = -1; g_score = (CostType)0.0; //f_score = (CostType)0.0;
    if (mode & CLEAR_NODE_SUCCESSORS) { successors.clear(); successors_created = false; }
    if (mode & CLEAR_NODE_LINEAGE) lineage_data = LineageDataType();
    backTrackVertex = 0; CameFromSimplex = NULL;
}

// -------------------------------------------------------------------------------------

template <class AlgDerived, class NodeType, class CostType>
void SStar::Algorithm<AlgDerived,NodeType,CostType>::clear (unsigned int mode)
{
    NodeType* thisNodeInHash_p;
    if (mode & CLEAR_NODES_HEAP) node_heap_p->clear();
        /* while (!(node_heap_p->empty()))
            node_heap_p->pop(); */
    
    if (mode & CLEAR_NODES_HASHTABLE)
        all_nodes_set_p->clear ();
    
    else if (mode & CLEAR_NODE_DATA)
        for (auto fit=all_nodes_set_p->begin(); fit != all_nodes_set_p->end(); ++fit) 
            fit->clear_search_data (mode);
    
    /*else if (resetNodesInHash)
        for (int a=0; a<all_nodes_set_p->hash_table_size; ++a)
            for (int b=0; b<all_nodes_set_p->HashTable[a].size(); ++b) {
                all_nodes_set_p->HashTable[a][b]->expanded = false;
                all_nodes_set_p->HashTable[a][b]->g_score = std::numeric_limits<CostType>::infinity();
                all_nodes_set_p->HashTable[a][b]->post_hash_insert_initiated = false;
        }*/
}

// -------------------------------------------------------------------------------------
// path reconstruction

#define DOUBLE_TO_INT_FACTOR 1000000

/* template <class NodeType, class DoubleType, bool convertWeightToInt=true>
NodeType computeWeightedSum (std::unordered_map <NodeType*, DoubleType> node_weight);

template <class NodeType, class DoubleType>
NodeType computeWeightedSum<NodeType,DoubleType,true> (std::unordered_map <NodeType*, DoubleType> node_weight) { }

template <class NodeType, class DoubleType>
NodeType computeWeightedSum<NodeType,DoubleType,false> (std::unordered_map <NodeType*, DoubleType> node_weight) { } */

// -----

template <class AlgDerived, class NodeType, class CostType>
    std::vector < std::unordered_map <NodeType*, CostType> > // node-weight pairs
        SStar::Algorithm<AlgDerived,NodeType,CostType>::reconstruct_weighted_pointer_path (NodeType n) // n should be a node that is in the S-star hash
{
    _dosl_verbose_head(1);
    
    // create a 0-simplex out of np
    NodeType* np = all_nodes_set_p->get(n);
    
    std::unordered_map <NodeType*, DoubleType>  thisPt;
    std::vector < std::unordered_map <NodeType*, DoubleType> >  ret;
    
    // isert the first point
    thisPt.insert (std::make_pair(np, 1.0));
    ret.push_back (thisPt);
    
    // next point
    typename MetricSimplexContainerType::PathPointType  path_pt =
            typename MetricSimplexContainerType::PathPointType 
                        (np->CameFromSimplex->w, np->CameFromSimplex->g_came_from_point, thisPt);
    
    if (_dosl_verbose_on(0))  n.print("Reconstruct path called");
    if (_dosl_verbose_on(1))  _dosl_printf("Node in hash: %x.", np);
    while (path_pt.p.size() > 0) {
        if (_dosl_verbose_on(1)) {
            path_pt.print("path point: ");
        }
        ret.push_back (path_pt.p);
        path_pt = AllSimplices.find_camefrom_point (path_pt);
    }
    
    return (ret);
}

template <class AlgDerived, class NodeType, class CostType>
    std::vector<NodeType*>  // convex combined using overloaded operators + and * of 'NodeType'
        SStar::Algorithm<AlgDerived,NodeType,CostType>::reconstruct_pointer_path (NodeType n, bool convertWeightsToInt)
{
    std::vector < std::unordered_map <NodeType*, DoubleType> >  path = reconstruct_weighted_pointer_path (n);
    std::vector<NodeType*>  ret;
    NodeType tn;
    
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
        ret.push_back (all_nodes_set_p->get(tn));
    }
    
    return (ret);
}

// =====================================================================================

#endif
