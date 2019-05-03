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

#ifndef __DOSL_AStar_TCC
#define __DOSL_AStar_TCC
// user-readable
#define DOSL_ALGORITHM_AStar

// includes

#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <limits>
#include <type_traits>
#include <string>
#include <memory>

#include "../utils/macros_constants.tcc"
#include "../utils/stl_utils.tcc"
#include "planner_bits.hpp"

// ====================================================================

class AStar {
public:
    declare_alg_name("AStar"); // macro from '_planner_bits'

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
    template <class NodeType, class _CostType=double> // CRTP
    class Node : public BinaryHeapElement
    {
    public:
        typedef _CostType  CostType;
        
        // Utility for looking up algorithm name:
        declare_alg_name("AStar"); // macro from '_planner_bits'
        
        // For keeping track of hash insertion (new node creation)
        bool post_hash_insert_initiated;
        
        // Node specific variables
        CostType f_score, g_score;
        bool expanded; // Whether in closed list or not
        LineageDataType lineage_data; // stores which start node the node came from
        
        // successors
        _DOSL_SMALL_MAP<NodeType*,CostType> successors;
        bool successors_created;
        
        // for reconstruction
        NodeType* came_from;
        
        // -------------------------------------
        // constructors
        Node() : post_hash_insert_initiated(false),
                      expanded(false), successors_created(false), came_from(NULL), lineage_data(LineageDataType()) { }
        
        void clear_search_data (unsigned int mode = CLEAR_NODE_SUCCESSORS) {
                    // to be used in getSuccessors. TODO: do this in copy constructor and operator== instead
            post_hash_insert_initiated = false; expanded = false; came_from = NULL;
            heap_pos = -1; g_score = (CostType)0.0; f_score = (CostType)0.0;
            if (mode & CLEAR_NODE_SUCCESSORS) { successors.clear(); successors_created = false; }
            if (mode & CLEAR_NODE_LINEAGE) lineage_data = LineageDataType();
        }
        
        // Define virtual functions of derived class
        inline CostType HeapKey() { return (f_score); }
        
        // ----------------------------------------------------------------------
        // Functions to be overwritten by user node type
        //      (need not be virtual since use is of only derived class members).
        // Need to have virtual members with same name in the problem class.
        inline int getHashBin (void) { return (0); }
        bool operator==(const NodeType& n)
            { _dosl_default_fun_warn("'AStar::Node::operator==' OR 'AStar::Algorithm::equalTo'"); }
        inline void getSuccessors (std::vector<NodeType>* s, std::vector<CostType>* c)
            { _dosl_default_fun_warn("AStar::[Algorithm|Node]::getSuccessors"); }
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
        inline CostType getHeapKey (double subopEps) { return (g_score + (CostType)(subopEps*getHeuristics())); }
    };


    template <class AlgDerived, class NodeType, class _CostType=double>  // CRTP
    class Algorithm
    {
    public:
        typedef _CostType  CostType;
        NodeType dummy_node; // forces compiler to generate code for the possibly template class 'NodeType'
        
        // Utility for looking up algorithm name:
        declare_alg_name("AStar"); // macro from '_planner_bits'
        
        // Constants
        #if _DOSL_EVENTHANDLER
        enum NodeEventType {
            // expanded: bit 0
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
        
        // user definable parameters for problem
        double subopt_eps;
        
        // counters
        int expand_count;
        // only for verbode:
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
            bool operator()(NodeType& n1, NodeType& n2) { return data_p->equalTo(n1,n2); }
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
        
        // =====================================================
        // ---------------------------
        // Constructors and initiators
        Algorithm (AlgDerived* shared_instance_p = NULL) : _this (static_cast<AlgDerived*>(this))
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
            
            subopt_eps = 1.0;
            progress_show_interval = 10000;
        }
        
        // ---------------------------
        // Main search functions
        void search (void);
        void clear (unsigned int mode = CLEAR_NODE_DATA | CLEAR_NODES_HEAP);
        
        // Functions for reading paths to arbitrary nodes
        std::vector < std::unordered_map <NodeType*, CostType> > reconstruct_weighted_pointer_path (NodeType n); // for compatibility with SStar
        std::vector<NodeType*> reconstruct_pointer_path (NodeType n);
        // access other node data
        inline NodeType* get_node_pointer (NodeType n) { return (all_nodes_set_p->get(n)); }
        inline CostType get_costs_to_nodes (NodeType n) { return (all_nodes_set_p->get(n)->g_score); }
        // bookmark nodes
        inline std::vector<NodeType*> get_bookmark_node_pointers (void) { return (bookmarked_node_pointers); }
        
        // =====================================================
        // -----------------------------------------------------
        // functions to be overwritten by user problem instance. Should be called with "_this->"
        // Also in node class
        unsigned int getHashBin (NodeType &n) { return (n.getHashBin()); }
        bool equalTo (NodeType &n1, NodeType &n2) { return (n1==n2); }
        void getSuccessors (NodeType &n, std::vector<NodeType>* const s, std::vector<CostType>* const c) 
                                            { return (n.getSuccessors(s,c)); }
        inline CostType getHeuristics (NodeType& n) { return (n.getHeuristics()); }
        inline std::vector<NodeType> getStartNodes (void)
            { _dosl_default_fun_warn("AStar::Algorithm::getStartNodes"); return (std::vector<NodeType>()); }
        inline bool bookmarkNode (NodeType &n) { return (n.bookmarkNode()); }
        inline bool stopSearch (NodeType &n) { return (n.stopSearch()); }
        #if _DOSL_EVENTHANDLER
        inline void nodeEvent (NodeType &n, unsigned int e) { n.nodeEvent(e); }
        #endif
        inline virtual void print (NodeType &n, std::string head) { n.print(head); }
        // Derived functions. Can also be directly overwritten.
        CostType getHeapKey (NodeType& n) { return (n.g_score + (CostType)(subopt_eps * _this->getHeuristics(n))); }
        
        // -------------------------------------------------
    private:
        AlgDerived* _this;
        
        // Derived and helper functions
        void generate_successors (NodeType* node_in_hash_p);
        
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
void AStar::Algorithm<AlgDerived,NodeType,CostType>::generate_successors (NodeType* node_in_hash_p)
{
    if ( !(node_in_hash_p->successors_created) ) // successors were not generated previously
    {
        this_successors.clear();
        this_transition_costs.clear();
        _this->getSuccessors (*node_in_hash_p, &this_successors, &this_transition_costs);
        
        #if _DOSL_DEBUG > 0
        if (this_successors.size()!=this_transition_costs.size())
            _dosl_err("Number of successors (%d) is different from numer of transition costs (%d) as returned by 'getSuccessors'.", this_successors.size(), this_transition_costs.size());
        #endif
        
        node_in_hash_p->successors.reserve (this_successors.size()); // reserve space for fast pushing
        for (a=0; a<this_successors.size(); a++) {
            this_successors[a].clear_search_data (CLEAR_NODE_SUCCESSORS); // in case the successors were created from copies of this
            // note: liniage data is set in the main search function
            node_in_hash_p->successors.insert (
                _DOSL_SMALL_MAP_pairfun (all_nodes_set_p->get(this_successors[a]), this_transition_costs[a]) );
        }
        node_in_hash_p->successors_created = true;
    }
}

// -------------------------------------------------------------------------------------

template <class AlgDerived, class NodeType, class CostType>
void AStar::Algorithm<AlgDerived,NodeType,CostType>::search (void)
{
    _dosl_verbose_head(1);
    
    start_nodes = _this->getStartNodes();
    
    #if _DOSL_DEBUG > 0
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
    for (int a=0; a<start_nodes.size(); a++) {
        start_nodes[a].clear_search_data (CLEAR_NODE_SUCCESSORS | CLEAR_NODE_LINEAGE); // in case this node is output of a previous search
        // put start node in hash
        thisNodeInHash_p = all_nodes_set_p->get (start_nodes[a]);
        // set member variables for start node
        thisNodeInHash_p->came_from = NULL;
        if (!(thisNodeInHash_p->lineage_data.is_set()))
            thisNodeInHash_p->lineage_data = LineageDataType(a); // lineage a, generation 0
        thisNodeInHash_p->expanded = false;
        thisNodeInHash_p->g_score = (CostType)0.0;
        thisNodeInHash_p->f_score = _this->getHeapKey(*thisNodeInHash_p); // TODO: this is not being used!
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
            if (progress_show_interval>0) {
                if (expand_count % progress_show_interval == 0) {
                    _dosl_printf("Number of states expanded: %d. Heap size: %d. Time elapsed: %f_score s.", 
                            expand_count, node_heap_p->size(), timer.read());
                }
            }
        expand_count++;
        
        // get the node with least f_score-value
        thisNodeInHash_p = node_heap_p->pop();
        if (_dosl_verbose_on(1)) {
            thisNodeInHash_p->print("Now expanding: ");
            _dosl_printf("g_score-value at expanding node: %f. f_score-value at expanding node: %f.", thisNodeInHash_p->g_score, thisNodeInHash_p->f_score);
        }
        
        // Generate the neighbours if they are already not generated
        generate_successors (thisNodeInHash_p);
        if (_dosl_verbose_on(1)) {
            _dosl_printf("Number of successors = %d.", thisNodeInHash_p->successors.size());
        }
        
        // -----------------------------------------------------
        // Expand
        
        thisNodeInHash_p->expanded = true; // Put in closed list
        
        #if _DOSL_EVENTHANDLER
        _this->nodeEvent (*thisNodeInHash_p, EXPANDED|POPPED);
        #endif
        
        // Check if we need to bookmark the node being Expanded
        if ( _this->bookmarkNode (*thisNodeInHash_p) )
        {
            bookmarked_node_pointers.push_back (thisNodeInHash_p);
            if (_dosl_verbose_on(0)) {
                thisNodeInHash_p->print ("Bookmarked a node: ");
                _dosl_printf("... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                        expand_count, node_heap_p->size(), timer.read());
            }
        }
        // Check if we need to stop furthur expansion
        if ( _this->stopSearch (*thisNodeInHash_p) ) {
            if (_dosl_verbose_on(0))
                if (progress_show_interval>0) {
                    thisNodeInHash_p->print ("Stopping search ('stopSearch' returned true) at: ");
                    _dosl_printf("... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                            expand_count, node_heap_p->size(), timer.read());
                }
            return;
        }
        
        // Initiate the neighbours (if required) and update their g_score & f_score values
        for (auto it = thisNodeInHash_p->successors.begin(); it!=thisNodeInHash_p->successors.end(); ++it)
        {
            this_neighbour_node_in_hash_p = it->first;
            thisTransitionCost = it->second;
            
            // An uninitiated neighbour node - definitely g_score & f_score values not set either.
            if (!(this_neighbour_node_in_hash_p->post_hash_insert_initiated)) {
                this_neighbour_node_in_hash_p->came_from = thisNodeInHash_p;
                this_neighbour_node_in_hash_p->lineage_data = thisNodeInHash_p->lineage_data.next_generation();
                this_neighbour_node_in_hash_p->g_score = thisNodeInHash_p->g_score + thisTransitionCost;
                this_neighbour_node_in_hash_p->f_score = _this->getHeapKey (*this_neighbour_node_in_hash_p);
                this_neighbour_node_in_hash_p->expanded = false;
                // Put in open list and continue to next neighbour
                this_neighbour_node_in_hash_p->post_hash_insert_initiated = true;
                                        // ^^^ Always set this when other variables have already been set
                
                if (_dosl_verbose_on(2)) {
                    this_neighbour_node_in_hash_p->print("Child (generated for the first time): ");
                    _dosl_printf("Transition cost to child: %f.", thisTransitionCost);
                }
                
                // push
                node_heap_p->push (this_neighbour_node_in_hash_p);
                #if _DOSL_EVENTHANDLER
                _this->nodeEvent (*this_neighbour_node_in_hash_p, PUSHED);
                #endif
                
                continue;
            }
            
            if (_dosl_verbose_on(2)) {
                this_neighbour_node_in_hash_p->print("Child (generated previously): ");
                _dosl_printf("Transition cost to child: %f.", thisTransitionCost);
            }
            
            // Neighbour that is in closed list is to be skipped
            if (this_neighbour_node_in_hash_p->expanded)
                continue;
            
            // Update came_from, g_score and f_score values if better
            test_g_val = thisNodeInHash_p->g_score + thisTransitionCost;
            if (test_g_val < this_neighbour_node_in_hash_p->g_score) {
                this_neighbour_node_in_hash_p->g_score = test_g_val;
                this_neighbour_node_in_hash_p->f_score = _this->getHeapKey (*this_neighbour_node_in_hash_p);
                this_neighbour_node_in_hash_p->came_from = thisNodeInHash_p;
                this_neighbour_node_in_hash_p->lineage_data = thisNodeInHash_p->lineage_data.next_generation();
                
                // Since thisNeighbourGraphNode->f_score is changed, re-arrange it in heap
                node_heap_p->update (this_neighbour_node_in_hash_p);
                #if _DOSL_EVENTHANDLER
                _this->nodeEvent (*this_neighbour_node_in_hash_p, UPDATED);
                #endif
            }
        }
    }
    
    if (_dosl_verbose_on(0)) {
        if (progress_show_interval>0 && node_heap_p->empty()) {
            _dosl_printf("Stopping search!! Heap is empty... Number of states expanded: %d. Heap size: %d. Time elapsed: %f s.", 
                       expand_count, node_heap_p->size(), timer.read());
        }
    }
}

// -------------------------------------------------------------------------------------

template <class AlgDerived, class NodeType, class CostType>
void AStar::Algorithm<AlgDerived,NodeType,CostType>::clear (unsigned int mode)
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
}

// -------------------------------------------------------------------------------------
// For reading results:

template <class AlgDerived, class NodeType, class CostType>
std::vector<NodeType*> AStar::Algorithm<AlgDerived,NodeType,CostType>::reconstruct_pointer_path (NodeType n)
{
    _dosl_verbose_head(1);
        
    std::vector<NodeType*> thisPath;
    // Reconstruct path
    NodeType* thisNodeInHash_p = all_nodes_set_p->get(n);
    if (_dosl_verbose_on(0))  n.print("Reconstruct path called");
    if (_dosl_verbose_on(1))  _dosl_printf("Node in hash: %x.", thisNodeInHash_p);
    while (thisNodeInHash_p) {
        if (_dosl_verbose_on(1)) thisNodeInHash_p->print();
        thisPath.push_back (thisNodeInHash_p);
        thisNodeInHash_p = thisNodeInHash_p->came_from;
    }
    return (thisPath);
}

// ===========================
// For compatibity with s-star

template <class AlgDerived, class NodeType, class CostType>
    std::vector < std::unordered_map <NodeType*, CostType> >
        AStar::Algorithm<AlgDerived,NodeType,CostType>::reconstruct_weighted_pointer_path (NodeType n)
{
    std::vector<NodeType*> npVec = reconstruct_pointer_path (n);
    
    std::vector < std::unordered_map <NodeType*, CostType> >  ret;
    for (auto it=npVec.begin(); it!=npVec.end(); ++it) {
        std::unordered_map <NodeType*, CostType> thisPt;
        thisPt[(*it)] = 1.0;
        ret.push_back (thisPt);
    }
    return (ret);
}

#endif
