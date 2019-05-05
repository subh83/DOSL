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

#ifndef __DOSL_emptyplanner_TCC
#define __DOSL_emptyplanner_TCC
// user-readable
#define DOSL_ALGORITHM_emptyplanner

// includes

#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <limits>
#include <string>
#include <memory>

#include "../utils/macros_constants.tcc"
#include "../utils/stl_utils.tcc"
#include "planner_bits.hpp"

/* Naming Conventions:
    User accessible:
    - Class or type names:                          'AbcXyzEfg'
    - Template class/typename parameters:           'abcXyzEfg'
    
    - Member variables:                             'abcXyzEfg'
    - Temporary/internal/local variables:           'abc_xyz_efg'
    
    - Member functions that user have access to:    'abcXyzEfg'
    - Internal functions that user will not use:    'abc_xyz_efg'
*/

class emptyplanner {
public:
    declare_alg_name("emptyplanner"); // macro from '_planner_bits'

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
        declare_alg_name("emptyplanner"); // macro from '_planner_bits'
        
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
        // pseudo-destructor
        void clear_search_data (unsigned int mode = CLEAR_NODE_SUCCESSORS);
        
        // ----------------------------------------------------------------------
        // Functions to be overwritten by user node type.
        // If functions in Algorithm class are overwritten, these will be ignored.
        inline int getHashBin (void) { return (0); }
        bool operator==(const NodeType& n)
            { _dosl_default_fun_warn("'emptyplanner::Node::operator==' OR 'emptyplanner::Algorithm::equalTo'"); }
        inline void getSuccessors (std::vector<NodeType>* s, std::vector<CostType>* c)
            { _dosl_default_fun_warn("emptyplanner::[Algorithm|Node]::getSuccessors"); }
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
        declare_alg_name("emptyplanner"); // macro from '_planner_bits'
        
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
        Algorithm (AlgDerived* shared_instance_p = NULL);
        
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
        // Functions to be overwritten by user problem instance. Should be called with "_this->"
        // The following defaults use the members of the node class.
        unsigned int getHashBin (NodeType &n) { return (n.getHashBin()); }
        bool equalTo (NodeType &n1, NodeType &n2) { return (n1==n2); }
        void getSuccessors (NodeType &n, std::vector<NodeType>* const s, std::vector<CostType>* const c) 
                                            { return (n.getSuccessors(s,c)); }
        inline CostType getHeuristics (NodeType& n) { return (n.getHeuristics()); }
        inline std::vector<NodeType> getStartNodes (void)
            { _dosl_default_fun_warn("emptyplanner::Algorithm::getStartNodes"); return (std::vector<NodeType>()); }
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

#endif
