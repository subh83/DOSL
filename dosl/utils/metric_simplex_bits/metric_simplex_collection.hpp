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

#ifndef __DOSL_METRIC_SIMPLEX_COLLECTION_HPP__ 
#define __DOSL_METRIC_SIMPLEX_COLLECTION_HPP__ 

#include "metric_simplex_base.tcc"

template < class _NodePointerType, class DoubleType, class DoubleVecType>
class AMetricSimplex;

// ---------------------------

template <class _NodePointerType, class DoubleType>
class PathPoint
{
public:
    _DOSL_SMALL_MAP <_NodePointerType, DoubleType>  p; // vertices whose convex combination is this point
    DoubleType g_score; // g-score at point
    void* containing_simplex_p;
    
    _DOSL_SMALL_MAP <_NodePointerType, DoubleType> last_p;
    
    // default onstructor
    PathPoint (_DOSL_SMALL_MAP <_NodePointerType, DoubleType> _p = _DOSL_SMALL_MAP <_NodePointerType, DoubleType>(), 
               DoubleType _G = std::numeric_limits<DoubleType>::max(), 
               _DOSL_SMALL_MAP <_NodePointerType, DoubleType> _last_p 
                                        = _DOSL_SMALL_MAP <_NodePointerType, DoubleType>(),
               void* _containing_simplex_p = NULL ) 
            : p(_p), g_score(_G), last_p(_last_p), containing_simplex_p(_containing_simplex_p) { }
    
    void print (std::string head="", std::string tail="") {
        _dosl_cout << _GREEN + head << " (" << this << ")" GREEN_ ": " << _dosl_endl;
        for (auto it=p.begin(); it!=p.end(); ++it) {
            it->first->print();
            _dosl_printf ("\tweight = %f", it->second);
        }
        _dosl_cout << tail << _dosl_endl;
    }
};

// ---------------------------

template < class _NodePointerType, class DoubleType=double, class DoubleVecType=_DOSL_SMALL_VECTOR<DoubleType> >
class MetricSimplexCollection
{
public:

    typedef AMetricSimplex<_NodePointerType,DoubleType,DoubleVecType>  MetricSimplexType;
    typedef std::unordered_set <MetricSimplexType*, MetricSimplexHasher<MetricSimplexType*>, 
                                    MetricSimplexEqualTo<MetricSimplexType*> >  MetricSimplexPointersUnorderedSetType;
    typedef PathPoint<_NodePointerType, DoubleType>  PathPointType;
    
    // Member variables
    MetricSimplexPointersUnorderedSetType  all_metric_simplex_pointers;
    
    inline size_t size(void) { return (all_metric_simplex_pointers.size()); }
    
    inline MetricSimplexType* find (MetricSimplexType* tmp) { // uses MetricSimplexEqualTo<MetricSimplexType*> to compare
        auto found_it = all_metric_simplex_pointers.find (tmp);
        if (found_it!=all_metric_simplex_pointers.end())
            return (*found_it);
        else
            return (NULL);
    }
    
    // ---------------------------------------
    
    // functions for creating new simplice
    
    MetricSimplexType* get_empty_simplex (void) { return (NULL); }
    
    MetricSimplexType* create_new_zero_simplex (_NodePointerType np);
            
    MetricSimplexType* create_new_one_simplex (_NodePointerType just_created, _NodePointerType came_from, 
                                                        DoubleType d=std::numeric_limits<DoubleType>::quiet_NaN() );
    
    MetricSimplexType* construct_simplex_from_vertices  // Use this as the primary simplex-construction routine
        (_DOSL_SMALL_VECTOR <_NodePointerType> nps, 
                unsigned int things_to_compute=COMPUTE_INCREMENTAL_QUANTITIES,
                bool force_recompute=false, int recursion_depth=0);
    
    // function for modifying simplex
    
    MetricSimplexType* check_connections_and_add_vertex 
                (MetricSimplexType* inSimplex, _NodePointerType np, 
                        unsigned int computation_steps=COMPUTE_INCREMENTAL_QUANTITIES, bool force_recompute=false);
    
    // Maximal simplex construction
    
    std::unordered_set<_NodePointerType> get_all_common_neighbors_of_simplex (MetricSimplexType* in_simplex_p);
    
    std::unordered_set <MetricSimplexType*> getAllMaximalSimplicesFromSet 
        (std::unordered_set<_NodePointerType>* node_set_p=NULL, // Will call get_all_common_neighbors_of_simplex if NULL
            MetricSimplexType* baseSimplex_p=NULL, unsigned int things_to_compute=COMPUTE_ALL,
            std::unordered_set<_NodePointerType> restricted_neighbor_set = std::unordered_set<_NodePointerType>() );
    
    std::unordered_set <MetricSimplexType*> get_all_attached_maximal_simplices 
        (MetricSimplexType* in_simplex_p,
            unsigned int things_to_compute=COMPUTE_ALL, // default: COMPUTE_ALL,
            bool return_base_if_maximal=true, // default: true
            bool force_compute=false, // default: false
            bool use_only_expnded_nodes=false, // default: false
                std::unordered_set<_NodePointerType>* neighbor_node_pointers_p=NULL);
    
    // ---------------------------------------
    // for path reconstruction
    PathPointType  find_camefrom_point (const PathPointType& in_path_point);
    
};

#endif

