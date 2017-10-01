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

template < class nodePointerType, class doubleType, class doubleVecType>
class AMetricSimplex;

// ---------------------------

template <class nodePointerType, class doubleType>
class PathPoint
{
public:
    _DOSL_SMALL_MAP <nodePointerType, doubleType>  p; // vertices whose convex combination is this point
    doubleType G; // g-score at point
    void* containing_simplex_p;
    
    _DOSL_SMALL_MAP <nodePointerType, doubleType> last_p;
    
    // default onstructor
    PathPoint (_DOSL_SMALL_MAP <nodePointerType, doubleType> _p = _DOSL_SMALL_MAP <nodePointerType, doubleType>(), 
               doubleType _G = std::numeric_limits<doubleType>::max(), 
               _DOSL_SMALL_MAP <nodePointerType, doubleType> _last_p 
                                        = _DOSL_SMALL_MAP <nodePointerType, doubleType>(),
               void* _containing_simplex_p = NULL ) 
            : p(_p), G(_G), last_p(_last_p), containing_simplex_p(_containing_simplex_p) { }
    
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

template < class nodePointerType, class doubleType=double, class doubleVecType=_DOSL_SMALL_VECTOR<doubleType> >
class MetricSimplexCollection
{
public:

    typedef AMetricSimplex<nodePointerType,doubleType,doubleVecType>  MetricSimplexType;
    typedef std::unordered_set <MetricSimplexType*, MetricSimplexHasher<MetricSimplexType*>, 
                                    MetricSimplexEqualTo<MetricSimplexType*> >  MetricSimplexPointersUnorderedSetType;
    typedef PathPoint<nodePointerType, doubleType>  PathPointType;
    
    // Member variables
    MetricSimplexPointersUnorderedSetType  AllMetricSimplexPointers;
    
    inline size_t size(void) { return (AllMetricSimplexPointers.size()); }
    
    inline MetricSimplexType* find (MetricSimplexType* tmp) { // uses MetricSimplexEqualTo<MetricSimplexType*> to compare
        auto found_it = AllMetricSimplexPointers.find (tmp);
        if (found_it!=AllMetricSimplexPointers.end())
            return (*found_it);
        else
            return (NULL);
    }
    
    // ---------------------------------------
    
    // functions for creating new simplice
    
    MetricSimplexType* getEmptySimplex (void) { return (NULL); }
    
    MetricSimplexType* createNewZeroSimplex (nodePointerType np);
            
    MetricSimplexType* createNewOneSimplex (nodePointerType just_created, nodePointerType came_from, 
                                                        doubleType d=std::numeric_limits<doubleType>::quiet_NaN() );
    
    MetricSimplexType* constructSimplexFromVertices  // Use this as the primary simplex-construction routine
        (_DOSL_SMALL_VECTOR <nodePointerType> nps, 
                unsigned int things_to_compute=COMPUTE_INCREMENTAL_QUANTITIES,
                bool force_recompute=false, int recursion_depth=0);
    
    // function for modifying simplex
    
    MetricSimplexType* checkConnectionsAndAddVertex 
                (MetricSimplexType* inSimplex, nodePointerType np, 
                        unsigned int computation_steps=COMPUTE_INCREMENTAL_QUANTITIES, bool force_recompute=false);
    
    // Maximal simplex construction
    
    std::unordered_set<nodePointerType> getAllCommonNeighborsOfSimplex (MetricSimplexType* inSimplex_p);
    
    std::unordered_set <MetricSimplexType*> getAllMaximalSimplicesFromSet 
        (std::unordered_set<nodePointerType>* nodeSet_p=NULL, // Will call getAllCommonNeighborsOfSimplex if NULL
            MetricSimplexType* baseSimplex_p=NULL, unsigned int things_to_compute=COMPUTE_ALL,
            std::unordered_set<nodePointerType> restricted_neighbor_set = std::unordered_set<nodePointerType>() );
    
    std::unordered_set <MetricSimplexType*> getAllAttachedMaximalSimplices 
        (MetricSimplexType* inSimplex_p,
            unsigned int things_to_compute=COMPUTE_ALL, // default: COMPUTE_ALL,
            bool return_base_if_maximal=true, // default: true
            bool forceCompute=false, // default: false
            bool use_only_expnded_nodes=false, // default: false
                std::unordered_set<nodePointerType>* neighbor_node_pointers_p=NULL);
    
    // ---------------------------------------
    // for path reconstruction
    PathPointType  findCamefromPoint (const PathPointType& inPathPoint);
    
};

#endif

