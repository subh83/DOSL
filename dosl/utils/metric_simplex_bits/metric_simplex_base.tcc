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

#ifndef __METRIC_SIMPLEX_BASE_TCC_
#define __METRIC_SIMPLEX_BASE_TCC_

#include <set>

#include "../macros_constants.tcc"
#include "../basic_math.tcc"

#define _MS_SMALL_VECTOR  std::vector  // _DOSL_SMALL_VECTOR //std::vector

// ===========================================================

// Comparison and hasher classes for metric simplex pointers stored in a big heap

enum SimplexCompareType {
    // Consider vertex 0 as special: bit 0
        POINTED = 1u,
    // Compare G-score as well: bit 1
        GSCORE = 2u
};

// --------------------------

template <class metricSimplexPointerType, unsigned int compareType=(POINTED|GSCORE)>
class MetricSimplexHasher
{
public:
    size_t operator()(const metricSimplexPointerType& sp) const {
        int compare_start_index = 0;
        size_t hash_val = 0u;
        if (compareType & POINTED) {
            hash_val = (((size_t)(sp->p[0])) << 8) * sp->n_vertices;
            compare_start_index = 1;
        }
        for (int a=compare_start_index; a<sp->p.size(); ++a) 
            hash_val += (size_t)(sp->p[a]);
        return (hash_val);
    }
};

template <class metricSimplexPointerType, unsigned int compareType=(POINTED|GSCORE)>
class MetricSimplexEqualTo
{
public:
    bool operator()(const metricSimplexPointerType& sp1, const metricSimplexPointerType& sp2) const {
        int i;
        // check vertex count
        if (sp1->p.size() != sp2->p.size())  return (false);
        int compare_start_index = 0;
        if (compareType & POINTED) {
            // Check vertex 0 identities
            if (sp1->p[0] != sp2->p[0])  return (false);
            compare_start_index = 1;
        }
        for (int a=compare_start_index; a<sp1->p.size(); ++a) {
            // Check vertex identities
            i = sp2->p.findi (sp1->p[a]);
            if (i<0) 
                return (false);
            // Check g-score
            if ( (compareType & GSCORE) && (fabs(sp2->gs[i]-sp1->gs[a]) > _MS_DOUBLE_EPS) ) return (false);
        }
        // No need to check edge distances (they will always be same)
        return (true);
    }
};

// --------------------------

template <class pair_type>
class CompareBySecond
{
public:
    bool operator()(const pair_type& p1, const pair_type& p2) const {
        if (p1.second==p2.second) return (p1.first<p2.first);
        else return (p1.second<p2.second);
    }
};

template <class F, class S>
using set_of_pairs_ordered_by_second = std::set < std::pair<F,S>, CompareBySecond< std::pair<F,S> > >; // c++11

// -----------------------------------------------------------------------------------------------

template < class doubleType=double, class doubleVecType=std::vector<doubleType> >
doubleType DistanceBetweenVertices (doubleVecType& v1, doubleVecType& v2)
{
    doubleType sum = 0.0;
    doubleType diff;
    for (int a=0; a<v1.size(); ++a) {
        diff = v2[a] - v1[a];
        sum += diff*diff;
    }
    return (sqrt(sum));
}

// ===========================================================

enum SIMPLEX_VERTEX_INSERTION_STEPS {
    NO_COMPUTATION = 0u,
    CHECK_CONNECTION = 1u, COMPUTE_LOCAL_COORDINATE = 2u, // bits 0 & 1
    SET_G_SCORE_OF_VERTEX = 4u, COMPUTE_COORDINATE_OF_O = 8u, // bits 2 & 3
    COMPUTE_G_SCORE_OF_APEX = 16u, // bit 4
    // compositions:
    EMBED_SIMPLEX = (CHECK_CONNECTION | COMPUTE_LOCAL_COORDINATE),
    EMBED_O = (CHECK_CONNECTION | SET_G_SCORE_OF_VERTEX | COMPUTE_COORDINATE_OF_O),
    EMBED_ALL = (EMBED_SIMPLEX | EMBED_O),
    DONT_COMPUTE_COORDINATE_OF_O = (CHECK_CONNECTION | COMPUTE_LOCAL_COORDINATE | SET_G_SCORE_OF_VERTEX),
    DONT_COMPUTE_G_SCORE_OF_APEX = EMBED_ALL,
    COMPUTE_DEFAULTS = EMBED_ALL,
    COMPUTE_INCREMENTAL_QUANTITIES = EMBED_ALL,
    COMPUTE_INCREMENTAL_QUANTITIES_EXCEPT_O = EMBED_SIMPLEX,
    COMPUTE_ALL = (EMBED_ALL | COMPUTE_G_SCORE_OF_APEX)
};

#endif
