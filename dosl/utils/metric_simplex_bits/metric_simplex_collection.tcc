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

#ifndef __DOSL_METRIC_SIMPLEX_COLLECTION_TCC__ 
#define __DOSL_METRIC_SIMPLEX_COLLECTION_TCC__ 

#include "metric_simplex_collection.hpp"

// ===================================================================
// ===================================================================
// MetricSimplexCollection members

// ----------------------------------------
// functions for creating new simplice

// Wrappers around constructors of AMetricSimplex:

template <class _NodePointerType, class DoubleType, class DoubleVecType>
    AMetricSimplex<_NodePointerType,DoubleType,DoubleVecType>* 
        MetricSimplexCollection<_NodePointerType,DoubleType,DoubleVecType>::create_new_zero_simplex 
            (_NodePointerType np) 
{
    MetricSimplexType* ret = new MetricSimplexType (np, this);
    
    // vs, vs_sq are empty; dist_mat and dist_sq_mat are 1x1 containg 0; w has one element with 0; -- no need to set these
    MetricSimplexType* found_p = find (ret); // check if simplex was already created.
    if (found_p) { delete ret; return (found_p); }
    
    all_metric_simplex_pointers.insert (ret);
    return (ret);
}


template <class _NodePointerType, class DoubleType, class DoubleVecType>
    AMetricSimplex<_NodePointerType,DoubleType,DoubleVecType>* 
        MetricSimplexCollection<_NodePointerType,DoubleType,DoubleVecType>::create_new_one_simplex 
            (_NodePointerType just_created, _NodePointerType came_from, DoubleType d) 
{
    MetricSimplexType* new_1_simplex = new MetricSimplexType (just_created, came_from, d, this);
    
    // check if simplex was already created. TODO: Make more efficient by chcking earlier
    MetricSimplexType* found_p = find (new_1_simplex);
    if (found_p) {
        delete new_1_simplex;
        return (found_p);
    }
    
    all_metric_simplex_pointers.insert (new_1_simplex);
    return (new_1_simplex);
}

// ---------------------------------------
// Incremental construction of simplex

template <class _NodePointerType, class DoubleType, class DoubleVecType>
    AMetricSimplex<_NodePointerType,DoubleType,DoubleVecType>* 
        MetricSimplexCollection<_NodePointerType,DoubleType,DoubleVecType>::check_connections_and_add_vertex 
            (AMetricSimplex<_NodePointerType,DoubleType,DoubleVecType>* in_simplex_p, _NodePointerType np, 
                unsigned int things_to_compute, bool force_recompute)
{
    _dosl_verbose_head(1);
    
    if (_dosl_verbose_on(0)) {
        in_simplex_p->print (_YELLOW "inSimplex" YELLOW_);
        np->print (_YELLOW "inNode" YELLOW_);
        _dosl_printf_nobreak ("check_connections_and_add_vertex: Requested things_to_compute=%s. computed: ", 
                                _uint_to_binary(things_to_compute).c_str());
    }
    
    // If inSimplex is an empty simplex (NULL)
    if (!in_simplex_p)
        return (create_new_zero_simplex (np));
    
    // need to return a copy (calls copy constryctors of DoubleType)
    MetricSimplexType* ret = new MetricSimplexType (*in_simplex_p);
    
    // Make room
    ret->MakeRoomForAnotherVertex (np);
    
    // Need to reset several variables in ret.
    // inherit: n_vertices, n_vertices_m1, n_vertices_m2, gs, AllSimplices_p
    ret->back_track_simplex = 0;
    for (int a=0; a<ret->face_simplices.size(); ++a)
        ret->face_simplices[a] = NULL;
    ret->face_simplices [ret->n_vertices-1] = in_simplex_p; // Face opposite to last vertex inserted.
    ret->simplex_computation_stage = 0u;
    ret->apex_came_from_face = NULL; // need to recompute after inserting the vertex
    ret->is_maximal = false;
    
    // clear stored data (TODO: can we make some of these incremental?
    // ret->face_simplices.clear(); // already taken care of
    ret->tw.clear();
    ret->w.clear();
    // TODO: The following are subsets of the quantities in in_simplex_p.
    ret->attached_absolute_maximal_simplices.clear();
    ret->all_common_neighbors.clear();
    
    
    if ( TO_BOOL(things_to_compute & CHECK_CONNECTION) ) {
        _MS_SMALL_VECTOR<DoubleType> dist_to_vs;
        if (TO_BOOL(in_simplex_p->simplex_computation_stage & CHECK_CONNECTION) &&
                in_simplex_p->checkConnections (np, &dist_to_vs, /*dist_to_p0,*/ true)) { // checks np->successors.
            ret->SetDistancesFromLastVertex (dist_to_vs);
            ret->simplex_computation_stage |= CHECK_CONNECTION;
            if (_dosl_verbose_on(0)) { printf ("CHECK_CONNECTION,"); }
        }
        else {
            if (_dosl_verbose_on(0)) {
                _dosl_printf (_RED "failed in 'CHECK_CONNECTION'." RED_ " in_simplex_p->simplex_computation_stage=%s(%d)", uint_to_binary(in_simplex_p->simplex_computation_stage), in_simplex_p->simplex_computation_stage);
                _dosl_printf ("np->successors.size()=%d", np->successors.size());
            }
            ret->simplex_computation_failure = CHECK_CONNECTION;
            return (ret);
        }
    }
    
    // Set g_score-score
    if (TO_BOOL(things_to_compute & SET_G_SCORE_OF_VERTEX)) {
        if (TO_BOOL(in_simplex_p->simplex_computation_stage & SET_G_SCORE_OF_VERTEX) &&
                    np->g_score!=std::numeric_limits<DoubleType>::max() ) {
            ret->SetGscoreOfLastVertex (np->g_score);
            ret->simplex_computation_stage |= SET_G_SCORE_OF_VERTEX;
            if (_dosl_verbose_on(0)) { printf("SET_G_SCORE_OF_VERTEX,"); }
        }
        else { // np->g_score is either not set or is invalid.
            if (_dosl_verbose_on(0)) {
                _dosl_printf ("failed in 'SetGscoreOfLastVertex'. in_simplex_p->simplex_computation_stage=%s(%d), np->g_score=%f", 
                        uint_to_binary(in_simplex_p->simplex_computation_stage), in_simplex_p->simplex_computation_stage, np->g_score);
            }
            ret->simplex_computation_failure = SET_G_SCORE_OF_VERTEX;
            return (ret);
        }
    }
    
    // +++++++++++++++++
    // Search
    if (!force_recompute) {
        MetricSimplexType* found_p = find (ret); // check if simplex was already created.
        if (found_p  &&  found_p->simplex_computation_stage==things_to_compute) {
            delete ret;
            found_p->CompleteComputation (things_to_compute); // TODO...
            if (_dosl_verbose_on(0)) { _dosl_printf("**found (%x) - stopping.", found_p); }
            return (found_p);
        }
        if (_dosl_verbose_on(0)) { _dosl_printf("(not found)"); }
    }
    // +++++++++++++++++
    
    // Phase 1: Compute embedding of vertices
    // ----  
    
    if (TO_BOOL(things_to_compute & COMPUTE_LOCAL_COORDINATE /*SET_G_SCORE_OF_VERTEX*/)) {
        if (TO_BOOL(in_simplex_p->simplex_computation_stage & COMPUTE_LOCAL_COORDINATE)) {
            if (!(ret->ComputeCoordinateOfLastVertex())) { // degenerate simplex
            	if (_dosl_verbose_on(0)) {
                	_dosl_warn("COMPUTE_LOCAL_COORDINATE: Degeneracy in metric simplex. Failed.");
                	ret->print("degenerate simplex (incomplete): ");
            	}
                ret->simplex_computation_failure = COMPUTE_LOCAL_COORDINATE;
                return (ret);
            }
            ret->simplex_computation_stage |= COMPUTE_LOCAL_COORDINATE;
            if (_dosl_verbose_on(0)) { printf("COMPUTE_LOCAL_COORDINATE,"); }
        }
        else {
            if (_dosl_verbose_on(0)) {
            	_dosl_warn("check_connections_and_add_vertex: Cannot compute local coordinates since D-scores not set.\n");
            	ret->print("degenerate simplex (incomplete): ");
        	}
        }
    }
    
    // Stage II computation:
    /* Compute coordinate of O - requires g_score of all other vertices - can changes */
    if ( TO_BOOL(things_to_compute & COMPUTE_COORDINATE_OF_O)) {
        if (TO_BOOL(in_simplex_p->simplex_computation_stage & COMPUTE_COORDINATE_OF_O) ) {
            // Need to initiate the C's if inSimplex was a 0-simplex
            if (in_simplex_p->p.size() == 1) {
                ret->C2 = 1.0;
                ret->C1 = -2 * ret->vs[1][0];
                ret->C0 = ret->vs[1][0] * ret->vs[1][0] - ret->gsSq[1];
            }
            // compute coordinate of O
            if (!(ret->ComputeCoordinateOfO())) { // degenerate simplex
            	#if _DOSL_DEBUG > 1
            	_dosl_warn("check_connections_and_add_vertex: Degeneracy in computing embedding of 'o'.\n");
            	#endif
                ret->simplex_computation_failure = COMPUTE_COORDINATE_OF_O;
                return (ret);
            }
            ret->simplex_computation_stage |= COMPUTE_COORDINATE_OF_O;
            if (_dosl_verbose_on(0)) { printf("COMPUTE_COORDINATE_OF_O,"); }
        }
    }
    
    /* Note: Typically will need to call g_score-score computation separately. */
    
    if (TO_BOOL(things_to_compute & COMPUTE_G_SCORE_OF_APEX)) { // won't be true in general. non-incremental
        if (!(ret->ComputeGScore ())) {
            ret->simplex_computation_failure = COMPUTE_G_SCORE_OF_APEX;
            return (ret);
        }
        ret->simplex_computation_stage |= COMPUTE_G_SCORE_OF_APEX;
        if (_dosl_verbose_on(0)) { printf("COMPUTE_G_SCORE_OF_APEX,"); }
    }
    
    all_metric_simplex_pointers.insert (ret);
    return (ret);
}

// ===================================================

template <class _NodePointerType, class DoubleType, class DoubleVecType>
    AMetricSimplex<_NodePointerType,DoubleType,DoubleVecType>* 
        MetricSimplexCollection<_NodePointerType,DoubleType,DoubleVecType>::construct_simplex_from_vertices 
            (_DOSL_SMALL_VECTOR <_NodePointerType> nps, 
                unsigned int things_to_compute, 
                bool force_recompute,
                int recursion_depth) // TODO: Taken in helper (sub/super) simplex?.
{
    _dosl_verbose_head(1);
    
    // First check if the full simplex exists in collection
    if (!force_recompute) {
        MetricSimplexType dummy_simplex;
        dummy_simplex.p = nps;
        for (int a=0; a<nps.size(); ++a)
            dummy_simplex.gs.push_back (nps[a]->g_score);
        MetricSimplexType* found_p = find (&dummy_simplex); // check if simplex was already created.
        if (found_p) {
            return (found_p);
        }
    }
    
    if (nps.size()==1) {
        MetricSimplexType* ret = create_new_zero_simplex (nps[0]);
        return (ret);
    }
    
    // recursion
    _NodePointerType last_node = nps.back();
    nps.pop_back();
    // recursive call first
    MetricSimplexType* ret = construct_simplex_from_vertices (nps, //distsFromP0, 
                                            (things_to_compute | COMPUTE_INCREMENTAL_QUANTITIES), force_recompute,
                                                recursion_depth+1);
    if (_dosl_verbose_on(0)) {
        _dosl_printf ("In construct_simplex_from_vertices: subsimplex created though recursion: simplex_computation_failure=%s.", 
                        uint_to_binary (ret->simplex_computation_failure) );
    }
    // push back last node
    ret = check_connections_and_add_vertex (ret, last_node, 
                                                (things_to_compute | COMPUTE_INCREMENTAL_QUANTITIES), force_recompute);
    
    if (_dosl_verbose_on(0)) {
        _dosl_printf ("In construct_simplex_from_vertices: simplex after pushing vertex: simplex_computation_failure=%s.", 
                        uint_to_binary (ret->simplex_computation_failure) );
    }
    
    // g_score-score
    if (recursion_depth==0  &&  TO_BOOL(things_to_compute & COMPUTE_G_SCORE_OF_APEX))
        ret->ComputeGScore();
    
    if (_dosl_verbose_on(0)) {
        _dosl_printf ("In construct_simplex_from_vertices: simplex after g_score-score computation: simplex_computation_failure=%s.", 
                        uint_to_binary (ret->simplex_computation_failure) );
    }
    
    return (ret);
    
}

// ========================================================

template <class _NodePointerType, class DoubleType, class DoubleVecType>
std::unordered_set< AMetricSimplex<_NodePointerType,DoubleType,DoubleVecType>* >
    MetricSimplexCollection<_NodePointerType,DoubleType,DoubleVecType>::get_all_attached_maximal_simplices
        (AMetricSimplex<_NodePointerType,DoubleType,DoubleVecType>* in_simplex_p,
            unsigned int things_to_compute, // default: COMPUTE_ALL,
            bool return_base_if_maximal, // default: true
            bool force_compute, // default: false
            bool use_only_expnded_nodes, // default: false
                std::unordered_set<_NodePointerType>* neighbor_node_pointers_p)
{
    _dosl_verbose_head(1);
    
    if (_dosl_verbose_on(0)) {
        in_simplex_p->print("inSimplex");
    }
    
    std::unordered_set<_NodePointerType> allCommonNeighbors, tmpCommonNeighbors;
    bool isAbsoluteMaximal = false;
    bool AllAttachedAbsoluteMaximalsComputed;
    
    if (!neighbor_node_pointers_p) { // Absolute maximal simplices computation
        isAbsoluteMaximal  = true;
        
        // If already computed, return that
        if (in_simplex_p->attached_absolute_maximal_simplices.size()>0  &&  !force_compute)
            return (in_simplex_p->attached_absolute_maximal_simplices);
        
        AllAttachedAbsoluteMaximalsComputed = true;
        
        // Get intersection of the successors.
        // Insert the successors of p[0]
        for (auto it = in_simplex_p->p[0]->successors.begin(); it != in_simplex_p->p[0]->successors.end(); ++it) {
            if (use_only_expnded_nodes && !(it->first->expanded)) {
                AllAttachedAbsoluteMaximalsComputed = false;
                continue;
            }
            allCommonNeighbors.insert (it->first);
        }
        // compute intersection with successors of p[i]
        for (int a=1; a<in_simplex_p->p.size(); ++a) {
            tmpCommonNeighbors = allCommonNeighbors; // so that iterator remains valid
            for (auto it = tmpCommonNeighbors.begin(); it != tmpCommonNeighbors.end(); ++it)
                if (in_simplex_p->p[a]->successors.find(*it) == in_simplex_p->p[a]->successors.end()) // not in intersection
                    allCommonNeighbors.erase (*it);
        }
        
        neighbor_node_pointers_p = &allCommonNeighbors;
    }
    
    if (_dosl_verbose_on(0)) {
        for (auto it=neighbor_node_pointers_p->begin(); it!=neighbor_node_pointers_p->end(); ++it) {
            (*it)->print("Common neighbor node: ");
        }
    }
    
    auto attached_simplices = getAllMaximalSimplicesFromSet (neighbor_node_pointers_p, in_simplex_p, things_to_compute);
    if (_dosl_verbose_on(0)) {
        _dosl_printf(_YELLOW "getAllMaximalSimplicesFromSet returned %d simplices." YELLOW_, attached_simplices.size());
    }
    
    // Remove self if other maximal simplices found
    if (attached_simplices.size() > 1)
        attached_simplices.erase (in_simplex_p);
    
    if (_dosl_verbose_on(0)) {
        for (auto it=attached_simplices.begin(); it!=attached_simplices.end(); ++it) {
            (*it)->print("Attached simplex: ");
        }
    }
    
    return (attached_simplices);
}

//----

template <class _NodePointerType, class DoubleType, class DoubleVecType>
std::unordered_set< AMetricSimplex<_NodePointerType,DoubleType,DoubleVecType>* >
    MetricSimplexCollection<_NodePointerType,DoubleType,DoubleVecType>::getAllMaximalSimplicesFromSet
        (std::unordered_set<_NodePointerType>* _node_set_p, // common neighbors of 'baseSimplex_p'
            AMetricSimplex<_NodePointerType,DoubleType,DoubleVecType>* baseSimplex_p,
                    unsigned int things_to_compute, /*, bool return_base_if_maximal*/ 
                        std::unordered_set<_NodePointerType> restricted_neighbor_set)
{
    /*
    Algorithm:
    Note: Three types of maximal simplices: 
        1. contains v1, not v2, : v1 x getAllMaximalSimplicesFromSet ({v \in node_set | v is connected to v1})
        2. contins v2, not v1, : v2 x getAllMaximalSimplicesFromSet ({v \in node_set | v is connected to v2})
        3. contains neither v1 nor v2 : getAllMaximalSimplicesFromSet (node_set - v1 - v2)
    
    Alternate Algorithm (not implemented):
    i.a. Identify any two vertices, v1 and v2, that are not neighbors.
        i.b. Return {node_set} if failed to identify such a pair of vertices.
    ii. Return: getAllMaximalSimplicesFromSet (node_set - v1)  U  getAllMaximalSimplicesFromSet (node_set - v2)
    */
    _dosl_verbose_head(1);
    
    COPY_IF_NOTNULL_ELSE_CREATE_POINTER_TO_LOCAL 
        (std::unordered_set<_NodePointerType>, _node_set_p, node_set_p, std::unordered_set<_NodePointerType>() );
    
    if (!_node_set_p) { // Request for absolute maximal simplices attached to baseSimplex
        // if it was already computed
        if (baseSimplex_p  &&  !(baseSimplex_p->isGScoreOfCommonNeighborsChanged()) ) { 
                               //^^^ the common neighbor were previosly generated and their g-scores did not change
                return (baseSimplex_p->attached_absolute_maximal_simplices);
        }
        // popuate 'node_set_p' (note: it's an unordered_set)
        for (auto it=baseSimplex_p->all_common_neighbors.begin(); it!=baseSimplex_p->all_common_neighbors.end(); ++it)
            node_set_p->insert (it->first);
    }
    
    
    std::unordered_set <MetricSimplexType*>  ret_simplices, tmp_simplices;
    unsigned int insertion_success;
    
    // Find two vertices in 'node_set' that are not neighbors of each other
    _MS_SMALL_VECTOR<_NodePointerType> v; // v1 = v[0], v2 = v[1]
    
    for (auto it=node_set_p->begin(); it!=node_set_p->end(); ++it) {
        auto next_it = it; next_it++;
        for (auto it2=next_it; it2!=node_set_p->end(); ++it2) // find something in 'node_set' that is not in (*it)->Successors
            if ((*it)->successors.find(*it2)==(*it)->successors.end()) { // not found
                v.push_back (*it);
                v.push_back (*it2);
                break;
            }
        if (v.size()) break;
    }
    
    if (_dosl_verbose_on(0)) {
        baseSimplex_p->print (_RED "baseSimplex:" RED_);
        _dosl_printf (_RED " node_set_p" RED_ " = (");
        for (auto itt=node_set_p->begin(); itt!=node_set_p->end(); ++itt)
            (*itt)->print();
        _dosl_printf (")");
        _dosl_printf ("v = (");
        for (auto itt=v.begin(); itt!=v.end(); ++itt)
            (*itt)->print();
        _dosl_printf (")");
    }
    
    // if not found, this is a complete simplex. 
    // TODO: Make more efficient by searcing in simplex collection first. -- Actually already done by 'check_connections_and_add_vertex'
    if (v.size()==0) {
        MetricSimplexType* ret_simplex = baseSimplex_p; // may be NULL. That's ok.
        for (auto it=node_set_p->begin(); it!=node_set_p->end(); ++it) {
            ret_simplex = check_connections_and_add_vertex (ret_simplex, *it, 
                                                            things_to_compute //(things_to_compute | COMPUTE_INCREMENTAL_QUANTITIES)
                                                            );
        }
        for (auto it2=restricted_neighbor_set.begin(); it2!=restricted_neighbor_set.end(); ++it2)
            if (ret_simplex->checkConnections (*it2))
                return (ret_simplices); // empty
        if (_dosl_verbose_on(0)) {
            _dosl_printf ("ret_simplex=%x: simplex_computation_stage=%s, simplex_computation_failure=%s, things_to_compute=%s", ret_simplex, 
                        _uint_to_binary(ret_simplex->simplex_computation_stage).c_str(),
                        _uint_to_binary(ret_simplex->simplex_computation_failure).c_str(),
                        _uint_to_binary(things_to_compute).c_str());
        }
        
        if (!TO_BOOL(ret_simplex->simplex_computation_failure & things_to_compute))
            ret_simplices.insert (ret_simplex);
        return (ret_simplices);
    }
    
    // Partitions 1 & 2: getAllMaximalSimplicesFromSet ({v \in node_set | v is connected to v[i]}, v[i] x baseSimplex)
    for (int i=0; i<v.size(); ++i) {
        // v[i]->successors  (intersection)  node_set
        std::unordered_set<_NodePointerType> vi_neighbors;
        for (auto it=node_set_p->begin(); it!=node_set_p->end(); ++it)
            if (v[i]->successors.find(*it) != v[i]->successors.end())
                vi_neighbors.insert (*it);
        // v[i] x baseSimplex
        MetricSimplexType* vi_cross_base = check_connections_and_add_vertex (baseSimplex_p, v[i], //NULL, 
                                                                    //(things_to_compute | COMPUTE_INCREMENTAL_QUANTITIES)
                                                                    things_to_compute
                                                                    );
        
        // check
        
        if (TO_BOOL(vi_cross_base->simplex_computation_failure & things_to_compute)) {
            if (_dosl_verbose_on(0)) {
                _dosl_printf("things_to_compute = %s", _uint_to_binary(things_to_compute).c_str());
                vi_cross_base->print("failure in vi_cross_base");
            }
            continue;
        }
        // recursive call
        tmp_simplices = getAllMaximalSimplicesFromSet (&vi_neighbors, vi_cross_base, 
                                                            things_to_compute, restricted_neighbor_set);
        ret_simplices.insert (tmp_simplices.begin(), tmp_simplices.end());
    }
    
    // Partition 3 : getAllMaximalSimplicesFromSet (node_set - {v[0],v[1]}, baseSimplex)
    std::unordered_set<_NodePointerType> remaining_nodes_set = *node_set_p;
    for (int i=0; i<v.size(); ++i) {
        remaining_nodes_set.erase (v[i]);
        restricted_neighbor_set.insert (v[i]);
    }
    
    // recursive call
    tmp_simplices = getAllMaximalSimplicesFromSet (&remaining_nodes_set, baseSimplex_p, 
                                                        //(things_to_compute | COMPUTE_INCREMENTAL_QUANTITIES)
                                                        things_to_compute, 
                                                        restricted_neighbor_set ); // <-- PROBLEM!!
    ret_simplices.insert (tmp_simplices.begin(), tmp_simplices.end());
    
    // --
    
    // Computation of g_score-score.
    if (TO_BOOL(things_to_compute & COMPUTE_G_SCORE_OF_APEX)) {
        auto ret_simplices_copy = ret_simplices; // to keep iterator valid
        for (auto it=ret_simplices_copy.begin(); it!=ret_simplices_copy.end(); ++it) {
            (*it)->ComputeGScore();
            if (TO_BOOL((*it)->simplex_computation_failure & things_to_compute))
                ret_simplices.erase (*it);
        }
    }
    
    // If these are maximl simplices.
    if (!_node_set_p) { // These simplices must be absolute maximal simplices attached to baseSimplex
        for (auto it=ret_simplices.begin(); it!=ret_simplices.end(); ++it) 
            (*it)->is_maximal = true;
        
    }
    
    return (ret_simplices);
}

// ===================================================================
// Path recpnstruction


template <class _NodePointerType, class DoubleType, class MetricSimplexType>
    _DOSL_SMALL_MAP <_NodePointerType, DoubleType>
        computeDistancesFromAllVertices
            (_DOSL_SMALL_MAP <_NodePointerType, DoubleType> point, MetricSimplexType* coord_simplex)
{
    /* Computes the distance of point from the vertices of coord_simplex.
       'point' is expressed as a weighted combination of the vertices.
       Assumes that the coordinates of the vrtices in the 'coord_simplex' are already computed
    */
    typedef typename MetricSimplexType::DoubleVecType DoubleVecType;
    
    // compute the coordinates of the point
    DoubleVecType point_coords (coord_simplex->n_vertices_m1, 0.0);
    typename _DOSL_SMALL_MAP<_NodePointerType,DoubleType>::const_iterator thisWeightIt;
    
    for (int a=0; a < coord_simplex->p.size(); ++a)
        if ((thisWeightIt = point.find (coord_simplex->p[a])) != point.end())
            point_coords = point_coords + (thisWeightIt->second) * coord_simplex->vs[a];
    
    // Compute the distances
    _DOSL_SMALL_MAP <_NodePointerType, DoubleType> ret;
    
    for (int a=0; a < coord_simplex->p.size(); ++a)
        ret[coord_simplex->p[a]] = norm (point_coords - coord_simplex->vs[a]);
    
    return (ret);
}

// --------------------------------------------

template <class _NodePointerType, class DoubleType>
void setAbstractVertexAsSuccessors (_NodePointerType  abstract_vertex_p, 
                                _DOSL_SMALL_MAP<_NodePointerType, DoubleType>  dist_to_other_vertices)
{
    abstract_vertex_p->successors = dist_to_other_vertices;
    for (auto it=dist_to_other_vertices.begin(); it!=dist_to_other_vertices.end(); ++it)
        it->first->successors [abstract_vertex_p] = dist_to_other_vertices [it->first];
}

template <class _NodePointerType, class DoubleType>
void removeAbstractVertexFromSuccessors (_NodePointerType  abstract_vertex_p, 
                                _DOSL_SMALL_MAP<_NodePointerType, DoubleType>  dist_to_other_vertices)
{
    abstract_vertex_p->successors.clear();
    for (auto it=dist_to_other_vertices.begin(); it!=dist_to_other_vertices.end(); ++it)
        it->first->successors.erase (abstract_vertex_p); // [abstract_vertex_p] = dist_to_other_vertices[it->first];
}

// --------------------------------------------

template <class _NodePointerType, class DoubleType, class DoubleVecType>
    PathPoint <_NodePointerType, DoubleType>
        MetricSimplexCollection<_NodePointerType,DoubleType,DoubleVecType>::find_camefrom_point
            (const PathPoint<_NodePointerType, DoubleType>& in_path_point)
{
    _dosl_verbose_head(1);
    
    _DOSL_SMALL_VECTOR<_NodePointerType>            ContainingSimplexVertices;
    _DOSL_SMALL_MAP <_NodePointerType, DoubleType>  inPointFiltered, inPoint = in_path_point.p;
    _NodePointerType inPointAbstractVertexPointer = new_p (_NodePointerType);
    
    // Get rid of zero wight vertices
    for (auto it=inPoint.begin(); it!=inPoint.end(); ++it) 
        if (fabs(it->second) > _MS_DOUBLE_EPS) {
            inPointFiltered [it->first] = it->second;
            ContainingSimplexVertices.push_back (it->first);
        }
    
    // Construct the simplex containing the point
    MetricSimplexType* ContainingSimplex = construct_simplex_from_vertices (ContainingSimplexVertices, //NULL,
                                                                            CHECK_CONNECTION); // won't compute embedding
    if (_dosl_verbose_on(0)) {
        _dosl_printf_nobreak ("Vertex-weight pairs: ");
        for (auto it=inPoint.begin(); it!=inPoint.end(); ++it)
            printf ("%x:%f, ", it->first, it->second);
        printf(". g_score-score at inPoint: %f \n", in_path_point.g_score);
        ContainingSimplex->print("ContainingSimplex: ");
    }
    if (TO_BOOL(ContainingSimplex->simplex_computation_failure)) { // ContainingSimplex is invalid
        if (_dosl_verbose_on(0)) {
            _dosl_printf("simplex_computation_failure=%d,%s", 
                ContainingSimplex->simplex_computation_failure, uint_to_binary(ContainingSimplex->simplex_computation_failure));
        }
        return (PathPointType());
    }
    
    std::unordered_set<MetricSimplexType*> attached_simplices = 
                    get_all_attached_maximal_simplices (ContainingSimplex, //COMPUTE_INCREMENTAL_QUANTITIES 
                                                                        //CHECK_CONNECTION
                                                                        EMBED_SIMPLEX
                                                                ); 
    
    /* For each of the attached maximal simplex,
        i. Compute distances from inPointFiltered to the vertices.
        ii. Extract the maximal subsimplices for which ContainingSimplex is not a subsimplex.
          (start with the simplex with vertices V_{attached maximal simplex} - V_{ContainingSimplex}, and then
           push every other vertex, except one)
       For each of those subsimplices,
        i. Compute distances from inPointFiltered (needs to be computed in the 'attached maximal simplex'(already done earlier).
        ii. Construct a metric simplex with inPointFiltered as p[0] and the subsimplex as the remaining face
            (do compute g_score-score for this).
        iii. Choose the one that gives lowest g_score-score.
    */
    
    DoubleType bestGScore = std::numeric_limits<DoubleType>::max(); 
    MetricSimplexType* bestSimplex = NULL;
    
    if (_dosl_verbose_on(0)) {
        _dosl_printf("'find_camefrom_point': Number of Attached Simplices: %d", attached_simplices.size());
    }
    
    for (auto it=attached_simplices.begin(); it!=attached_simplices.end(); ++it) {
        if (_dosl_verbose_on(0)) {
            (*it)->print("Attached Simplex: ");
        }
        
        if (TO_BOOL((*it)->simplex_computation_failure)) {
            if (_dosl_verbose_on(0)) {
                _dosl_printf("Not processed since it's a failure.");
            }
            continue;
        }
        
        // compute full simplex, *it, 
        
        // Compute distances from inPointFiltered to the vertices
        _DOSL_SMALL_MAP <_NodePointerType, DoubleType>  this_attached_simplex_dists = 
                                                                computeDistancesFromAllVertices (inPointFiltered, *it);
        if (_dosl_verbose_on(0)) {
            _dosl_printf_nobreak ("inPointFiltered wights: ");
            for (auto itt=inPointFiltered.begin(); itt!=inPointFiltered.end(); ++itt)
                printf ("%x (%f), ", itt->first, itt->second);
            printf("\n");
            _dosl_cout << "Distances from vertices of attached simplex: ";
            for (auto itt=this_attached_simplex_dists.begin(); itt!=this_attached_simplex_dists.end(); ++itt)
                printf ("%x (%f), ", itt->first, itt->second);
            printf("\n");
        }
        
        // inPointAbstractVertexPointer->successors = this_attached_simplex_dists;
        setAbstractVertexAsSuccessors (inPointAbstractVertexPointer, this_attached_simplex_dists); // temporary action
        
        // Extract the maximal subsimplices for which ContainingSimplex is not a subsimplex.
        //  These are simplices made up of this_attached_simplex.p - (a vertex from ContainingSimplex)
        for (auto it2=ContainingSimplexVertices.begin(); it2!=ContainingSimplexVertices.end(); ++it2) {
            _DOSL_SMALL_VECTOR<_NodePointerType> VerticesWithPointAsApex = (*it)->p; // the complete attache simplex
            for (int b=0; b<VerticesWithPointAsApex.size(); ++b) { // find the vertex *it2, and replace/remove it
                if (VerticesWithPointAsApex[b] == *it2) {
                    VerticesWithPointAsApex[b] = VerticesWithPointAsApex[0]; // remove the b-th vertex
                    VerticesWithPointAsApex[0] = inPointAbstractVertexPointer;
                }
            }
            
            MetricSimplexType* SimplexWithPointAsApex = construct_simplex_from_vertices 
                                        (VerticesWithPointAsApex, COMPUTE_ALL); // need to compute g_score-score of apex
            if (_dosl_verbose_on(0)) {
                SimplexWithPointAsApex->print ("Potential came-from simplex: ");
                _dosl_printf ("is better [(%f>=) %f > %f]? %d", in_path_point.g_score, 
                                            bestGScore, SimplexWithPointAsApex->g_came_from_point,
                                                (bestGScore > SimplexWithPointAsApex->g_came_from_point));
            }
            
            if (  // !TO_BOOL(SimplexWithPointAsApex->simplex_computation_failure) && 
                    TO_BOOL(SimplexWithPointAsApex->simplex_computation_stage & COMPUTE_G_SCORE_OF_APEX) && 
                      (bestGScore > SimplexWithPointAsApex->g_score) &&
                          (SimplexWithPointAsApex->g_came_from_point < in_path_point.g_score)  ) {
                bestGScore = SimplexWithPointAsApex->g_score;
                bestSimplex = SimplexWithPointAsApex;
            }
        }
        
        removeAbstractVertexFromSuccessors (inPointAbstractVertexPointer, this_attached_simplex_dists);
        
    }
    
    if (bestSimplex)
        return ( PathPointType(bestSimplex->w, bestSimplex->g_came_from_point, in_path_point.p, ContainingSimplex) ); // (bestSimplex->w);
    else
        return (PathPointType());
}

// -------------------------



#endif
