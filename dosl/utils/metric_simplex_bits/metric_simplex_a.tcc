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

#ifndef __METRIC_SIMPLEX_A_TCC_
#define __METRIC_SIMPLEX_A_TCC_

#include <string>
#include <unordered_set>
#include <set>
#include <limits>

#include "metric_simplex_base.tcc"
#include "metric_simplex_collection.hpp"

// ===========================================================

template <class _NodePointerType, class _CostType=double>
class MetricSimplexVertex
{
public:
    typedef _CostType  CostType;
    typedef _NodePointerType  NodePointerType;
    
    CostType g_score;
    bool expanded; // indicates if the g_score-score is usable
    _DOSL_SMALL_MAP<NodePointerType,CostType> successors; // pointer-distance pair
    std::unordered_set<NodePointerType> children_influenced; // 3D
    
    // default constructor
    MetricSimplexVertex()  : g_score (std::numeric_limits<CostType>::max()), expanded(false) { }
    
    // other functions
    
    bool isGScoreValid (bool check_expanded=false) {
        if (!check_expanded)
            return ( g_score != std::numeric_limits<CostType>::max() );
        return ( (g_score != std::numeric_limits<CostType>::max())  &&  expanded);
    }
    
    CostType getDistanceToSuccessor (NodePointerType np) {
        // returns nan if np is not a successor or if the distance to it is nan
        auto found_neighbor_it = successors.find (np);
        if ( found_neighbor_it==successors.end())
            return (std::numeric_limits<CostType>::quiet_NaN());
        return (found_neighbor_it->second); // returns nan if the distance is nan
    }
};

// ------------------------------------------------------------------


template < class _NodePointerType, // '_NodePointerType' should be a pointer type of a class derived from MetricSimplexVertex
                class _DoubleType=double, class _DoubleVecType=_DOSL_SMALL_VECTOR<_DoubleType> >
class AMetricSimplex  // pointed
{
public:
    typedef _NodePointerType NodePointerType;
    typedef _DoubleType DoubleType;
    typedef _DoubleVecType DoubleVecType;
    typedef std::unordered_set <AMetricSimplex*, MetricSimplexHasher<AMetricSimplex*>, 
                                    MetricSimplexEqualTo<AMetricSimplex*> >  MetricSimplexPointersUnorderedSetType;
    typedef MetricSimplexCollection <_NodePointerType,DoubleType,DoubleVecType>  MetricSimplexContainerType;
    
    
    //============================================================================
    // Embedded simplex:
    
    _DOSL_SMALL_VECTOR <_NodePointerType> p; // size n (vertices 0,1,2,...,n-1)
    _DOSL_SMALL_MAP <_NodePointerType, unsigned int> i; // indices.
    
    // ------------------------------
    
    _MS_SMALL_VECTOR <DoubleVecType> vs; // size n x (n-1): vs[0], vs[1], vs[2], ..., vs[n]: Each a (n-1)-vector.
    _MS_SMALL_VECTOR <DoubleVecType> vs_sq; // size n x (n-1)
    
    _MS_SMALL_VECTOR <DoubleVecType> dist_mat; // size n x n
    _MS_SMALL_VECTOR <DoubleVecType> dist_sq_mat; // size n x n
    /*  v[0] is the primary vertex (one that is expanded last in this simplex).
        v[0].p->expanded is false, v[i].p->expanded is true for all other i. */
    
    // ------------------------------
    
    int n_vertices, n_vertices_m1, n_vertices_m2; // n_vertices == p.size()
    
    // -------------------------------------------------------------
    
    void SetDistancesFromLastVertex (const _MS_SMALL_VECTOR<DoubleType>& dists) {
        DoubleType this_dist_sq;
        for (int vertex_num=0; vertex_num < n_vertices_m1; ++vertex_num) {
            dist_mat[n_vertices_m1][vertex_num] = dists[vertex_num];
            dist_mat[vertex_num][n_vertices_m1] = dists[vertex_num];
            this_dist_sq = dists[vertex_num] * dists[vertex_num];
            dist_sq_mat[n_vertices_m1][vertex_num] = this_dist_sq;
            dist_sq_mat[vertex_num][n_vertices_m1] = this_dist_sq;
        }
    }
    
    bool ComputeCoordinateOfLastVertex (void) {
        // Assume enough room, and that 'SetDistancesFromLastVertex' has already been called.
        int j = n_vertices_m1; // Compute coordinates of the j-th vertex. 0-th vertex has all 0 coordinates.
        for (int k=0; k<j-1; ++k) {
            // we are determining vs_sq[j][k-1] // coordsq(k-1,j)
            DoubleType sum = vs_sq[k+1][k]; // coordsq(k-1,k);
            for (int p=0; p<k; ++p)
                sum += vs_sq[k+1][p] - 2*vs[j][p]*vs[k+1][p]; // coordsq(p,k) - 2*coord(p,j)*coord(p,k);
            vs[j][k] = (dist_sq_mat[j][0] - dist_sq_mat[j][k+1] + sum) / (2.0*vs[k+1][k]);
            vs_sq[j][k] = vs[j][k] * vs[j][k];
        }
        // we are determining coordsq(j-1,j)
        vs_sq[j][j-1] = dist_sq_mat[j][0]; //coordsq(j-1,j) = DistMat_sq(0,j);
        for (int p=0; p<j-1; ++p)
            vs_sq[j][j-1] -= vs_sq[j][p]; // coordsq(j-1,j) -= coordsq(p,j);
        
        if (vs_sq[j][j-1] < (DoubleType)(0.0)) // _MS_DOUBLE_EPS degenerate simplex. Don't add!
            return (false);
        vs[j][j-1] = sqrt(vs_sq[j][j-1]); //coord(j-1,j) = sqrt(coordsq(j-1,j));
        return (true);
    }
    
    //============================================================================
    
    // ------------------------------
    int back_track_simplex;
    AMetricSimplex* apex_came_from_face;
    
    // ------------------------------
    // Maximality
    bool is_maximal; // is self maximal
    // has all the absolute maximal simplices attached computed?
    std::unordered_set<AMetricSimplex*> attached_absolute_maximal_simplices;
    // all common neighborss
    std::unordered_map<_NodePointerType,DoubleType> all_common_neighbors; // vertex and g_score-score pairs
    /* Note: If 'attached_absolute_maximal_simplices' were previously computed and the 
             'all_common_neighbors' is exactly the same, we can return the attached_absolute_maximal_simplices. */
    
    // ------------------------------
    
    DoubleVecType o; // origin/start. size n-1
    DoubleVecType gs; // size n-1 (g-score of vertices 1,2,...,n-1). Use gs[1,2,...,n-1], ignore gs[0]
    DoubleVecType gsSq;
    DoubleVecType A, B; // will use indices 1,2,...,n-2
    DoubleType C2, C1, C0;
    
    DoubleType g_score; // of p[0]
    // NOTE: p[0]->g_score is the best g_score-score of the vertex across all simplices, while g_score-score of simplex (member g_score) is the g_score-score for path through this simplex.
    
    SetOfPairsOrderedBySecond <int, DoubleType> tw; // index-weight pairs (temporary variable)
    _DOSL_SMALL_MAP <_NodePointerType, DoubleType> w; // final computed weights
    DoubleType g_came_from_point;
    
    _MS_SMALL_VECTOR<AMetricSimplex*> face_simplices; // use indices 1,2,...,n-1. storing for fast computation of w  // non-incremental
    
    // ---------------------
    // computation tracking and error recording
    unsigned int simplex_computation_stage; // stage up to which computation has been completed successfully.
    unsigned int simplex_computation_failure; // stage at which a failure occurred. 0u if no failure yet.
    
    MetricSimplexContainerType*  all_simplices_p;
    
    // ---------------------------------
    
    void print_pgs (std::string head="", std::string tail="") {
        _dosl_cout << _RED + head << " (" << this << ")" RED_ ":" << " n_vertices=" << n_vertices << "; " << _dosl_endl;
        _dosl_cout << "p[0]: " << p[0] << ", ";
        for (int a=1; a<n_vertices; ++a) {
            std::cout << "p[" << a << "]: (" << p[a] << "," << gs[a] << "), ";
        }
        _dosl_cout << tail;
    }
    
    void print (std::string head="", std::string tail="") {
        _dosl_cout << _RED + head << " (" << this << ")" RED_ ":" << " n_vertices=" << n_vertices << _dosl_endl;
        for (int a=0; a<n_vertices; ++a) {
                p[a]->print ("vertex "+std::to_string(a)+((a==0)?" (owner)":""));
                _dosl_cout << "Embedded coordinate = (";
                for (int b=0; b<vs[a].size(); ++b)
                    std::cout << vs[a][b] << ", ";
                std::cout << "); g-score in simplex = " << ((a==0)?g_score:gs[a]) << _dosl_endl;
                if ((this==p[0]->CameFromSimplex)  &&  fabs(((a==0)?g_score:gs[a]) - p[a]->g_score) > _MS_DOUBLE_EPS)
                    std::cout << " (g-score in came-from simplex is different!!! back_track_simplex=" << back_track_simplex << ")" << _dosl_endl;
        }
        
        _dosl_printf(_GREEN "Other simplex data:" GREEN_);
        _dosl_cout << "o: " << "Embedded coordinate = (";
        for (int b=0; b<o.size(); ++b)
            std::cout << o[b] << ", ";
        std::cout << ")" << _dosl_endl;
        
        _dosl_printf("simplex_computation_stage=%s, simplex_computation_failure=%s", 
                uint_to_binary(simplex_computation_stage), uint_to_binary(simplex_computation_failure) );
        if (TO_BOOL(simplex_computation_stage & COMPUTE_G_SCORE_OF_APEX)) {
            _dosl_printf_nobreak ("weights = [");
            for (int a=1; a<p.size(); ++a)
                printf("(%d,%x): %f, ", a, p[a], w[p[a]]);
            printf("]\n");
            _dosl_cout << "g_came_from_point = ";
            (g_came_from_point==std::numeric_limits<double>::max())? printf("INF") : printf("%f", g_came_from_point);
            _dosl_cout << _dosl_endl;
        }
        
        std::cout << tail;
        std::cout.flush();
    }
    
    // ===========================================================================
    
    bool isValid (void) { // checks if the g_score-scores stored at the base vertices is same as 'gs'
        for (int a=1; a<gs.size(); ++a)
            if (fabs(gs[a] - p[a]->g_score) > _MS_DOUBLE_EPS)
                return (false);
        return (true);
    }
    
    // ------------------------------------------------------------------
    
    void MakeRoomForAnotherVertex (_NodePointerType newVertexPtr=NULL) {
        p.push_back (newVertexPtr);
        for (int vertex_num=0; vertex_num<n_vertices; ++vertex_num) {
            vs[vertex_num].push_back (0.0); // now the length of each vertex is 'n_vertices'
            vs_sq[vertex_num].push_back (0.0);
            dist_mat[vertex_num].push_back (0.0); // Now each column is 'n_vertices + 1' long
            dist_sq_mat[vertex_num].push_back (0.0);
        }
        vs.push_back (DoubleVecType(n_vertices,0.0));
        vs_sq.push_back (DoubleVecType(n_vertices,0.0));
        A.resize (n_vertices);
        B.resize (n_vertices);
        n_vertices_m2 = n_vertices_m1;
        n_vertices_m1 = n_vertices;
        ++n_vertices;
        dist_mat.push_back (DoubleVecType(n_vertices,0.0));
        dist_sq_mat.push_back (DoubleVecType(n_vertices,0.0));
        
        o.resize (n_vertices-1); //push_back (0.0);
        gs.resize (n_vertices); // ignore index 0 //push_back (0.0);
        gsSq.resize (n_vertices); // ignore index 0 //push_back (0.0);
        
        face_simplices.resize (n_vertices, NULL);
        //w.resize (n_vertices, 0.0);
    }
    
    // -------------------------------------------------------------
    
    void SetGscoreOfLastVertex (const DoubleType& gg) {
        gs[n_vertices_m1] = gg;
        gsSq[n_vertices_m1] = gg*gg;
    }
    
    // --------------------------------
    
    bool ComputeCoordinateOfO (void) {
        // need to compute A[n_vertices-2], B[n_vertices-2] and o
        // Both A and B is of size 'n_vertices_m1' with indices 0,1,...,n_vertices-2. Discard index 0.
        if (n_vertices > 2) {
            int k = n_vertices_m2; // >= 1
            int kp1 = n_vertices_m1;
            A[k] = vs[1][0] - vs[kp1][0]; // \alpha_{k,0} * vs[k+1][k]
            B[k] = gsSq[1] - gsSq[kp1] -vs_sq[1][0]; // \beta_k * vs[k+1][k]
            for (int p=0; p<kp1; ++p)  B[k] += vs_sq[kp1][p]; // p=0,1,...,k
            B[k] /= 2;
            for (int p=1; p<k; ++p) { // p=1,...,k-1
                A[k] -= vs[kp1][p] * A[p];
                B[k] -= vs[kp1][p] * B[p];
            }
            A[k] /= vs[kp1][k]; // already checked that vs[kp1][k] is non-zero
            B[k] /= vs[kp1][k];
        }
        C2 += A[n_vertices_m2] * A[n_vertices_m2];
        C1 += 2 * A[n_vertices_m2] * B[n_vertices_m2];
        C0 += B[n_vertices_m2] * B[n_vertices_m2];
        
        DoubleType C1_2=C1/C2, C0_2=C0/C2;
        DoubleType CR = C1_2*C1_2 - 4*C0_2;
        //printf ("C0=%f, C1=%f, C2=%f, CR=%f\n", C0, C1, C2, CR);
        
        //printf("\nCR = %e\n", CR);
        
        if (CR < (DoubleType)(-_MS_DOUBLE_EPS_SQ)) // _MS_DOUBLE_EPS degenerate simplex
        //if (CR < (DoubleType)(0.0)) // degenerate simplex
            return (false);
        else if (CR < (DoubleType)(_MS_DOUBLE_EPS_SQ))
            CR = (DoubleType)0.0;
        o[0] = (sqrt(CR) - C1_2) / 2; //(-C1 + sqrt(C1*C1 - 4.0*C2*C0)) / (2.0*C2);
        for (int p=1; p<n_vertices_m1; ++p) // 1,2,...,n-2
            o[p] = A[p]*o[0] + B[p];
        
        return (true);
    }
    
    // ==============================================================
    
    // constructors
    AMetricSimplex (MetricSimplexContainerType* sp=NULL) : 
            back_track_simplex(0),
            n_vertices(0), n_vertices_m1(-1), n_vertices_m2(-2), g_score(std::numeric_limits<DoubleType>::max()),
            gs(DoubleVecType(1,-1.0)), all_simplices_p(NULL),
            simplex_computation_stage(0u), simplex_computation_failure(0u), apex_came_from_face(NULL),
            is_maximal (false), g_came_from_point(std::numeric_limits<DoubleType>::max()) { }
    
    AMetricSimplex (_NodePointerType np, MetricSimplexContainerType* sp=NULL) 
            : back_track_simplex(0),
              n_vertices(0), n_vertices_m1(-1), n_vertices_m2(-2), g_score(std::numeric_limits<DoubleType>::max()),
              gs(DoubleVecType(1,-1.0)), all_simplices_p(sp), 
              simplex_computation_stage(COMPUTE_ALL), simplex_computation_failure(0u), // everything computed
              apex_came_from_face(NULL), is_maximal (false), g_came_from_point(std::numeric_limits<DoubleType>::max())
    // np = p[0]
    {
        MakeRoomForAnotherVertex (np);
        // vs, vs_sq are empty; dist_mat and dist_sq_mat are 1x1 containg 0; w has one element with 0; -- no need to set these
    }
    
    AMetricSimplex (_NodePointerType just_created, _NodePointerType came_from, 
                        DoubleType d=std::numeric_limits<DoubleType>::quiet_NaN(), MetricSimplexContainerType* sp=NULL) 
            : back_track_simplex(0),
              n_vertices(0), n_vertices_m1(-1), n_vertices_m2(-2), g_score(std::numeric_limits<DoubleType>::max()),
              gs(DoubleVecType(1,-1.0)), all_simplices_p(sp), 
              simplex_computation_stage(0u), simplex_computation_failure(0u), apex_came_from_face(NULL),
              is_maximal (false), g_came_from_point(std::numeric_limits<DoubleType>::max())
    // just_created = p[1], came_from = p[0]
    {
        MakeRoomForAnotherVertex (just_created);
        // vs, vs_sq are empty; dist_mat and dist_sq_mat are 1x1 containg 0; w has one element with 0; -- no need to set these
        
        MakeRoomForAnotherVertex (came_from);
        
        // Phase 1: Compute embedding of vertices
        // ----
        
        if (std::isnan(d)) {
            if ( std::isnan(d = came_from->getDistanceToSuccessor (just_created)) ) {
                simplex_computation_failure = CHECK_CONNECTION;
                return;
            }
        }
        
        SetDistancesFromLastVertex (_MS_SMALL_VECTOR<DoubleType>(1,d));
        simplex_computation_stage |= CHECK_CONNECTION;
        
        // --
        if (ComputeCoordinateOfLastVertex()) 
            simplex_computation_stage |= COMPUTE_LOCAL_COORDINATE;
        else {
            simplex_computation_failure = COMPUTE_LOCAL_COORDINATE;
            return;
        }
        
        // Phase 2: Compute g_score-score
        // ----
        
        SetGscoreOfLastVertex (came_from->g_score);
        simplex_computation_stage |= SET_G_SCORE_OF_VERTEX;
        
        // --
        C2 = 1.0;
        C1 = -2*vs[1][0];
        C0 = vs[1][0]*vs[1][0] - gsSq[1];
        
        if (ComputeCoordinateOfO ()) 
            simplex_computation_stage |= COMPUTE_COORDINATE_OF_O;
        else {
            simplex_computation_failure = COMPUTE_COORDINATE_OF_O;
            return;
        }
        
        
        g_score = came_from->g_score + d;
        w[p[1]] = 1.0;
        simplex_computation_stage |= COMPUTE_G_SCORE_OF_APEX; // no failure possible here
        
        apex_came_from_face = this;
    }
    
    // ==============================================================
    
    bool checkConnections (_NodePointerType np, 
                                _MS_SMALL_VECTOR<DoubleType>* dists_p=NULL, // return distances if not NULL
                                // DoubleType const* dist_to_p0=NULL,
                                bool stopIfConnectionFails=false)
    {
        _dosl_verbose_head(1);
        
        // sets NaN for whicherver connection does not exist.
        // 'stopIfConnectionFails' is redundent if 'dists_p' is NULL
        bool allConnectionsExist = true;
        
        // Assume 'np->successors' have been generated.
        for (int a=0; a<n_vertices; ++a) {
            // skip checking with p[0]
            /* if (a==0 && dist_to_p0) {
                if (dists_p)
                    dists_p->push_back (*dist_to_p0);
                continue;
            } */
            
            auto found_it = np->successors.find (p[a]);
            
            if (found_it==np->successors.end() || std::isnan(found_it->second)) { // np cannot be inserted in this simplex
                if (_dosl_verbose_on(0)) {
                    np->print("np = ");
                    p[a]->print("p[a] = ");
                    //printf ("(found_it==np->successors.end())=%d, found_it->second=%f", (found_it==np->successors.end()), found_it->second);
                }
                if (!stopIfConnectionFails  &&  dists_p) {
                    dists_p->push_back ( std::numeric_limits<DoubleType>::quiet_NaN() );
                    allConnectionsExist = false;
                }
                else // if 'stopIfConnectionFails' is true or 'dists_p' is NULL
                    return (false);
            }
            else if (dists_p)
                dists_p->push_back (found_it->second);
        }
        
        return (allConnectionsExist);
    }
    
    // -------------------------------------------
    
    
    bool isGScoreOfCommonNeighborsChanged (bool forceRecompute=false)
    {
        // Returns true if updated with the current g_score-scores, 
        //      otherwise (if returning previously computed list) return false
        
        if (all_common_neighbors.size()) {
            // check if g_score-scores changed
            bool GScoresUpdated = false;
            for (auto it=all_common_neighbors.begin(); it!=all_common_neighbors.end();++it)
                if (it->first->g_score != it->second) {
                    if (forceRecompute) {
                        GScoresUpdated = true;
                        it->second = it->first->g_score; // update 'all_common_neighbors' if it changed.
                    }
                    else
                        return (true);
                    /* NOTE: g_score-score of a vertex is the best g_score-score across all simplices, 
                        while g_score-score of simplex is the g_score-score of apex for a path through that simplex. */
                }
            // check if the g_score-scores are to be updated.
            return (GScoresUpdated);
        }

        // not previously computed. recompute from scratch.
        all_common_neighbors.clear();
        std::unordered_map<_NodePointerType,DoubleType>  tmpCommonNeighbors;
        // The successors of p[0] that are not part of inSimplex
        for (auto it=p[0]->successors.begin(); it!=p[0]->successors.end(); ++it)
            all_common_neighbors[it->first] = it->second; // Note: successors is an unordered map
        for (int a=1; a<p.size(); ++a) {
            // remove p[a] if it existed in allCommonNeighbors
            all_common_neighbors.erase (p[a]);
            // now remove everything else that's not part of p[a].Successor
            tmpCommonNeighbors = all_common_neighbors; // so that iterator remains valid
            for (auto it = tmpCommonNeighbors.begin(); it != tmpCommonNeighbors.end(); ++it)
                if (p[a]->successors.find(it->first) == p[a]->successors.end())
                    all_common_neighbors.erase (it->first);
        }
        
        return (true);
    }
    
    // ----------------------------------------------------------------------------------------
    // computing shortest distance between p[0] and o, and corresponding weights, within the 
    
    AMetricSimplex* GetFaceSimplex (int i, SetOfPairsOrderedBySecond<int, DoubleType>* FaceComputationOrder) {
        // compute face simplex coordinates and weights (recursive call to 'ComputeGScore')
        // FaceComputationOrder will contain i as well
        // 'FaceComputationOrder' is used to prioritize the order in which the sub-sub-simplices are created.
        _dosl_verbose_head(1);
        
        int a, b; // temporary variables
        if (_dosl_verbose_on(0)) {
            _dosl_printf ("GetFaceSimplex(%d) for simplex %x: ", i, this); // std::cout.flush();
        }
        
        // Test face_simplices[i]
        if (face_simplices[i]) {
            // No need to check if g-scores changed, because if they did, this would be a different simplex.
            if (_dosl_verbose_on(0)) {
                _dosl_printf ("Already computed. Returning %x.", face_simplices[i]);
            }
            return (face_simplices[i]);
        }
        
        // --------------
        // Search in 'all_metric_simplex_pointers'
        
        AMetricSimplex* found_simplex = NULL;
        AMetricSimplex* search_simplex = new AMetricSimplex; // temporary metric simplex for comparison
        search_simplex->p = p;
        search_simplex->gs = gs;
        search_simplex->n_vertices = p.size();
        
        if (_dosl_verbose_on(0)) {
            _dosl_printf ("Will do cascading search in AllSimplices table (size=%d). ", all_simplices_p->size());
        }
        
        int ind, ind_in_face;
        _MS_SMALL_VECTOR<int> removed_vertex_indices;
        // search for the starting elementary simplex (note: VertexInsertionOrderReversed.size() == n_vertices-1 )
        typename SetOfPairsOrderedBySecond<int,DoubleType>::iterator ind_it;
        if (FaceComputationOrder)
            ind_it = FaceComputationOrder->begin();
        for (a=0; a<n_vertices; ++a) {
            // Remove:  i, {FaceComputationOrder[0], FaceComputationOrder[1], ..., FaceComputationOrder[n-2]} - i
            if (a==0)
                ind = i;
            else {
                if (FaceComputationOrder) {
                    ind = ind_it->first; //VertexInsertionOrderReversed [removed_vertex_count];
                    ++ind_it;
                }
                else
                    ind = a;
                if (ind==i) continue;
            }
            // remove p[ind] from 'search_simplex'
            
            removed_vertex_indices.push_back (ind);
            if (a == n_vertices-1) // no point in searching for a 0-simplex!
                break;
            
            ind_in_face = search_simplex->p.findi (p[ind]);
            search_simplex->p.erase (ind_in_face);
            search_simplex->gs.erase (ind_in_face);
            --search_simplex->n_vertices;
            
            if (_dosl_verbose_on(0)) {
                search_simplex->print_pgs ("search_simplex: ");
            }
            found_simplex = all_simplices_p->find (search_simplex);
            if (found_simplex) break;
        }
        
        delete search_simplex; // TODO: can use this for efficiency?
        
        if (ind==i && found_simplex) { // found the exact required face simplex
            face_simplices[i] = found_simplex;
            return (face_simplices[i]);
        }
        
        // Build the face simplex (insert one vertex at a time)
        
        AMetricSimplex* face_simplex;
        /* int ind, ind2, b;
        DoubleVecType dists; */
        
        if (found_simplex) {
            if (_dosl_verbose_on(0)) {
                _dosl_printf ("Found in AllSimplices table.\n");
                //found_simplex->print ("", 1);
            }
            face_simplex = found_simplex; // removed_vertex_count < VertexInsertionOrderReversed.size().
                                          // VertexInsertionOrderReversed [0, 1, 2, ..., removed_vertex_count-1] removed
        }
        else {
            if (_dosl_verbose_on(0)) {
                _dosl_printf ("Not found! Will construct new.\n");
            }
            ind = removed_vertex_indices.back();
            removed_vertex_indices.pop_back();
            face_simplex = all_simplices_p->create_new_one_simplex (p[0], p[ind], dist_mat[0][ind]);
            if (_dosl_verbose_on(0)) {
                face_simplex->print_pgs ("new 1-simplex created");
            }
        }
        
        for (a=removed_vertex_indices.size()-1; a>0; --a) { // removed_vertex_indices[0] == i. So don't insert that.
            ind = removed_vertex_indices [a]; // p[ind] to be inserted
            
            face_simplex = all_simplices_p->check_connections_and_add_vertex (face_simplex, p[ind]);
            if (_dosl_verbose_on(0)) {
                face_simplex->print_pgs ("new simplex created");
            }
            
            if (!face_simplex) {
                // face_simplex is NULL -- simplex with imaginary coordinates or degenerate
                return (NULL);
            }
        }
        
        return (face_simplex);
    }
    
    // =========================================================
    
    bool ComputeGScore (int depth=0) {
        _dosl_verbose_head(1);
        
        if (TO_BOOL(simplex_computation_stage & COMPUTE_G_SCORE_OF_APEX)) // already computed.
            return (true);
        
        DoubleVecType tw0 (n_vertices); // ignore index 0. use 1,2,...,n_vertices-1
        
        // Compute naieve weights
        DoubleType wvSum, wSum=0.0;
        for (int j=n_vertices-1; j>0; --j) {
            wvSum = 0.0;
            for (int i=j+1; i<n_vertices; ++i)
                wvSum += tw0[i] * vs[i][j-1];
            tw0[j] = (o[j-1] - wvSum) / vs[j][j-1];
            wSum += tw0[j];
        }
        // ordered set.
        
        DoubleType twt;
        for (int j=1; j<n_vertices; ++j) {
            twt = tw0[j]/wSum;
            tw.insert ( std::make_pair(j,twt) );
            w[p[j]] = twt;
        }
        DoubleType wo = 1.0/wSum; // g_score of cameFromPoint = g_score * (1 - wo);
        
        if (_dosl_verbose_on(0)) {
            _dosl_printf_nobreak ("call to ComputeGScore (this=%x, depth=%d): n_vertices = %d, naive weights = [", 
                                                                    this, depth, n_vertices);
            for (auto it=tw.begin(); it!=tw.end(); ++it)
                printf(" (%x,%d,%f) ", p[it->first], it->first, it->second);
            printf("]\n");
        }
        // tw contains the weights in increasing order. tw.begin()->second is the smallest weight.
        
        if (tw.begin()->second >= 0.0) { // all weights are positive
            g_score = distance_between_vertices (o, vs[0]);
            g_came_from_point = (1.0 - wo) * g_score;
            apex_came_from_face = this;
            simplex_computation_stage |= COMPUTE_G_SCORE_OF_APEX;
            return (true);
        }
        
        //--
        // create subsimplex by removing a negative weight vertex.
        AMetricSimplex* face_simplex;
        g_score = std::numeric_limits<DoubleType>::max();
        
        // int best_face_search_depth = max_depth_to_search-1, this_face_search_depth;
        
        for (auto it=tw.begin(); it!=tw.end(); ++it) {
            if (it->second < 0.0) { // need to check all faces that have negative weights
                
                if (face_simplices[it->first]) {
                    if ( !TO_BOOL(face_simplices[it->first]->simplex_computation_stage & COMPUTE_G_SCORE_OF_APEX)  && 
                           !TO_BOOL(face_simplices[it->first]->simplex_computation_failure) )
                       face_simplices[it->first]->ComputeGScore (/*distsFromP0,*/ depth+1);
                   face_simplex = face_simplices[it->first];
                }
                else {
                    // Polulate the list of vertices that make up this face
                    _DOSL_SMALL_VECTOR <_NodePointerType> pFace = p;
                    pFace.erase (it->first);
                    face_simplex = all_simplices_p->construct_simplex_from_vertices (pFace, /*distsFromP0,*/ COMPUTE_ALL); 
                                                                                        // compute g_score-score as well
                    face_simplices[it->first] = face_simplex;
                }
                
                // Compare
                if (face_simplex  &&  TO_BOOL(face_simplex->simplex_computation_stage & COMPUTE_G_SCORE_OF_APEX)  &&
                        face_simplex->g_score < g_score ) {
                    g_score = face_simplex->g_score;
                    g_came_from_point = face_simplex->g_came_from_point;
                    // Update weight assignment.
                    w[p[it->first]] = 0.0;
                    for (auto itW=face_simplex->w.begin(); itW!=face_simplex->w.end(); ++itW) // copy weights
                        w[itW->first] = itW->second;
                    apex_came_from_face = face_simplex->apex_came_from_face;
                    // simplex_computation_stage >= COMPUTE_G_SCORE_OF_APEX;
                    simplex_computation_stage |= COMPUTE_G_SCORE_OF_APEX;
                }
                
            }
            else
                break;
        }
        
        return (TO_BOOL(simplex_computation_stage & COMPUTE_G_SCORE_OF_APEX));
    }
    
    // ---------------------------
    
    void CompleteComputation (unsigned int things_to_compute) {
        // TODO: Other steps?
        if ( (simplex_computation_stage & EMBED_ALL) && !(simplex_computation_stage & COMPUTE_G_SCORE_OF_APEX) 
                && (things_to_compute & COMPUTE_G_SCORE_OF_APEX) )
            ComputeGScore();
    }
    
    // ==============================================================
    
    void SetChildInfluence (void) { // to be called on new camefrom simple, when updating g_score-score
        for (int a=1; a<p.size(); ++a)
            p[a]->children_influenced.insert (p[0]);
    }
    
    void UnsetChildInfluence (void) { // to be called on previous camefrom simple, when updating g_score-score
        for (int a=1; a<p.size(); ++a)
            p[a]->children_influenced.erase (p[0]);
    }
};

#endif
