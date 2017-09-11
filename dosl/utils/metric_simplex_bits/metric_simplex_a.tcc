/** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal Search Library (DOSL)                                              *
*    A template-based C++ library for discrete search                                    *
*    Version 3.1                                                                         *
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

template <class nodePointerType, class costType=double>
class MetricSimplexVertex
{
public:
    typedef costType  CostType;
    typedef nodePointerType  NodePointerType;
    
    CostType G;
    bool Expanded; // indicates if the G-score is usable
    _DOSL_SMALL_MAP<NodePointerType,CostType> Successors; // pointer-distance pair
    std::unordered_set<NodePointerType> ChildrenInfluenced; // 3D
    
    // default constructor
    MetricSimplexVertex()  : G (std::numeric_limits<costType>::max()), Expanded(false) { }
    
    // other functions
    
    bool isGScoreValid (bool check_expanded=false) {
        if (!check_expanded)
            return ( G != std::numeric_limits<costType>::max() );
        return ( (G != std::numeric_limits<costType>::max())  &&  Expanded);
    }
    
    CostType getDistanceToSuccessor (NodePointerType np) {
        // returns nan if np is not a successor or if the distance to it is nan
        auto found_neighbor_it = Successors.find (np);
        if ( found_neighbor_it==Successors.end())
            return (std::numeric_limits<CostType>::quiet_NaN());
        return (found_neighbor_it->second); // returns nan if the distance is nan
    }
};

// ------------------------------------------------------------------


template < class nodePointerType, // 'nodePointerType' should be a pointer type of a class derived from MetricSimplexVertex
                class doubleType=double, class doubleVecType=_DOSL_SMALL_VECTOR<doubleType> >
class AMetricSimplex  // pointed
{
public:
    typedef nodePointerType NodePointerType;
    typedef std::unordered_set <AMetricSimplex*, MetricSimplexHasher<AMetricSimplex*>, 
                                    MetricSimplexEqualTo<AMetricSimplex*> >  MetricSimplexPointersUnorderedSetType;
    typedef MetricSimplexCollection <nodePointerType,doubleType,doubleVecType>  MetricSimplexContainerType;
    typedef doubleVecType DoubleVecType;
    
    //============================================================================
    // Embedded simplex:
    
    _DOSL_SMALL_VECTOR <nodePointerType> p; // size n (vertices 0,1,2,...,n-1)
    _DOSL_SMALL_MAP <nodePointerType, unsigned int> i; // indices.
    
    // ------------------------------
    
    _MS_SMALL_VECTOR <doubleVecType> vs; // size n x (n-1): vs[0], vs[1], vs[2], ..., vs[n]: Each a (n-1)-vector.
    _MS_SMALL_VECTOR <doubleVecType> vsSq; // size n x (n-1)
    
    _MS_SMALL_VECTOR <doubleVecType> distMat; // size n x n
    _MS_SMALL_VECTOR <doubleVecType> distSqMat; // size n x n
    /*  v[0] is the primary vertex (one that is expanded last in this simplex).
        v[0].p->Expanded is false, v[i].p->Expanded is true for all other i. */
    
    // ------------------------------
    
    int n_vertices, n_vertices_m1, n_vertices_m2; // n_vertices == p.size()
    
    // -------------------------------------------------------------
    
    void SetDistancesFromLastVertex (const _MS_SMALL_VECTOR<doubleType>& dists) {
        doubleType thisDistSq;
        for (int vertex_num=0; vertex_num < n_vertices_m1; ++vertex_num) {
            distMat[n_vertices_m1][vertex_num] = dists[vertex_num];
            distMat[vertex_num][n_vertices_m1] = dists[vertex_num];
            thisDistSq = dists[vertex_num] * dists[vertex_num];
            distSqMat[n_vertices_m1][vertex_num] = thisDistSq;
            distSqMat[vertex_num][n_vertices_m1] = thisDistSq;
        }
    }
    
    bool ComputeCoordinateOfLastVertex (void) {
        // Assume enough room, and that 'SetDistancesFromLastVertex' has already been called.
        int j = n_vertices_m1; // Compute coordinates of the j-th vertex. 0-th vertex has all 0 coordinates.
        for (int k=0; k<j-1; ++k) {
            // we are determining vsSq[j][k-1] // coordsq(k-1,j)
            doubleType sum = vsSq[k+1][k]; // coordsq(k-1,k);
            for (int p=0; p<k; ++p)
                sum += vsSq[k+1][p] - 2*vs[j][p]*vs[k+1][p]; // coordsq(p,k) - 2*coord(p,j)*coord(p,k);
            vs[j][k] = (distSqMat[j][0] - distSqMat[j][k+1] + sum) / (2.0*vs[k+1][k]);
            vsSq[j][k] = vs[j][k] * vs[j][k];
        }
        // we are determining coordsq(j-1,j)
        vsSq[j][j-1] = distSqMat[j][0]; //coordsq(j-1,j) = DistMat_sq(0,j);
        for (int p=0; p<j-1; ++p)
            vsSq[j][j-1] -= vsSq[j][p]; // coordsq(j-1,j) -= coordsq(p,j);
        
        if (vsSq[j][j-1] < (doubleType)(0.0)) // _MS_DOUBLE_EPS degenerate simplex. Don't add!
            return (false);
        vs[j][j-1] = sqrt(vsSq[j][j-1]); //coord(j-1,j) = sqrt(coordsq(j-1,j));
        return (true);
    }
    
    //============================================================================
    
    // ------------------------------
    int backTrackSimplex;
    AMetricSimplex* ApexCameFromFace;
    
    // ------------------------------
    // Maximality
    bool isMaximal; // is self maximal
    // has all the absolute maximal simplices attached computed?
    std::unordered_set<AMetricSimplex*> AttachedAbsoluteMaximalSimplices;
    // all common neighborss
    std::unordered_map<nodePointerType,doubleType> AllCommonNeighbors; // vertex and G-score pairs
    /* Note: If 'AttachedAbsoluteMaximalSimplices' were previously computed and the 
             'AllCommonNeighbors' is exactly the same, we can return the AttachedAbsoluteMaximalSimplices. */
    
    // ------------------------------
    
    doubleVecType o; // origin/start. size n-1
    doubleVecType gs; // size n-1 (g-score of vertices 1,2,...,n-1). Use gs[1,2,...,n-1], ignore gs[0]
    doubleVecType gsSq;
    doubleVecType A, B; // will use indices 1,2,...,n-2
    doubleType C2, C1, C0;
    
    doubleType G; // of p[0]
    // NOTE: p[0]->G is the best G-score of the vertex across all simplices, while G-score of simplex (member G) is the G-score for path through this simplex.
    
    set_of_pairs_ordered_by_second <int, doubleType> tw; // index-weight pairs (temporary variable)
    _DOSL_SMALL_MAP <nodePointerType, doubleType> w; // final computed weights
    doubleType G_cameFromPoint;
    
    _MS_SMALL_VECTOR<AMetricSimplex*> FaceSimplices; // use indices 1,2,...,n-1. storing for fast computation of w  // non-incremental
    
    // ---------------------
    // computation tracking and error recording
    unsigned int simplex_computation_stage; // stage up to which computation has been completed successfully.
    unsigned int simplex_computation_failure; // stage at which a failure occurred. 0u if no failure yet.
    
    MetricSimplexContainerType*  AllSimplices_p;
    
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
                std::cout << "); g-score in simplex = " << ((a==0)?G:gs[a]) << _dosl_endl;
                if ((this==p[0]->CameFromSimplex)  &&  fabs(((a==0)?G:gs[a]) - p[a]->G) > _MS_DOUBLE_EPS)
                    std::cout << " (g-score in came-from simplex is different!!! backTrackSimplex=" << backTrackSimplex << ")" << _dosl_endl;
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
            _dosl_cout << "G_cameFromPoint = ";
            (G_cameFromPoint==std::numeric_limits<double>::max())? printf("INF") : printf("%f", G_cameFromPoint);
            _dosl_cout << _dosl_endl;
        }
        
        std::cout << tail;
        std::cout.flush();
    }
    
    // ===========================================================================
    
    bool isValid (void) { // checks if the G-scores stored at the base vertices is same as 'gs'
        for (int a=1; a<gs.size(); ++a)
            if (fabs(gs[a] - p[a]->G) > _MS_DOUBLE_EPS)
                return (false);
        return (true);
    }
    
    // ------------------------------------------------------------------
    
    void MakeRoomForAnotherVertex (nodePointerType newVertexPtr=NULL) {
        p.push_back (newVertexPtr);
        for (int vertex_num=0; vertex_num<n_vertices; ++vertex_num) {
            vs[vertex_num].push_back (0.0); // now the length of each vertex is 'n_vertices'
            vsSq[vertex_num].push_back (0.0);
            distMat[vertex_num].push_back (0.0); // Now each column is 'n_vertices + 1' long
            distSqMat[vertex_num].push_back (0.0);
        }
        vs.push_back (doubleVecType(n_vertices,0.0));
        vsSq.push_back (doubleVecType(n_vertices,0.0));
        A.resize (n_vertices);
        B.resize (n_vertices);
        n_vertices_m2 = n_vertices_m1;
        n_vertices_m1 = n_vertices;
        ++n_vertices;
        distMat.push_back (doubleVecType(n_vertices,0.0));
        distSqMat.push_back (doubleVecType(n_vertices,0.0));
        
        o.resize (n_vertices-1); //push_back (0.0);
        gs.resize (n_vertices); // ignore index 0 //push_back (0.0);
        gsSq.resize (n_vertices); // ignore index 0 //push_back (0.0);
        
        FaceSimplices.resize (n_vertices, NULL);
        //w.resize (n_vertices, 0.0);
    }
    
    // -------------------------------------------------------------
    
    void SetGscoreOfLastVertex (const doubleType& gg) {
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
            B[k] = gsSq[1] - gsSq[kp1] -vsSq[1][0]; // \beta_k * vs[k+1][k]
            for (int p=0; p<kp1; ++p)  B[k] += vsSq[kp1][p]; // p=0,1,...,k
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
        
        doubleType C1_2=C1/C2, C0_2=C0/C2;
        doubleType CR = C1_2*C1_2 - 4*C0_2;
        //printf ("C0=%f, C1=%f, C2=%f, CR=%f\n", C0, C1, C2, CR);
        
        //printf("\nCR = %e\n", CR);
        
        if (CR < (doubleType)(-_MS_DOUBLE_EPS_SQ)) // _MS_DOUBLE_EPS degenerate simplex
        //if (CR < (doubleType)(0.0)) // degenerate simplex
            return (false);
        else if (CR < (doubleType)(_MS_DOUBLE_EPS_SQ))
            CR = (doubleType)0.0;
        o[0] = (sqrt(CR) - C1_2) / 2; //(-C1 + sqrt(C1*C1 - 4.0*C2*C0)) / (2.0*C2);
        for (int p=1; p<n_vertices_m1; ++p) // 1,2,...,n-2
            o[p] = A[p]*o[0] + B[p];
        
        return (true);
    }
    
    // ==============================================================
    
    // constructors
    AMetricSimplex (MetricSimplexContainerType* sp=NULL) : 
            backTrackSimplex(0),
            n_vertices(0), n_vertices_m1(-1), n_vertices_m2(-2), G(std::numeric_limits<doubleType>::max()),
            gs(doubleVecType(1,-1.0)), AllSimplices_p(NULL),
            simplex_computation_stage(0u), simplex_computation_failure(0u), ApexCameFromFace(NULL),
            isMaximal (false), G_cameFromPoint(std::numeric_limits<doubleType>::max()) { }
    
    AMetricSimplex (nodePointerType np, MetricSimplexContainerType* sp=NULL) 
            : backTrackSimplex(0),
              n_vertices(0), n_vertices_m1(-1), n_vertices_m2(-2), G(std::numeric_limits<doubleType>::max()),
              gs(doubleVecType(1,-1.0)), AllSimplices_p(sp), 
              simplex_computation_stage(COMPUTE_ALL), simplex_computation_failure(0u), // everything computed
              ApexCameFromFace(NULL), isMaximal (false), G_cameFromPoint(std::numeric_limits<doubleType>::max())
    // np = p[0]
    {
        MakeRoomForAnotherVertex (np);
        // vs, vsSq are empty; distMat and distSqMat are 1x1 containg 0; w has one element with 0; -- no need to set these
    }
    
    AMetricSimplex (nodePointerType just_created, nodePointerType came_from, 
                        doubleType d=std::numeric_limits<doubleType>::quiet_NaN(), MetricSimplexContainerType* sp=NULL) 
            : backTrackSimplex(0),
              n_vertices(0), n_vertices_m1(-1), n_vertices_m2(-2), G(std::numeric_limits<doubleType>::max()),
              gs(doubleVecType(1,-1.0)), AllSimplices_p(sp), 
              simplex_computation_stage(0u), simplex_computation_failure(0u), ApexCameFromFace(NULL),
              isMaximal (false), G_cameFromPoint(std::numeric_limits<doubleType>::max())
    // just_created = p[1], came_from = p[0]
    {
        MakeRoomForAnotherVertex (just_created);
        // vs, vsSq are empty; distMat and distSqMat are 1x1 containg 0; w has one element with 0; -- no need to set these
        
        MakeRoomForAnotherVertex (came_from);
        
        // Phase 1: Compute embedding of vertices
        // ----
        
        if (std::isnan(d)) {
            if ( std::isnan(d = came_from->getDistanceToSuccessor (just_created)) ) {
                simplex_computation_failure = CHECK_CONNECTION;
                return;
            }
        }
        
        SetDistancesFromLastVertex (_MS_SMALL_VECTOR<doubleType>(1,d));
        simplex_computation_stage |= CHECK_CONNECTION;
        
        // --
        if (ComputeCoordinateOfLastVertex()) 
            simplex_computation_stage |= COMPUTE_LOCAL_COORDINATE;
        else {
            simplex_computation_failure = COMPUTE_LOCAL_COORDINATE;
            return;
        }
        
        // Phase 2: Compute G-score
        // ----
        
        SetGscoreOfLastVertex (came_from->G);
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
        
        
        G = came_from->G + d;
        w[p[1]] = 1.0;
        simplex_computation_stage |= COMPUTE_G_SCORE_OF_APEX; // no failure possible here
        
        ApexCameFromFace = this;
    }
    
    // ==============================================================
    
    bool checkConnections (nodePointerType np, 
                                _MS_SMALL_VECTOR<doubleType>* dists_p=NULL, // return distances if not NULL
                                // doubleType const* dist_to_p0=NULL,
                                bool stopIfConnectionFails=false)
    {
        _dosl_verbose_head(1);
        
        // sets NaN for whicherver connection does not exist.
        // 'stopIfConnectionFails' is redundent if 'dists_p' is NULL
        bool allConnectionsExist = true;
        
        // Assume 'np->Successors' have been generated.
        for (int a=0; a<n_vertices; ++a) {
            // skip checking with p[0]
            /* if (a==0 && dist_to_p0) {
                if (dists_p)
                    dists_p->push_back (*dist_to_p0);
                continue;
            } */
            
            auto found_it = np->Successors.find (p[a]);
            
            if (found_it==np->Successors.end() || std::isnan(found_it->second)) { // np cannot be inserted in this simplex
                if (_dosl_verbose_on(0)) {
                    np->print("np = ");
                    p[a]->print("p[a] = ");
                    //printf ("(found_it==np->Successors.end())=%d, found_it->second=%f", (found_it==np->Successors.end()), found_it->second);
                }
                if (!stopIfConnectionFails  &&  dists_p) {
                    dists_p->push_back ( std::numeric_limits<doubleType>::quiet_NaN() );
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
        // Returns true if updated with the current G-scores, 
        //      otherwise (if returning previously computed list) return false
        
        if (AllCommonNeighbors.size()) {
            // check if G-scores changed
            bool GScoresUpdated = false;
            for (auto it=AllCommonNeighbors.begin(); it!=AllCommonNeighbors.end();++it)
                if (it->first->G != it->second) {
                    if (forceRecompute) {
                        GScoresUpdated = true;
                        it->second = it->first->G; // update 'AllCommonNeighbors' if it changed.
                    }
                    else
                        return (true);
                    /* NOTE: G-score of a vertex is the best G-score across all simplices, 
                        while G-score of simplex is the G-score of apex for a path through that simplex. */
                }
            // check if the G-scores are to be updated.
            return (GScoresUpdated);
        }

        // not previously computed. recompute from scratch.
        AllCommonNeighbors.clear();
        std::unordered_map<nodePointerType,doubleType>  tmpCommonNeighbors;
        // The successors of p[0] that are not part of inSimplex
        for (auto it=p[0]->Successors.begin(); it!=p[0]->Successors.end(); ++it)
            AllCommonNeighbors[it->first] = it->second; // Note: Successors is an unordered map
        for (int a=1; a<p.size(); ++a) {
            // remove p[a] if it existed in allCommonNeighbors
            AllCommonNeighbors.erase (p[a]);
            // now remove everything else that's not part of p[a].Successor
            tmpCommonNeighbors = AllCommonNeighbors; // so that iterator remains valid
            for (auto it = tmpCommonNeighbors.begin(); it != tmpCommonNeighbors.end(); ++it)
                if (p[a]->Successors.find(it->first) == p[a]->Successors.end())
                    AllCommonNeighbors.erase (it->first);
        }
        
        return (true);
    }
    
    // ----------------------------------------------------------------------------------------
    // computing shortest distance between p[0] and o, and corresponding weights, within the 
    
    AMetricSimplex* GetFaceSimplex (int i, set_of_pairs_ordered_by_second<int, doubleType>* FaceComputationOrder) {
        // compute face simplex coordinates and weights (recursive call to 'ComputeGScore')
        // FaceComputationOrder will contain i as well
        // 'FaceComputationOrder' is used to prioritize the order in which the sub-sub-simplices are created.
        _dosl_verbose_head(1);
        
        int a, b; // temporary variables
        if (_dosl_verbose_on(0)) {
            _dosl_printf ("GetFaceSimplex(%d) for simplex %x: ", i, this); // std::cout.flush();
        }
        
        // Test FaceSimplices[i]
        if (FaceSimplices[i]) {
            // No need to check if g-scores changed, because if they did, this would be a different simplex.
            if (_dosl_verbose_on(0)) {
                _dosl_printf ("Already computed. Returning %x.", FaceSimplices[i]);
            }
            return (FaceSimplices[i]);
        }
        
        // --------------
        // Search in 'AllMetricSimplexPointers'
        
        AMetricSimplex* found_simplex = NULL;
        AMetricSimplex* search_simplex = new AMetricSimplex; // temporary metric simplex for comparison
        search_simplex->p = p;
        search_simplex->gs = gs;
        search_simplex->n_vertices = p.size();
        
        if (_dosl_verbose_on(0)) {
            _dosl_printf ("Will do cascading search in AllSimplices table (size=%d). ", AllSimplices_p->size());
        }
        
        int ind, ind_in_face;
        _MS_SMALL_VECTOR<int> removed_vertex_indices;
        // search for the starting elementary simplex (note: VertexInsertionOrderReversed.size() == n_vertices-1 )
        typename set_of_pairs_ordered_by_second<int,doubleType>::iterator ind_it;
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
            found_simplex = AllSimplices_p->find (search_simplex);
            if (found_simplex) break;
        }
        
        delete search_simplex; // TODO: can use this for efficiency?
        
        if (ind==i && found_simplex) { // found the exact required face simplex
            FaceSimplices[i] = found_simplex;
            return (FaceSimplices[i]);
        }
        
        // Build the face simplex (insert one vertex at a time)
        
        AMetricSimplex* face_simplex;
        /* int ind, ind2, b;
        doubleVecType dists; */
        
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
            face_simplex = AllSimplices_p->createNewOneSimplex (p[0], p[ind], distMat[0][ind]);
            if (_dosl_verbose_on(0)) {
                face_simplex->print_pgs ("new 1-simplex created");
            }
        }
        
        for (a=removed_vertex_indices.size()-1; a>0; --a) { // removed_vertex_indices[0] == i. So don't insert that.
            ind = removed_vertex_indices [a]; // p[ind] to be inserted
            
            face_simplex = AllSimplices_p->checkConnectionsAndAddVertex (face_simplex, p[ind]);
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
        
        doubleVecType tw0 (n_vertices); // ignore index 0. use 1,2,...,n_vertices-1
        
        // Compute naieve weights
        doubleType wvSum, wSum=0.0;
        for (int j=n_vertices-1; j>0; --j) {
            wvSum = 0.0;
            for (int i=j+1; i<n_vertices; ++i)
                wvSum += tw0[i] * vs[i][j-1];
            tw0[j] = (o[j-1] - wvSum) / vs[j][j-1];
            wSum += tw0[j];
        }
        // ordered set.
        
        doubleType twt;
        for (int j=1; j<n_vertices; ++j) {
            twt = tw0[j]/wSum;
            tw.insert ( std::make_pair(j,twt) );
            w[p[j]] = twt;
        }
        doubleType wo = 1.0/wSum; // G of cameFromPoint = G * (1 - wo);
        
        if (_dosl_verbose_on(0)) {
            _dosl_printf_nobreak ("call to ComputeGScore (this=%x, depth=%d): n_vertices = %d, naive weights = [", 
                                                                    this, depth, n_vertices);
            for (auto it=tw.begin(); it!=tw.end(); ++it)
                printf(" (%x,%d,%f) ", p[it->first], it->first, it->second);
            printf("]\n");
        }
        // tw contains the weights in increasing order. tw.begin()->second is the smallest weight.
        
        if (tw.begin()->second >= 0.0) { // all weights are positive
            G = DistanceBetweenVertices (o, vs[0]);
            G_cameFromPoint = (1.0 - wo) * G;
            ApexCameFromFace = this;
            simplex_computation_stage |= COMPUTE_G_SCORE_OF_APEX;
            return (true);
        }
        
        //--
        // create subsimplex by removing a negative weight vertex.
        AMetricSimplex* face_simplex;
        G = std::numeric_limits<doubleType>::max();
        
        // int best_face_search_depth = max_depth_to_search-1, this_face_search_depth;
        
        for (auto it=tw.begin(); it!=tw.end(); ++it) {
            if (it->second < 0.0) { // need to check all faces that have negative weights
                
                if (FaceSimplices[it->first]) {
                    if ( !TO_BOOL(FaceSimplices[it->first]->simplex_computation_stage & COMPUTE_G_SCORE_OF_APEX)  && 
                           !TO_BOOL(FaceSimplices[it->first]->simplex_computation_failure) )
                       FaceSimplices[it->first]->ComputeGScore (/*distsFromP0,*/ depth+1);
                   face_simplex = FaceSimplices[it->first];
                }
                else {
                    // Polulate the list of vertices that make up this face
                    _DOSL_SMALL_VECTOR <nodePointerType> pFace = p;
                    pFace.erase (it->first);
                    face_simplex = AllSimplices_p->constructSimplexFromVertices (pFace, /*distsFromP0,*/ COMPUTE_ALL); 
                                                                                        // compute G-score as well
                    FaceSimplices[it->first] = face_simplex;
                }
                
                // Compare
                if (face_simplex  &&  TO_BOOL(face_simplex->simplex_computation_stage & COMPUTE_G_SCORE_OF_APEX)  &&
                        face_simplex->G < G ) {
                    G = face_simplex->G;
                    G_cameFromPoint = face_simplex->G_cameFromPoint;
                    // Update weight assignment.
                    w[p[it->first]] = 0.0;
                    for (auto itW=face_simplex->w.begin(); itW!=face_simplex->w.end(); ++itW) // copy weights
                        w[itW->first] = itW->second;
                    ApexCameFromFace = face_simplex->ApexCameFromFace;
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
    
    void SetChildInfluence (void) { // to be called on new camefrom simple, when updating G-score
        for (int a=1; a<p.size(); ++a)
            p[a]->ChildrenInfluenced.insert (p[0]);
    }
    
    void UnsetChildInfluence (void) { // to be called on previous camefrom simple, when updating G-score
        for (int a=1; a<p.size(); ++a)
            p[a]->ChildrenInfluenced.erase (p[0]);
    }
};

#endif
