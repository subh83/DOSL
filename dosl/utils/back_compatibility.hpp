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

#ifndef __DOSL_BACK_COMP_HPP
#define __DOSL_BACK_COMP_HPP

// v3.1 -> 3.2
#define AStarNode               AStar::Node
#define AStarProblem            AStar::Algorithm
#define SStarNode               SStar::Node
#define SStarProblem            SStar::Algorithm
#define ThetaStarNode           ThetaStar::Node
#define ThetaStarProblem        ThetaStar::Algorithm
#define getPointerPathToNode    reconstructPointerPath

// v3.2x -> 3.3
//#define AllNodesSet             all_nodes_set_p.operator*()
//#define NodeHeap                node_heap_p.operator*()
// renames
#define set_hash_table_size              reserve
#define G                                g_score
#define F                                f_score
#define PostHashInsertInitiated          post_hash_insert_initiated
#define lineageData                      lineage_data
#define Successors                       successors
#define SuccessorsCreated                successors_created
#define CameFrom                         came_from
#define Expanded                         expanded
#define Node_hasher                      NodeHasherFunc
#define Node_equal_to                    NodeEqualToFunc
#define AllNodesSet_p                    all_nodes_set_p
#define Node_key_less_than               NodeKeyLessThanFunc
#define NodeHeap_p                       node_heap_p
#define chrono_timer                     ChronoTimer
#define HeapPos                          heap_pos
#define ExpandCount                      expand_count
#define SubopEps                         subopt_eps
#define ProgressShowInterval             progress_show_interval
#define StartNodes                       start_nodes
#define nodeEventType                    NodeEventType
#define BookmarkNodePointers             bookmarked_node_pointers
#define reconstructPath                  reconstruct_path
#define reconstruct_path                 reconstruct_weighted_pointer_path
#define reconstructPointerPath           reconstruct_pointer_path
#define getNodePointer                   get_node_pointer
#define getCostsToNodes                  get_costs_to_nodes
#define getBookmarkNodePointers          get_bookmark_node_pointers
#define costType                         _CostType
#define GenerateSuccessors               generate_successors
#define thisSuccessors                   this_successors
#define thisTransitionCosts              this_transition_costs
#define nodeInHash_p                     node_in_hash_p
#define thisNeighbourNodeInHash_p        this_neighbour_node_in_hash_p
#define nodePointerType                  _NodePointerType
#define ChildrenInfluenced               children_influenced
#define doubleType                       DoubleType
#define doubleVecType                    DoubleVecType
#define vsSq                             vs_sq
#define distMat                          dist_mat
#define distSqMat                        dist_sq_mat
#define thisDistSq                       this_dist_sq
#define backTrackSimplex                 back_track_simplex
#define ApexCameFromFace                 apex_came_from_face
#define isMaximal                        is_maximal
#define AttachedAbsoluteMaximalSimplices attached_absolute_maximal_simplices
#define AllCommonNeighbors               all_common_neighbors
#define G_cameFromPoint                  g_came_from_point
#define FaceSimplices                    face_simplices
#define AllSimplices_p                   all_simplices_p
#define getEmptySimplex                  get_empty_simplex
#define createNewZeroSimplex             create_new_zero_simplex
#define createNewOneSimplex              create_new_one_simplex
#define constructSimplexFromVertices     construct_simplex_from_vertices
#define checkConnectionsAndAddVertex     check_connections_and_add_vertex
#define getAllCommonNeighborsOfSimplex   get_all_common_neighbors_of_simplex
#define getAllAttachedMaximalSimplices   get_all_attached_maximal_simplices
#define findCamefromPoint                find_camefrom_point
#define inPathPoint                      in_path_point
#define forceCompute                     force_compute
#define nodeSet_p                        node_set_p
#define inSimplex_p                      in_simplex_p
#define AllMetricSimplexPointers         all_metric_simplex_pointers
#define metricSimplexPointerType         MetricSimplexPointerType
#define compareType                      CompareType
#define pair_type                        PairType
#define set_of_pairs_ordered_by_second   SetOfPairsOrderedBySecond
#define DistanceBetweenVertices          distance_between_vertices
#define HashFunctorInstance              hash_functor_instance
#define EqualToFunctorInstance           equal_to_functor_instance
#define HashTableSize                    hash_table_size






#endif
