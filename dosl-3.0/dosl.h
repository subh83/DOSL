/** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal Search Library (DOSL)                                              *
*        [A branch of DOSL]                                                           *
*    A template-based C++ library for discrete (graph) search                            *
*    Version 3.0                                                                         *
*    ----------------------------------------------------------                          *
*    Copyright (C) 2016  Subhrajit Bhattacharya                                          *
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
*    Contact:  subhrajit@gmail.com ,  http://subhrajit.net/                              *
*                                                                                        *
*                                                                                        *
*************************************************************************************** **/
#ifndef __DOSL_H
#define __DOSL_H

/* *** Helper macro for selecting planner ***

Set 'DOSL_ALGORITHM' before including this file. Otherwise, multiple algorithm files will be included.
Ex:
    #define _DOSL_ALGORITHM  AStar
    #include "dosl.h"

If set, also provides macro 'DOSL_ALGORITHM(str)'
Ex:
    DOSL_ALGORITHM(Node)
expands to
    AStarNode
*/

#ifdef _DOSL_ALGORITHM

    #define QMAKESTR(x) #x
    #define MAKESTR(x) QMAKESTR(x)
    #define EVAL(x) x
    #define MAKEINC(x) planners/EVAL(x).h

    // include:
    #include MAKESTR(MAKEINC(_DOSL_ALGORITHM))
    
    // macro 'DOSL_ALGORITHM'
    #define QJOIN(x, y) x ## y
    #define JOIN(x, y) QJOIN(x, y)
    #define DOSL_ALGORITHM(trail) JOIN(_DOSL_ALGORITHM,trail)
    
#else
    
    #define _DOSL_ALGORITHM  UndefinedAlgorithm
    
    #include "planners/AStar.h"
    // #include "planners/SStar.h"
    
    // this cannot be used as usual:
    #define DOSL_ALGORITHM(trail) JOIN(,trail)
    
#endif

#endif
