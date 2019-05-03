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
    
    
};

#endif
