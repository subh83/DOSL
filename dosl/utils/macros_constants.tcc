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

#ifndef __DOSL_MACROS_CONSTANTS_TCC
#define __DOSL_MACROS_CONSTANTS_TCC

// =================================================

// non-standard libraries
// #define _DOSL_LIB_ARMADILLO      0

// =================================================

#include <stdio.h>
#include <string>
#include <string.h>
#include <iostream>
#include<stdexcept>
#include<string>
#include <stdarg.h>
#include <chrono>

#ifndef _DOSL_DEBUG
#define _DOSL_DEBUG 1
#endif

#ifndef _DOSL_AUTOCORRECT
#define _DOSL_AUTOCORRECT 1
#endif

#ifndef _DOSL_PRINT_COLORS
#define _DOSL_PRINT_COLORS 1
#endif

#ifndef _DOSL_EVENTHANDLER
#define _DOSL_EVENTHANDLER 1
#endif

// =================================================

// print styles codes:

#if _DOSL_PRINT_COLORS
    #define _PRE            "\033["
    #define CLEAR_          _PRE "0m"
    // weights and styles
    #define _BOLD           _PRE "1m"
    #define BOLD_           CLEAR_ // _PRE "21m"
    //colors
    #define _RED            _PRE "31m"
    #define RED_            CLEAR_ // _PRE "39m"
    #define _GREEN          _PRE "32m"
    #define GREEN_          CLEAR_ // _PRE "39m"
    #define _YELLOW         _PRE "33m"
    #define YELLOW_         CLEAR_ // _PRE "39m"
    #define _BLUE           _PRE "34m"
    #define BLUE_           CLEAR_ // _PRE "39m"
#else
    #define _BOLD
    #define BOLD_
    #define _PRE
    #define _RED
    #define RED_
    #define _GREEN
    #define GREEN_
    #define _YELLOW
    #define YELLOW_
    #define _BLUE
    #define BLUE_
#endif


// Errors and warnings:

#ifndef _dosl_err
#define _dosl_err(...) { std::cout << std::flush << _RED _BOLD "ERROR: " BOLD_ ; char tmpstr[1024]; sprintf(tmpstr, __VA_ARGS__); std::cout << RED_ << std::endl; throw std::runtime_error(tmpstr); }
#endif

#ifndef _dosl_warn
#define _dosl_warn(...) { std::cout << std::flush << _YELLOW _BOLD "WARNING: " BOLD_; printf(__VA_ARGS__); std::cout << YELLOW_  << std::endl; }
#endif

#ifndef _dosl_info
#define _dosl_info(...) { std::cout << std::flush << _BLUE _BOLD "NOTE: " BOLD_; printf(__VA_ARGS__); std::cout << BLUE_  << std::endl; }
#endif

#if _DOSL_DEBUG
    #ifndef _dosl_warn_once
        #define _dosl_warn_once(warn_name,...) { \
            static std::unordered_map <std::string, int> _dosl_warn_counts; \
            if (_dosl_warn_counts[warn_name] == 0) { _dosl_warn(__VA_ARGS__); _dosl_warn_counts[warn_name] = 1; } \
        }
        #define _dosl_default_fun_warn(fun_name) { \
            _dosl_warn_once(fun_name, "Member '" \
            fun_name "' has not been overwritten. Using default member (this will most likely produce undesirable results)!"); \
        }
    #endif
#else
     #define _dosl_warn_once(warn_name,...) { }
     #define _dosl_default_fun_warn(fun_name) { }
#endif

// =================================================

#ifndef _DOSL_VERBOSE_LEVEL
#define _DOSL_VERBOSE_LEVEL 0 // -1: Turn ALL verbose off. 0: No verbose, except by '_DOSL_VERBOSE_ITEMS'; 1,2,...: verbose levels
#endif

#ifndef _DOSL_VERBOSE_ITEMS
#define _DOSL_VERBOSE_ITEMS "BY_LEVEL" // Comma-separated list starting with "ALL", "NONE" or "BY_LEVEL" (default), and
                                       //                     a list of functions with items 'func_name', '!func_name'
/* If 'func_name' appears in list, its verbose will be turned on irrespective of value of _DOSL_VERBOSE_LEVEL.
   If '!func_name' appears in list, its verbose will be turned off irrespective of value of _DOSL_VERBOSE_LEVEL. */
#endif

#include "verbose_tools.tcc"

// =================================================

class ChronoTimer {
public:
    std::chrono::high_resolution_clock::time_point t;
    void start(void) { t = std::chrono::high_resolution_clock::now(); }
    double read(void) {
        return (double)(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-t).count()) / 1e9;
    }
};

// =================================================

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#define TO_BOOL(v) ((v)!=0)

// =================================================
// ---------------------------
// std::vector or equivalent

#ifndef _DOSL_FAST_VECTOR
    
    // fast_vector<..., 1u>
    /* #include "fast_vector.tcc"
    #define _DOSL_FAST_VECTOR  fast_vector */

    // fast_vector<..., 0u>
    /* #include "fast_vector.tcc"
    template <class T>
    using _DOSL_FAST_VECTOR = fast_vector<T, 0u>; // C++11 */

    // std::vector
    #include "std_vector.tcc"
    #define _DOSL_FAST_VECTOR  std_vector
    
#endif

// ---------------------------
// std::vector or equivalent
// Should also provide find, rfind, findi, rfindi, erase

#ifndef _DOSL_SMALL_VECTOR

    // fast_vector<..., 1u>
    #include "fast_vector.tcc"
    #define _DOSL_SMALL_VECTOR  fast_vector

    // fast_vector<..., 0u>
    /* #include "fast_vector.tcc"
    template <class T>
    using _DOSL_FAST_VECTOR = fast_vector<T, 0u>; // C++11 */

    // std::vector
    /* #include "std_vector.tcc"
    #define _DOSL_FAST_VECTOR  std_vector */
    
#endif

// ---------------------------
// std::unordered_map or equivalent

#ifndef _DOSL_SMALL_MAP

    // Unordered map
    #include <unordered_map>
    #define _DOSL_SMALL_MAP              std::unordered_map // 'map', 'unordered_map'
    #define _DOSL_SMALL_MAP_pairfun      std::make_pair

#endif

// ---------------------------
// std::unordered_set or equivalent
// Should also provide: get (find and create if not exist)

#ifndef _DOSL_LARGE_UNORDERED_SET

    // HashTableContainer
    #include "fast_unordered_set.tcc"
    #define _DOSL_LARGE_UNORDERED_SET        fast_unordered_set
    
    /*// unfortunately following does not work since, once inserted, elements cannot be modified.
    // todo in future: Make search data 'mutable', while user data non-mutable? Bad idea?
    #include "std_unordered_set.tcc"
    #define _DOSL_LARGE_UNORDERED_SET       std_unordered_set */

#endif

// ---------------------------

#ifndef _DOSL_HEAP

    #include "binary_heap.tcc"
    #define _DOSL_HEAP                   binary_heap  //std::priority_queue
    // _DOSL_HEAP <key, CompareFunctor, HeapPosFunctor>
    // Provides: insert (key&);

#endif


// =================================================
// MEMORY

#define COPY_IF_NOTNULL_ELSE_CREATE_POINTER_TO_LOCAL(type,inPointer,copyPointer,defaultVal) \
                type _localvalue_ ## copyPointer; \
                type* copyPointer = & _localvalue_ ## copyPointer; \
                if (inPointer) copyPointer = inPointer; \
                else *copyPointer = defaultVal;

// ----------------------------

// To allocate memory when only pointer_type is known.
template <class T>
T* _new_p (T* dummy) {
    return (new T);
}
#define new_p(pointer_type)  _new_p(static_cast<pointer_type>(NULL))

#endif
