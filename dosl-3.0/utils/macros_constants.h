#ifndef __DOSL_MACROS_CONSTANTS_H
#define __DOSL_MACROS_CONSTANTS_H

// =================================================

// non-standard libraries
// #define _DOSL_LIB_ARMADILLO      0

// =================================================

#include<stdexcept>
#include<string>

#ifndef _DOSL_DEBUG
#define _DOSL_DEBUG 0
#endif

#ifndef _DOSL_AUTOCORRECT
#define _DOSL_AUTOCORRECT 1
#endif

#ifndef _DOSL_PRINT_COLORS
#define _DOSL_PRINT_COLORS 1
#endif

#ifndef _DOSL_VERBOSE
#define _DOSL_VERBOSE 1
#endif

#ifndef _DOSL_EVENTHANDLER
#define _DOSL_EVENTHANDLER 1
#endif

#ifndef _dosl_err
#define _dosl_err(...) { std::cout << std::flush; char tmpstr[1024]; sprintf(tmpstr, __VA_ARGS__); throw std::runtime_error(tmpstr); }
#endif

#ifndef _dosl_warn
#define _dosl_warn(...) { std::cout << "WARNING: "; printf(__VA_ARGS__); std::cout << std::flush; }
#endif

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
    /* #include "utils/fast_vector.h"
    #define _DOSL_FAST_VECTOR  fast_vector */

    // fast_vector<..., 0u>
    /* #include "utils/fast_vector.h"
    template <class T>
    using _DOSL_FAST_VECTOR = fast_vector<T, 0u>; // C++11 */

    // std::vector
    #include "utils/std_vector.h"
    #define _DOSL_FAST_VECTOR  std_vector
    
#endif

// ---------------------------
// std::vector or equivalent
// Should also provide find, rfind, findi, rfindi, erase

#ifndef _DOSL_SMALL_VECTOR

    // fast_vector<..., 1u>
    #include "utils/fast_vector.h"
    #define _DOSL_SMALL_VECTOR  fast_vector

    // fast_vector<..., 0u>
    /* #include "utils/fast_vector.h"
    template <class T>
    using _DOSL_FAST_VECTOR = fast_vector<T, 0u>; // C++11 */

    // std::vector
    /* #include "utils/std_vector.h"
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
    #include "utils/fast_unordered_set.h"
    #define _DOSL_LARGE_UNORDERED_SET        fast_unordered_set

#endif

// ---------------------------

#ifndef _DOSL_HEAP

    #include "utils/binary_heap.h"
    #define _DOSL_HEAP                   binary_heap  //std::priority_queue
    // _DOSL_HEAP <key, CompareFunctor, HeapPosFunctor>
    // Provides: insert (key&);

#endif

// =================================================

int RUNTIME_VERBOSE_SWITCH = 1;

#ifndef _DOSL_VERBOSE_ITEMS
#define _DOSL_VERBOSE_ITEMS "ALL" // Comma-separated list containing "ALL", "NONE", 'func_name' or '!func_name'
#endif

constexpr int _is_func_in_list (char const * needle, char const * haystack=_DOSL_VERBOSE_ITEMS,
                            int const needle_pos=0, int const haystack_pos=0) {
    return ( ((haystack [haystack_pos] == '\0') && (needle [needle_pos] != '\0'))  // not in list
             ? 0
             : ( ((needle [needle_pos] == '\0')) 
                 ? ( (haystack_pos>needle_pos  &&  haystack[haystack_pos-needle_pos-1]=='!')
                     ? -1   // in list, but negated
                     : 1 )  // in list
                 : ( (needle [needle_pos] == haystack [haystack_pos])
                     ? _is_func_in_list (needle, haystack, needle_pos+1, haystack_pos+1) // check next character
                     : _is_func_in_list (needle, haystack, 0, haystack_pos+1) ) ) ); // start afresh
}

constexpr bool _is_verbose_on (char const * funcname, unsigned int const level=0) {
    return ( (_DOSL_VERBOSE==0 || _DOSL_VERBOSE>=level)  && 
             ( ( _is_func_in_list ("ALL")>0  &&  _is_func_in_list (funcname)>=0 )   ||
               ( _is_func_in_list ("ALL")<=0 &&  _is_func_in_list ("NONE")>=0  &&  _is_func_in_list (funcname)>0 )
             )
           );
}

#define _if_verbose_on(f,l)  if(RUNTIME_VERBOSE_SWITCH && _is_verbose_on(f,l))

// ---------------------

#include <string>
#include <iostream>

// Tab Printing
unsigned int GLOBAL_INDENTATION = 0; // to keep track of tab indentation in printing
void increase_indentation (void) { ++GLOBAL_INDENTATION; }
void decrease_indentation (void) { if (GLOBAL_INDENTATION>0) --GLOBAL_INDENTATION; }
void print_indentation (int relTabs = 0) {
    int nTabs = GLOBAL_INDENTATION + relTabs;
    std::string tabStr = "";
    for (int a=0; a<nTabs; ++a)
        tabStr += "\t";
    std::cout << tabStr;
}
#define _INDENT     increase_indentation(); // print_indentation();
#define _NEWLINE_   printf("\n"); print_indentation();
#define INDENT_     decrease_indentation();

#define return_INDENT_(x)  { INDENT_; return((x)); }

// print styles codes:
#if _DOSL_PRINT_COLORS
    #define _PRE "\033["
    #define _RED        _PRE "31m"
    #define RED_        _PRE "39m"
    #define _GREEN        _PRE "32m"
    #define GREEN_        _PRE "39m"
    #define _YELLOW        _PRE "33m"
    #define YELLOW_        _PRE "39m"
    #define _BLUE        _PRE "34m"
    #define BLUE_        _PRE "39m"
#else
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


std::string _uint_to_binary (unsigned int x) {
    std::string  binStr="";
    for (unsigned int z=1; z<=x; z<<=1)
        binStr = (((x & z)==0)?"0":"1") + binStr;
    return binStr;
}

#define uint_to_binary(x) (_uint_to_binary(x).c_str())


#endif
