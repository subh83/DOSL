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

#ifndef __DOSL_VERBOSE_TOOLS_TCC
#define __DOSL_VERBOSE_TOOLS_TCC

#include "macros_constants.tcc"

// Tab Printing


static int _DOSL_VERBOSE_FUN_DEPTH = 0;
static int DOSL_RUNTIME_VERBOSE_SWITCH = 1;
static bool _dosl_is_tab_printed = false;


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
    return ( 
             // "ALL" or "ALL,func_name", but not "ALL,!func_name":
             ( _is_func_in_list ("ALL")>0 && _is_func_in_list (funcname)>=0 ) || 
             // "NONE,func_name"
             ( _is_func_in_list ("NONE")>0 && _is_func_in_list (funcname)>0 ) || 
             // verbose by level:
             ( ( _is_func_in_list ("BY_LEVEL")>0 || (_is_func_in_list ("ALL")==0 && _is_func_in_list ("NONE")==0) ) && 
                  ( level<=_DOSL_VERBOSE_LEVEL || _is_func_in_list (funcname)>0 ) && // with level<=_DOSL_VERBOSE_LEVEL or "BY_LEVEL,func_name"
                 !(_is_func_in_list (funcname)<0) ) // but not "BY_LEVEL,!func_name"
           );
}


class dosl_verbose_function {
public:
    char funcname[256];
    bool isInitiated;
    int verbose_level, level_increment, subblock_increment;
    
    dosl_verbose_function() : level_increment(0), subblock_increment(0), isInitiated(false) {}
    dosl_verbose_function (const char* fn, int dl=1) { init (fn, dl); }
    
    void init (const char* fn, int dl=1) {
        strcpy(funcname,fn);
        level_increment = dl;
        subblock_increment = 0;
        _DOSL_VERBOSE_FUN_DEPTH += level_increment;
        verbose_level = _DOSL_VERBOSE_FUN_DEPTH;
        isInitiated = true;
    }
    
    bool set_subblock_increment (int dl=0) {
        subblock_increment = dl;
        _DOSL_VERBOSE_FUN_DEPTH = verbose_level;
        _DOSL_VERBOSE_FUN_DEPTH += subblock_increment;
        return (true);
    }
    
    ~dosl_verbose_function() {
        if (level_increment) {
            _DOSL_VERBOSE_FUN_DEPTH -= level_increment + subblock_increment;
            printf ("\n");
        }
    }
    
    bool to_true (void) { return true; }
};

// ----

// Typical use: "_dosl_verbose_head(1)" at the start of the function.
// #define _dosl_verbose_head(dl) false;
#if _DOSL_VERBOSE_LEVEL<0
#define _dosl_verbose_head(dl) false;
#else
#define _dosl_verbose_head(dl) \
                dosl_verbose_function  fun_verbose; \
                if ( DOSL_RUNTIME_VERBOSE_SWITCH  &&  _is_verbose_on(__func__, _DOSL_VERBOSE_FUN_DEPTH+(dl)) ) { \
                    fun_verbose.init (__func__,dl); \
                    printf("\n"); print_indentation(0); \
                    printf(_BLUE "function: %s. level: %d." BLUE_, fun_verbose.funcname, fun_verbose.verbose_level); \ 
                    std::cout << std::endl; \
                }
#endif

// Typical use: "if (_dosl_verbose_on(0)) { ... }".
#if _DOSL_VERBOSE_LEVEL<0
#define _dosl_verbose_on(dl) false
#else
#define _dosl_verbose_on(dl) (DOSL_RUNTIME_VERBOSE_SWITCH && fun_verbose.isInitiated \
                                && _is_verbose_on(fun_verbose.funcname, fun_verbose.verbose_level+(dl)) \
                                && fun_verbose.set_subblock_increment(dl) )
#endif

// ---------------------------

template <class IntType>
std::string get_indentation_string (IntType relTabs = 0) {
    int nTabs = _DOSL_VERBOSE_FUN_DEPTH + relTabs;
    std::string tabStr = "";
    for (int a=0; a<nTabs; ++a)
        tabStr += "\t";
    return (tabStr);
}

template <class IntType>
void print_indentation (IntType relTabs = 0) {
    std::cout << get_indentation_string(relTabs);
}

/* void increase_indentation (void) { ++GLOBAL_INDENTATION; }
void decrease_indentation (void) { if (GLOBAL_INDENTATION>0) --GLOBAL_INDENTATION; }
#define _INDENT     increase_indentation(); // print_indentation();
#define INDENT_     decrease_indentation(); */

#define DOSL_INDENT    print_indentation(0);
#define DOSL_NEWLINE   printf("\n"); print_indentation(0);

// --

#ifndef _dosl_printf
#define _dosl_printf(...) { if(!_dosl_is_tab_printed) DOSL_INDENT; \
                            printf(__VA_ARGS__); std::cout << std::endl; \
                            _dosl_is_tab_printed = false; }
#endif

#ifndef _dosl_printf_nobreak
#define _dosl_printf_nobreak(...) { if(!_dosl_is_tab_printed) { DOSL_INDENT; _dosl_is_tab_printed=true;} \
                                    printf(__VA_ARGS__); }
#endif

#ifndef _dosl_linebreak
#define _dosl_linebreak { printf("\n"); DOSL_INDENT; _dosl_is_tab_printed=true; }
#endif

// --

/*class TAB_TRACKED_OSTREAM : public std::ostream {
public:
    bool isTabPrinted;
    TAB_TRACKED_OSTREAM (std::ostream oo) : std::ostream(oo) { isTabPrinted= false; }
};

TAB_TRACKED_OSTREAM _dosl_cout (std::cout); */

/*std::ostream _dosl_cout_ (void) {
    std::ostream tmp;
    if (!_dosl_is_tab_printed) { tmp << ; _dosl_is_tab_printed=true;};
    
} */

#ifndef _dosl_cout
#define _dosl_cout  { if(!_dosl_is_tab_printed) { DOSL_INDENT; _dosl_is_tab_printed=true;} } std::cout 
#endif

#ifndef _dosl_endl
#define _dosl_endl  std::endl; { _dosl_is_tab_printed=false; }
#endif

#ifndef _dosl_newl
#define _dosl_newl  std::endl << get_indentation_string(0)
#endif

// ====================================================

template <class U>
std::string _uint_to_binary (U _x) {
    unsigned int x = (unsigned int) _x;
    std::string  binStr="";
    for (unsigned int z=1; z<=x; z<<=1)
        binStr = (((x & z)==0)?"0":"1") + binStr;
    return binStr;
}

#define uint_to_binary(x) (_uint_to_binary(x).c_str())

#endif

