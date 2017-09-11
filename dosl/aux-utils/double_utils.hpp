/** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal search Library (DOSL)                                              *
*    A template-based C++ library for discrete search                                    *
*    Version 3.x                                                                         *
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
*    Contact:  subhrajit@gmail.com                                                       *
*              https://www.lehigh.edu/~sub216/ , http://subhrajit.net/                   *
*                                                                                        *
*                                                                                        *
*************************************************************************************** **/
#ifndef _DOSL_DOUBLE_UTILS_HPP
#define _DOSL_DOUBLE_UTILS_HPP

#include <math.h>

// ------------------------------------------------------------
// generic macros

#define isEqual_i(x,y)  ((x)==(y))
#define isLess_i(x,y)  ((x)<(y)) // x < y
#define isGreater_i(x,y)  ((x)>(y)) // x > y
#define isLessEq_i(x,y)  ((x)<=(y)) // x <= y
#define isGreaterEq_i(x,y)  ((x)>=(y)) // x <= y

#define sign(x)    (((x)>0.0)?1.0:(((x)<0.0)?-1.0:0.0))
#define iround(x)  ((int)round(x))

// ------------------------------------------------------------
// relaxed comparisons and other macros (appropriate if x and y take discrete values)

#ifndef INFINITESIMAL_DOUBLE
#define INFINITESIMAL_DOUBLE  1e-8
#endif

#define isEqual_d(x,y)  ( fabs((x)-(y)) < INFINITESIMAL_DOUBLE ) // x == y
#define isLess_d(x,y)  ( (x)+INFINITESIMAL_DOUBLE < (y) ) // x < y
#define isGreater_d(x,y)  ( (x) > (y)+INFINITESIMAL_DOUBLE ) // x > y
#define isLessEq_d(x,y)  (!isGreater_d(x,y)) // x <= y
#define isGreaterEq_d(x,y)  (!isLess_d(x,y)) // x <= y

#define sign_d(x)   ((isGreater_d((x),0.0))?1.0:((isLess_d((x),0.0))?-1.0:0.0))

// ------------------------------------------------------------
// math constnts
#define PI       3.1415926535897931
#define PI_BY_3  1.0471975511965976
#define SQRT3BY2 0.8660254037844386
#define SQRT2    1.4142135623730951
#define SQRT3    1.7320508075688772

// ------------------------------------------------------------
// Functions

int approx_floor (double x, double tol=INFINITESIMAL_DOUBLE) {
    double ret = floor (x);
    if (ret+1.0-x<tol) ++ret;
    return ((int)ret);
}

#endif
