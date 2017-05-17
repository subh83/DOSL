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
#ifndef __DOSL_BASIC_MATH_TCC_
#define __DOSL_BASIC_MATH_TCC_

#include "macros_constants.tcc"

// ============================================

// Vector operations
// T = int, double, etc.

#include<vector>
#include<cmath>

// Vector addition/substraction

template <class T> 
_DOSL_SMALL_VECTOR<T> operator+ (const _DOSL_SMALL_VECTOR<T>& l, const _DOSL_SMALL_VECTOR<T>& r) {
    _DOSL_SMALL_VECTOR<T> ret = l;
    for (int a=0; a<l.size(); ++a)
        ret[a] += r[a];
    return (ret);
}

template <class T> 
_DOSL_SMALL_VECTOR<T> operator- (const _DOSL_SMALL_VECTOR<T>& l, const _DOSL_SMALL_VECTOR<T>& r) {
    _DOSL_SMALL_VECTOR<T> ret = l;
    for (int a=0; a<l.size(); ++a)
        ret[a] -= r[a];
    return (ret);
}

// Scalar multiplication

template <class T> 
_DOSL_SMALL_VECTOR<T> operator* (const T& l, const _DOSL_SMALL_VECTOR<T>& r) {
    _DOSL_SMALL_VECTOR<T> ret = r;
    for (int a=0; a<r.size(); ++a)
        ret[a] *= l;
    return (ret);
}

template <class T> 
_DOSL_SMALL_VECTOR<T> operator* (const _DOSL_SMALL_VECTOR<T>& l, const T& r) {
    _DOSL_SMALL_VECTOR<T> ret = l;
    for (int a=0; a<l.size(); ++a)
        ret[a] *= r;
    return (ret);
}

// Norm

template <class T>
T norm (const _DOSL_SMALL_VECTOR<T>& v, int p=2) {
    T ret = ((T)0.0);
    if (p==2) {
        for (int a=0; a<v.size(); ++a)
            ret += v[a]*v[a];
        return (sqrt(ret));
    }
}

// ============================================

#endif
