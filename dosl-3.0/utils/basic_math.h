#ifndef __DOSL_BASIC_MATH_H_
#define __DOSL_BASIC_MATH_H_

#include "utils/macros_constants.h"

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
