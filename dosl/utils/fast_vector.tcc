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

#ifndef __DOSL_FAST_VECTOR_TCC_
#define __DOSL_FAST_VECTOR_TCC_

#include<vector>

// ================================================================================================
// Efficient vector
template<class T, unsigned int level=1u> class fast_vector;

// Level 0:
//    - Constructors and destructors of T are not called explicitly.
//    - To construct newly reserved items, you need to call construct()
//    - Need to call reserve() manually before calling push_back(). No bound check is done.
//    - No bound check done by at() or operator[].
//    - If capacity is reduced below size by reserve(), size is not changed.
template<class T>
class fast_vector<T, 0u>
{
public:
    T* data;
    size_t _capacity, old_capacity;
    size_t _size;
    
    fast_vector() : data(NULL), _capacity(0u), old_capacity(0u), _size(0u) { }
    fast_vector(size_t n) : old_capacity(0u), _size(n) { _capacity=MAX(n,1u); data = (T*)malloc(sizeof(T)*_capacity); }
    fast_vector(size_t n, T const& item) : old_capacity(0u), _size(n) 
        { _capacity=MAX(n,1u); data = (T*)malloc(sizeof(T)*_capacity); T* _end=data+_size; for (T* dp=data; dp!=_end; ++dp) *dp=item; }
    
    inline void reserve(size_t n) { data = (T*)realloc(data, sizeof(T)*n); old_capacity = _capacity; _capacity = n; 
                                    for(int a=old_capacity;a<_capacity;++a) data[a]=T(); /*initialize*/ }
    inline void construct(T const& val = T()) { for (int a=old_capacity; a<_capacity; a++) data[a] = val; }
    inline void destroy(void) { free(data); data = NULL; _capacity = 0u; _size = 0u; }
    
    inline size_t size(void) const { return (_size); }
    inline size_t capacity(void) const { return (_capacity); }
    inline T* begin(void) const { return (data); }
    inline T* end(void) const { return (data+_size); }
    
    inline void push_back(T const& item) { data[_size++] = item; }
    // inline void insert(T const& item) { data[_size++] = item; } // alias of push_back
    inline void pop_back(void) { --_size; }
    inline void clear(void) { _size = 0u; }
    inline void resize(size_t n) { _size = n; }
    inline void resize(size_t n, T const& item) { for (int a=_size; a<n; ++a) data[a]=item;  _size = n; }
    inline T& at(size_t n) const { return(data[n]); }
    inline T& operator[](size_t n) const { return(data[n]); } // does not check bound
    inline T& back(void) const { return(data[_size-1]); } // does not check bound
    
    // Functions not native to std::vector
    template <class Pred=std::equal_to<T> >
    inline T* find(T const& item) const 
        { T* _end=end(); Pred _pred; for (T* dp=data; dp!=_end; ++dp) { if (_pred(*dp,item)) return (dp); } return (_end); }
    template <class Pred=std::equal_to<T> >
    inline T* rfind(T const& item) const 
        { T* _end=end(); T* dp=_end; Pred _pred; do { --dp; if (_pred(*dp,item)) return (dp); } while (dp!=data); return (_end); }
    
    template <class Pred=std::equal_to<T> >
    inline int findi(T const& item) const 
        { Pred _pred; for (int a=0; a<size(); ++a) { if (_pred(operator[](a),item)) return (a); } return (-1); }
    template <class Pred=std::equal_to<T> >
    inline int rfindi(T const& item) const 
        { Pred _pred; for (int a=size()-1; a>=0; --a) { if (_pred(operator[](a),item)) return (a); } return (-1); }
    
    inline void erase (size_t n)
        { if (n<_size) { T* _end=data+(_size-1); for (T* dp=data+n; dp!=_end; ++dp) *dp=*(dp+1); --_size; } }
    
    // Constructors, destructor, etc
    ~fast_vector() { free(data); }
    
    fast_vector(fast_vector const& other) : data(NULL), _capacity(0u), _size(0u)
        { if (other.data) { reserve(other._capacity); old_capacity = other.old_capacity; for(_size=0;_size<other._size;++_size) data[_size]=other.data[_size]; } }
    fast_vector& operator=(fast_vector const& rhs) 
        { if (rhs.data) { data=NULL; reserve(rhs._capacity); old_capacity = rhs.old_capacity; for(_size=0;_size<rhs._size;++_size) data[_size]=rhs.data[_size]; } return *this; }
    void copy_shared(fast_vector const& rhs) // shared data copy. warning: destroying one will invalidate other.
        { free(data); data = rhs.data; old_capacity = rhs.old_capacity; _capacity = rhs._capacity; _size = rhs._size; }
};

// Level 1:
// All features of Level 0. In addition,
//    - push_back() and resize() automatically calls reserve() if required.
template<class T>
class fast_vector<T, 1u> : public fast_vector<T, 0u>
{
public:
    using fast_vector<T, 0u>::reserve;
    using fast_vector<T, 0u>::_size;
    using fast_vector<T, 0u>::size;
    using fast_vector<T, 0u>::_capacity;
    using fast_vector<T, 0u>::capacity;
    using fast_vector<T, 0u>::data;
    //using fast_vector<T, 0u>::operator[];

    inline fast_vector() : fast_vector<T,0u>() { reserve(4); }
    inline fast_vector(size_t n) : fast_vector<T,0u>(n) { }
    fast_vector(size_t n, T const& item) : fast_vector<T,0u>(n, item) { } //_capacity(n), old_capacity(0u), _size(n) 
    //    { data = (T*)malloc(sizeof(T)*n); T* _end=data+_size; for (T* dp=data; dp!=_end; ++dp) *dp=item; }
    
    void push_back(T const& item) { if(_size==_capacity) reserve(_capacity+MIN(_capacity,1048576u)); data[_size++] = item; }
    void resize(size_t n) { if(_capacity<n) reserve(MAX(n,_capacity+MIN(_capacity,1048576u))); _size = n; }
    inline void resize(size_t n, T const& item) 
        { if(_capacity<n) reserve(MAX(n,_capacity+MIN(_capacity,1048576u))); for (int a=_size; a<n; ++a) data[a]=item;  _size = n; }
};

/* // Level 2:
// All features of Level 1. In addition,
//      - operator[] and at() checks for bounds, and calls reserve() if required.
template<class T>
class fast_vector<T, 2u> : public fast_vector<T, 1u>
{
public:
    
}; */

/* // ================================================================================================
// Efficient vector with find() and remove() operations

template<class T, unsigned int level=0u>
class fast_search_vector : public fast_vector<T,level>
{
public:
    
}; */

#endif
