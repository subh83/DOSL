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

#ifndef __DOSL_SAFE_SHARED_PTR_TCC
#define __DOSL_SAFE_SHARED_PTR_TCC

#include <mutex>
#include <memory>

/*template <class T>
class safe_shared_ptr_data {
public:
    std::shared_ptr<T> data;
    std::mutex mtx;
};*/

template <class T>
class safe_shared_ptr : public std::shared_ptr<T> {
public:
    
    std::mutex* mtx_p;
    
    // constructors
    safe_shared_ptr() : std::shared_ptr<T>() { mtx_p = new std::mutex; };
    /*constexpr shared_ptr(nullptr_t) : shared_ptr() {}
    template <class U> explicit shared_ptr (U* p);
    template <class U, class D> shared_ptr (U* p, D del);
    template <class D> shared_ptr (nullptr_t p, D del);
    template <class U, class D, class Alloc> shared_ptr (U* p, D del, Alloc alloc);
    template <class D, class Alloc> shared_ptr (nullptr_t p, D del, Alloc alloc);
    shared_ptr (const shared_ptr& x) noexcept;
    template <class U> shared_ptr (const shared_ptr<U>& x) noexcept;
    template <class U> explicit shared_ptr (const weak_ptr<U>& x);
    shared_ptr (shared_ptr&& x) noexcept;
    template <class U> shared_ptr (shared_ptr<U>&& x) noexcept;
    template <class U> shared_ptr (auto_ptr<U>&& x);
    template <class U, class D> shared_ptr (unique_ptr<U,D>&& x);
    template <class U> shared_ptr (const shared_ptr<U>& x, element_type* p) noexcept; */

    
    
    // Access
    // "behave like pointer"
    T* operator->() { return (instance_p); }
    T& operator*() { return (*instance_p); }
    // convert to pointer of type T*
    //operator T*() const { return (instance_p); }
    
};

#endif
