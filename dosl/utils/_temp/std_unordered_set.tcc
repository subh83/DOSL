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

#ifndef __DOSL_STD_UNORDERED_SET_TCC
#define __DOSL_STD_UNORDERED_SET_TCC

#include <unordered_set>


template <class key, // 'key' must have copy constructer defined
          class HashFunctor = std::hash<key>,
          class EqualToFunctor = std::equal_to<key> >
class std_unordered_set : public std::unordered_set<key,HashFunctor,EqualToFunctor>
{
public:
    using std::unordered_set<key,HashFunctor,EqualToFunctor>::begin;
    using std::unordered_set<key,HashFunctor,EqualToFunctor>::end;
    using std::unordered_set<key,HashFunctor,EqualToFunctor>::size;
    using std::unordered_set<key,HashFunctor,EqualToFunctor>::empty;
    using std::unordered_set<key,HashFunctor,EqualToFunctor>::clear;
    using std::unordered_set<key,HashFunctor,EqualToFunctor>::insert;
    
    std_unordered_set () { }
    std_unordered_set (size_t n=1024 , const HashFunctor& hf = HashFunctor(), const EqualToFunctor& eql = EqualToFunctor())
         : std::unordered_set<key,HashFunctor,EqualToFunctor>(n, hf, eql) {  }
     
    key* get (key& n) {
        //auto it = find(n);
        //if (it != end()) return &(*it);
        return &(*(insert(n).first));
    }
};

#endif


