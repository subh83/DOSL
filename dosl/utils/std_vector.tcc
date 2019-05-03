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

#ifndef __DOSL_STD_VECTOR_TCC_
#define __DOSL_STD_VECTOR_TCC_

#include<vector>

// ================================================================================================
// Wrapper on std::vector
template<class T, unsigned int level=1u> class std_vector;

template<class T>
class std_vector<T, 1u> : public std::vector<T>
{
public:
    typedef typename std::vector<T>::iterator iterator;
    using std::vector<T>::begin;
    using std::vector<T>::end;
    using std::vector<T>::size;
    using std::vector<T>::operator[];
    
    // Functions not native to std::vector
    template <class Pred=std::equal_to<T> >
    inline iterator find(T const& item) 
        { Pred _pred; for (iterator it=begin(); it!=end(); ++it) { if (_pred(*it,item)) return (it); } return (end()); }
    template <class Pred=std::equal_to<T> >
    inline iterator rfind(T const& item) 
        { Pred _pred; iterator it=end(); do { --it; if (_pred(*it,item)) return (it); } while (it!=begin()); return (end()); }
    
    template <class Pred=std::equal_to<T> >
    inline int findi(T const& item) 
        { Pred _pred; for (int a=0; a<size(); ++a) { if (_pred(operator[](a),item)) return (a); } return (-1); }
    template <class Pred=std::equal_to<T> >
    inline int rfindi(T const& item) 
        { Pred _pred; for (int a=size()-1; a>=0; --a) { if (_pred(operator[](a),item)) return (a); } return (-1); }
    
    inline void erase (size_t n)
        { erase(begin()+n); }
};

#endif
