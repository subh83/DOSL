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

#ifndef _PLANNER_BITS_HPP
#define _PLANNER_BITS_HPP

#define declare_alg_name(s)    static constexpr const char* AlgorithmName = s; \
                               std::string algorithm_name(void) { return (s); }


#define FUNCTOR_BEGIN(functor_name,data_class_name)     class functor_name { \
                                                        public: \
                                                            data_class_name* data_p; \
                                                            functor_name () : data_p(NULL) { } \
                                                            functor_name (data_class_name* p) : data_p(p) { }
                                                            
                                                            // RTYPE operator()(...) { return _this->fun(...); }
                                                            
#define FUNCTOR_END(instance_name)                      } instance_name;

#endif
