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

#ifndef __DOSL_THREAD_UTILS_TCC
#define __DOSL_THREAD_UTILS_TCC

#include <thread>
#include <mutex>

#ifndef _DOSL_MULTITHREAD
#define _DOSL_MULTITHREAD 1
#endif

#if _DOSL_MULTITHREAD
    #define CREATE_MUTEX(name)  std::mutex _ ## name ## _mutex;
    #define LOCK_MUTEX(name)    _ ## name ## _mutex.lock();
    #define UNLOCK_MUTEX(name)  _ ## name ## _mutex.unlock();
#else
    #define CREATE_MUTEX(name)  
    #define LOCK_MUTEX(name)    
    #define UNLOCK_MUTEX(name)  
#endif

#define MUTEXED_MEMBER(type,name)   type _ ## name; \
                                    CREATE_MUTEX(name) \
                                    void set_ ## name (type& v) { \
                                        LOCK_MUTEX(name) \
                                        _ ## name = v; \
                                        UNLOCK_MUTEX(name) \
                                    }
/* Use:

*/

#endif
