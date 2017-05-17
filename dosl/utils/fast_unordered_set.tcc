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
#ifndef __DOSL_FAST_UNORDERED_SET_TCC
#define __DOSL_FAST_UNORDERED_SET_TCC

#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <functional>

#include "macros_constants.tcc"
#include "misc_wrappers.tcc"

#ifndef _DOSL_HASH_BIN_CHECK
#define _DOSL_HASH_BIN_CHECK 1
#endif

// ==========================================================================================

template <class T>
class simple_equal_to {
public:
    bool operator() (const T& x, const T& y) {return x==y;}
};

// -----------------

template <class key, // 'key' must have copy constructer defined
          class HashFunctor = std::hash<key>,
          class EqualToFunctor = simple_equal_to<key> >
class fast_unordered_set
{
/* 'Key' must have "bool operator==(key const &)" overloaded. */
/* 'HashFunctor' should be a class with "int operator()(key const &)" overloaded */
private:
    // hash function
    typedef unsigned int (HashFunctor::*HashFunctionPointerType)(key &);
    ExplicitFunctor<HashFunctor, HashFunctionPointerType>  HashFunctorInstance;
    
    // equal-to function
    typedef bool (EqualToFunctor::*EqualToFunctionPointerType)(key const &, key const &);
    ExplicitFunctor<EqualToFunctor, EqualToFunctionPointerType>  EqualToFunctorInstance;

public:
    typedef key Key;
    
    // parameters (TODO: make private)
    int HashTableSize;
    
    // Functions for setting private variables
    
    void set_hash_table_size (int s) 
        { if (HashTable) _dosl_err("Hash table already initiated!"); HashTableSize = s; }
    
    void set_hash_function (HashFunctor* hash_functor_instance_pointer=NULL, 
                                    HashFunctionPointerType hash_function_pointer=NULL) {
        if (HashTable) _dosl_err("Hash table already initiated!");
        HashFunctorInstance.set_pointers (hash_functor_instance_pointer, hash_function_pointer);
    }
    
    void set_equal_to_function (EqualToFunctor* equalto_functor_instance_pointer=NULL, 
                                        EqualToFunctionPointerType equalto_function_pointer=NULL) {
        if (HashTable) _dosl_err("Hash table already initiated!");
        EqualToFunctorInstance.set_pointers (equalto_functor_instance_pointer, equalto_function_pointer);
    }
    
    // variables
    _DOSL_FAST_VECTOR <Key*>* HashTable; // contsins pointers to the items
    int Size;
    
    // initializor
    void init (void);
    
    // Constructors
    fast_unordered_set() : HashTable(NULL), Size(-1), HashTableSize(1024) { }
    
    // Main interface function
    inline bool empty (void) { return (HashTable==NULL); }
    inline int size() { return (Size); }
    Key* get (key & n); // Returns pointer to already-existing item, else creates one
    
    // Clear
    void clear (bool destroyHashTable=false) {
        if (HashTable) {
            for (int a=0; a<HashTableSize; a++)
                for (int b=0; b<HashTable[a].size(); b++)
                    delete HashTable[a][b];
            Size = 0;
            if (destroyHashTable) {
                delete[] HashTable;
                HashTable=NULL;
                Size=-1;
            }
        }
    }
    
    // Destructor
    ~fast_unordered_set() {
        if (HashTable) {
            for (int a=0; a<HashTableSize; a++)
                for (int b=0; b<HashTable[a].size(); b++)
                    delete HashTable[a][b];
            delete[] HashTable;
            HashTable = NULL;
        }
    }
};

// ==========================================================================================


template <class key, class HashFunctor, class EqualToFunctor>
void fast_unordered_set<key,HashFunctor,EqualToFunctor>::init (void)
{
    if (HashTable) _dosl_err("Hash table already initiated! Call clear to reset and re-initiate.");
    
    HashTable = new _DOSL_FAST_VECTOR <Key*> [ HashTableSize ];
    Size = 0;
}

template <class key, class HashFunctor, class EqualToFunctor>
key* fast_unordered_set<key,HashFunctor,EqualToFunctor>::get (key & n)
{
    if (!HashTable) init();
    
    // Search in bin
    unsigned int hashBin = HashFunctorInstance(n);
    #if _DOSL_HASH_BIN_CHECK
    hashBin %= HashTableSize;
    #endif
    
    // A SIGSEGV signal generated from here most likely 'GetHashBin' returned a bin index larger than (HashTableSize-1).
    for (int a=0; a<HashTable[hashBin].size(); ++a)
        if ( EqualToFunctorInstance (*(HashTable[hashBin][a]), n) )
            return (HashTable[hashBin][a]);
    
    // If new node, create it!
    Key* newHashItem_p = new Key (n); // Invokes copy constructor of Key
    HashTable[hashBin].push_back (newHashItem_p);
    ++Size;
    return (newHashItem_p);
}

#endif
