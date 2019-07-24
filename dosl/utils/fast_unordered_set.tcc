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

#ifndef __DOSL_FAST_UNORDERED_SET_TCC
#define __DOSL_FAST_UNORDERED_SET_TCC

#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <functional>

#include "macros_constants.tcc"
//#include "misc_wrappers.tcc"

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
    HashFunctor hash_functor_instance; // hash function
    EqualToFunctor equal_to_functor_instance; // equal-to function
    size_t hash_table_size;
    
public:
    typedef key Key;
    
    // Functions for setting private variables
    
    void reserve (size_t s) 
        { if (HashTable) _dosl_err("Hash table already initiated!"); hash_table_size = s; }
    
    // variables
    _DOSL_FAST_VECTOR <Key*>* HashTable; // contsins pointers to the items
    int Size;
    
    // initializor
    void init (void);
    
    // Constructors
    fast_unordered_set () : HashTable(NULL), Size(-1), hash_table_size(1024) { }
    fast_unordered_set (size_t n=1024 , const HashFunctor& hf = HashFunctor(), const EqualToFunctor& eql = EqualToFunctor())
         : HashTable(NULL), Size(-1), hash_table_size(n), hash_functor_instance(hf), equal_to_functor_instance(eql) {  }
    
    // Main interface function
    inline bool empty (void) { return (HashTable==NULL); }
    inline int size() { return (Size); }
    Key* get (key & n); // Returns pointer to already-existing item, else creates one
    template<class FindEqFunc=EqualToFunctor> std::vector<Key*> getall (key & n, FindEqFunc eqfunc_instance, bool search_other_bins=false);
    
    // Clear
    void clear (bool destroyHashTable=false) {
        if (HashTable) {
            for (int a=0; a<hash_table_size; a++)
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
            for (int a=0; a<hash_table_size; a++)
                for (int b=0; b<HashTable[a].size(); b++)
                    delete HashTable[a][b];
            delete[] HashTable;
            HashTable = NULL;
        }
    }
    
    // ----------------------------------
    
    class forward_iterator {
    public:
        fast_unordered_set* p;
        int bin, pos;
        forward_iterator (fast_unordered_set* pp, int bb=0, int poss=-1) : p(pp), bin(bb), pos(poss) { if (pos<0) operator++(); }
        
        // increment (++fit)
        forward_iterator& operator++() {
            do {
                if (bin >= p->hash_table_size) { bin = p->hash_table_size; pos=0; } // end
                else if (pos < p->HashTable[bin].size() - 1) ++pos;
                else { ++bin; pos=-1; }
            } while (pos<0);
            return *this;
        }
        // compare
        bool operator==(const forward_iterator& other) { return (bin==other.bin && pos==other.pos); }
        bool operator!=(const forward_iterator& other) { return (bin!=other.bin || pos!=other.pos); }
        // pointer
        key& operator*() const { return *(p->HashTable[bin][pos]); }
        key*& operator->() const { return (p->HashTable[bin][pos]); }
    };
    
    forward_iterator begin(void) { return forward_iterator(this); }
    forward_iterator end(void) { return forward_iterator(this, hash_table_size); }
    
};

// ==========================================================================================


template <class key, class HashFunctor, class EqualToFunctor>
void fast_unordered_set<key,HashFunctor,EqualToFunctor>::init (void)
{
    if (HashTable) _dosl_err("Hash table already initiated! Call clear to reset and re-initiate.");
    
    HashTable = new _DOSL_FAST_VECTOR <Key*> [ hash_table_size ];
    Size = 0;
}

template <class key, class HashFunctor, class EqualToFunctor>
key* fast_unordered_set<key,HashFunctor,EqualToFunctor>::get (key & n)
{
    if (!HashTable) init(); // TODO: remove this check?
    
    // Search in bin
    unsigned int hashBin = hash_functor_instance(n);
    #if _DOSL_HASH_BIN_CHECK
    hashBin %= hash_table_size;
    #endif
    
    // A SIGSEGV signal generated from here most likely 'GetHashBin' returned a bin index larger than (hash_table_size-1).
    for (int a=0; a<HashTable[hashBin].size(); ++a)
        if ( equal_to_functor_instance (*(HashTable[hashBin][a]), n) )
            return (HashTable[hashBin][a]);
    
    // If new node, create it!
    Key* newHashItem_p = new Key (n); // Invokes copy constructor of Key
    HashTable[hashBin].push_back (newHashItem_p);
    ++Size;
    return (newHashItem_p);
}

template <class key, class HashFunctor, class EqualToFunctor>
template<class FindEqFunc>
std::vector<key*> fast_unordered_set<key,HashFunctor,EqualToFunctor>::getall (key& n, FindEqFunc eqfunc_instance, bool search_other_bins) {
    if (!HashTable) init(); // TODO: remove this check?
    
    unsigned int startHashBin=0, endHashBinP1=hash_table_size;
    // Search in bin
    if (!search_other_bins) {
    startHashBin = hash_functor_instance(n);
    #if _DOSL_HASH_BIN_CHECK
    startHashBin %= hash_table_size;
    #endif
    endHashBinP1 = startHashBin+1;
    }
    
    std::vector<key*> ret;
    for (int hashBin=startHashBin; hashBin<endHashBinP1; ++hashBin)
        for (int a=0; a<HashTable[hashBin].size(); ++a)
            if ( eqfunc_instance (*(HashTable[hashBin][a]), n) )
                ret.push_back (HashTable[hashBin][a]);
                
    return (ret);
}

#endif
