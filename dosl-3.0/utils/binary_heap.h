#ifndef __DOSL_BINARY_HEAP_H
#define __DOSL_BINARY_HEAP_H

#include <stdio.h>
#include <cstdlib>
#include <limits>
#include <unordered_map>
#include <functional> // std::less

#include "utils/macros_constants.h"
#include "utils/misc_wrappers.h"

// ==========================================================================================

template <class key>
class DefautHeapPosFunctor
{
public:
    std::unordered_map <key,int>  keyBoolMap;
    int& operator() (key const & k) {
        if (keyBoolMap.find(k)==keyBoolMap.end())
            keyBoolMap[k] = -1;
        return(keyBoolMap[k]);
    }
};

// -----------------

template <class T>
class simple_less {
public:
    bool operator() (const T& x, const T& y) {return x<y;}
};

// -----------------

#define _less_than(k1,k2)  ( (HeapPosFunctorInstance(k2)<0)  ||  CompareFunctorInstance((k1),(k2)) ) 
#define _in_heap(k)  (HeapPosFunctorInstance(k)>=0)

template <class key,  // use pointer here directly
          class CompareFunctor = simple_less<key>,
          class HeapPosFunctor = DefautHeapPosFunctor<key> // Provides "int& HeapPosFunctor::operator() (key&)"
          /*, class ComparatorClass = less<key>*/ >
class binary_heap
{
private:
    // cmpare function
    typedef bool (CompareFunctor::*CompareFunctionPointerType)(key const &, key const &);
    ExplicitFunctor<CompareFunctor, CompareFunctionPointerType>  CompareFunctorInstance;
    
    // heappos function
    typedef int& (HeapPosFunctor::*HeapPosFunctionPointerType)(key const &); 
                                        // ^^^ TODO: 'const': what if we need to change key?
    ExplicitFunctor<HeapPosFunctor, HeapPosFunctionPointerType>  HeapPosFunctorInstance;
    
public:
    typedef key  Key;
    
    // Variables for storing heap
    _DOSL_FAST_VECTOR <Key> heapArray;
    
    typedef enum {UNKNOWN_BUBDIR, UP_BUBDIR, DOWN_BUBDIR} BubbleDirection;
    
    // set the functions
    void set_compare_function (CompareFunctor* compare_functor_instance_pointer=NULL, 
                                    CompareFunctionPointerType compare_function_pointer=NULL)
        { CompareFunctorInstance.set_pointers (compare_functor_instance_pointer, compare_function_pointer); }
    
    void set_heappos_function (HeapPosFunctor* heappos_functor_instance_pointer=NULL, 
                                    HeapPosFunctionPointerType heappos_function_pointer=NULL)
        { HeapPosFunctorInstance.set_pointers (heappos_functor_instance_pointer, heappos_function_pointer); }
    
    
    // Main interface functions
    // push
    void insert (key const & np);
    void push (key const & np) { insert(np); }
    // TODO: void push (Key n);
    // erase
    void erase (key const & np);
    void remove (key const & np) { erase(np); }
    // update, pop
    Key pop (void);
    void update (Key np, BubbleDirection bubdir=UNKNOWN_BUBDIR);
    
    inline void clear (void) { while (heapArray.size()>0) pop(); }
    inline bool empty (void) { return (heapArray.size()==0); }
    inline int size (void) { return (heapArray.size()); }
    #if _DOSL_DEBUG
    bool HeapCheck (int heapNodePos=0, int depth=0);
    #endif
    
    // destructor
    // ~binary_heap () { delete (ComparatorClassInstancePointer); }
};

// =====================================================================================

template <class key, class CompareFunctor, class HeapPosFunctor>
void binary_heap<key,CompareFunctor,HeapPosFunctor>::update (Key np, BubbleDirection bubdir)
{
    // Does both bubble up and down
    if (!(_in_heap(np))) { _dosl_err("attempted to heap-update node not in heap"); } //return;
    int currentPos = HeapPosFunctorInstance(np);
    
    // Bubble up
    if (currentPos>0 && bubdir!=DOWN_BUBDIR) { // either up or unknown
        int parentPos = (currentPos - 1) / 2; // position of parent
        if ( _less_than (np, heapArray[parentPos]) ) { // parent is larger. bubble up
            heapArray[currentPos] = heapArray[parentPos];
            heapArray[parentPos] = np;
            HeapPosFunctorInstance (heapArray[currentPos]) = currentPos;
            HeapPosFunctorInstance (heapArray[parentPos]) = parentPos;
            update(heapArray[parentPos], UP_BUBDIR); // recursive call. If heap, downward bubbling won't be needed.
            return; // no need to check children, since chilren will have higher values
        }
        else if (bubdir == UP_BUBDIR) // no need to check for bubble down if it was given that it would be up
            return;
    }
    
    // Bubble down
    int leftChildPos = 2*currentPos + 1;
    int rightChildPos = 2*currentPos + 2;
    int exchangeChildPos;
    if (leftChildPos<size() && rightChildPos<size())
        exchangeChildPos = ( _less_than (heapArray[rightChildPos], heapArray[leftChildPos]) ) ? rightChildPos : leftChildPos;
    else if (leftChildPos<size())
        exchangeChildPos = leftChildPos;
    else if (rightChildPos<size())
        exchangeChildPos = rightChildPos;
    else
        return; // no child exists
        
    if ( _less_than (heapArray[exchangeChildPos], np) ) { // the smallest child is smaller than self. bubble down
        heapArray[currentPos] = heapArray[exchangeChildPos];
        heapArray[exchangeChildPos] = np;
        HeapPosFunctorInstance (heapArray[currentPos]) = currentPos;
        HeapPosFunctorInstance (heapArray[exchangeChildPos]) = exchangeChildPos;
        update(heapArray[exchangeChildPos], DOWN_BUBDIR); // recursive call. If heap, upward bubbling won't be needed.
        return;
    }
}

// SB_BINHEAP -------------------------------------------------------


template <class key, class CompareFunctor, class HeapPosFunctor>
void binary_heap<key,CompareFunctor,HeapPosFunctor>::insert (key const & np)
{
    // Assumes np is not in heap. No check is done for that.
    // Add as leaf
    // Check if already exists. Don't push then.
    if (_in_heap(np)) return;
    
    heapArray.push_back(np);
    HeapPosFunctorInstance(np) = size() - 1;
    // Bubble up
    update(np, UP_BUBDIR);
}


template <class key, class CompareFunctor, class HeapPosFunctor>
key binary_heap<key,CompareFunctor,HeapPosFunctor>::pop (void)
{
    #if _DOSL_AUTOCORRECT
    if (heapArray.size()==0)  return (NULL);
    #endif
    Key ret = heapArray[0];
    HeapPosFunctorInstance(ret) = -1; // not in heap
    // Take last leaf and place it at root
    heapArray[0] = heapArray[size()-1];
    HeapPosFunctorInstance (heapArray[0]) = 0;
    heapArray.pop_back();
    // Bubble down
    if (heapArray.size()>0)
        update(heapArray[0], DOWN_BUBDIR);
    
    return (ret);
}


template <class key, class CompareFunctor, class HeapPosFunctor>
void binary_heap<key,CompareFunctor,HeapPosFunctor>::erase (key const & np)
{
    #if _DOSL_AUTOCORRECT
    if (!(_in_heap(np)))  return;
    #endif
    HeapPosFunctorInstance(np) = -1;
    update(np);
    pop();
}

// ================================

#if _DOSL_DEBUG
template <class key, class CompareFunctor, class HeapPosFunctor>
bool binary_heap<key,CompareFunctor,HeapPosFunctor>::HeapCheck (int heapNodePos, int depth)
{
    if (heapNodePos >= size()) return (true);
    
    #if _DOSL_DEBUG > 2
    for (int a=0; a<depth; ++a) std::cout << " ";
    std::cout << heapArray[heapNodePos] << std::endl;
    #endif
    
    int leftChildPos = 2*heapNodePos + 1;
    int rightChildPos = 2*heapNodePos + 2;
    if ( (leftChildPos < size()  &&  !_less_than(heapArray[heapNodePos], heapArray[leftChildPos]) )  ||
            (rightChildPos < size()  &&  !_less_than (heapArray[heapNodePos], heapArray[rightChildPos]) ) ) {
        #if _DOSL_DEBUG > 2
        //for (int a=0; a<size(); ++a)
            //printf("%d: %f, ", a, heapArray[a]);
        #endif
        std::cout << "heapArray[" << heapNodePos << "]=" << heapArray[heapNodePos] << ", heapArray[" << leftChildPos << "]=" << heapArray[leftChildPos] << ", heapArray[" << rightChildPos << "]=" << heapArray[rightChildPos] << "\n"; 
        _dosl_err("Heap property violated!!");
        return (false);
    }
    return ( HeapCheck(leftChildPos, depth+1) && HeapCheck(rightChildPos, depth+1) );
}
#endif

#endif
