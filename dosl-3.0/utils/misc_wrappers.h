#ifndef __DOSL_MISC_WRAPPERS_H
#define __DOSL_MISC_WRAPPERS_H

#include <type_traits>
#include <utility>
#include<stdexcept>
#include <typeinfo>
#include <sstream>
#include <string>

template <class T>
class AutoPointer {
// AutoPointer<T> -- maintains a pointer of type T*, including creating and destroying it if necessary.

private:
    T* instance_p;
    bool new_created;

public:
    
    AutoPointer () : new_created(false) { }
    AutoPointer (T* in_pointer) : new_created(false) { init (in_pointer); }
    
    // explicit init
    
    void init (T* in_pointer=NULL) { 
        if (in_pointer) {
            clear(); // clears current pointer
            instance_p = in_pointer;
            new_created = false;
        }
        else if (!instance_p) { // init() -- will initiate if already not initiated.
            instance_p = new T;
            new_created = true;
        }
    }
    
    void operator() (T* in_pointer) { init(in_pointer); }
    
    // Delete / destructor
    
    void clear (void) {
        if (new_created)
            delete instance_p;
        instance_p = NULL;
        new_created = false;
    }
    
    ~AutoPointer () {
        clear();
    }
    
    // Access
    
    T* operator->() { init(); return (instance_p); }
    
    T& operator*() { init(); return (*instance_p); }
};

// ====================================================
// Simplified macro-based auto-pointer

#define COPY_IF_NOTNULL_ELSE_CREATE_POINTER_TO_LOCAL(type,inPointer,copyPointer,defaultVal) \
                type _localvalue_ ## copyPointer; \
                type* copyPointer = & _localvalue_ ## copyPointer; \
                if (inPointer) copyPointer = inPointer; \
                else *copyPointer = defaultVal;

// ====================================================

template <class F>
struct return_type;

template <class R, class C, class... A>
struct return_type <R (C::*)(A...)> {
  typedef R type;
};

template <class R, class C, class... A>
struct return_type <R (C::*)(A...) const> {
    typedef R type;             // ^^^^^
};


// SFINAE for getting operator()
template <typename fun_ptr_type, fun_ptr_type fun_ptr> struct type_check;

template <typename functor_type, typename fun_ptr_type>
    fun_ptr_type  get_bracket_operator (type_check<fun_ptr_type, &functor_type::operator()>*)
        { return ( static_cast<fun_ptr_type>(&functor_type::operator()) ); }
        
/* template <typename functor_type, typename fun_ptr_type>
    fun_ptr_type  get_bracket_operator (type_check<fun_ptr_type const, &functor_type::operator()>*)
        { return ( static_cast<fun_ptr_type>(&functor_type::operator()) ); } */
        
template <typename functor_type, typename fun_ptr_type>
    fun_ptr_type  get_bracket_operator (...) { return (NULL); }

// --------------

template <class FunctorType, typename FunctionPointerType> // FunctionPointerType is a pointer type to a member function of F
class ExplicitFunctor {
    /* Declaration Syntax: ExplicitFunctor <functor_class, ... (functor_class::*)(...)>  ef;
       Use: ef (...)
    */
private:
    FunctorType* functor_instance_p;
    bool new_created;
    FunctionPointerType function_pointer;
    
public:

    // function_pointer instantiation: function_pointer = &F::f
    // ------------------------------------------
    
    ExplicitFunctor () : functor_instance_p(NULL), new_created(false), function_pointer(NULL) { /*init();*/ }
    ExplicitFunctor (FunctorType* in_pointer, FunctionPointerType in_func=NULL) : 
                        functor_instance_p(NULL), new_created(false), function_pointer(NULL)  { set_pointers (in_pointer, in_func); }
    
    // check if empty
    bool empty (void) { return (functor_instance_p==NULL); }
    
    // ------------
    // explicit initializations
    
    void set_functor_instance (FunctorType* in_pointer=NULL) {
        if (in_pointer) {
            clear(); // clears current pointer
            functor_instance_p = in_pointer;
            new_created = false;
        }
        else if (!functor_instance_p) { // init() -- will initiate if already not initiated.
            functor_instance_p = new FunctorType;
            new_created = true;
        }
    }
    
    void set_function_pointer (FunctionPointerType in_func=NULL) {
        if (in_func)
            function_pointer = in_func;
        else
            function_pointer = get_bracket_operator<FunctorType,FunctionPointerType>(0);
    }
    
    void set_pointers (FunctorType* in_pointer=NULL, FunctionPointerType in_func=NULL) {
        set_functor_instance (in_pointer);
        set_function_pointer (in_func);
    }
    
    // ------------
    // operator()
    
    /*template <class F>
    struct return_type;

    template <class R, class... A>
    struct return_type <R (FunctorType::*)(A...)> {
      typedef R type;
    };*/
    
    // references
    /* template <class... A>
    typename return_type<FunctionPointerType>::type  operator() (A & ... args) { // receive everything by reference
        if (!functor_instance_p) set_functor_instance();
        if (!function_pointer)  set_function_pointer();
        return ( (functor_instance_p->*function_pointer) (args...) );
    }
    
    // const references
    template <class... A>
    typename return_type<FunctionPointerType>::type  operator() (const A & ... args) { // receive everything by reference
        if (!functor_instance_p) set_functor_instance();
        if (!function_pointer)  set_function_pointer();
        return ( (functor_instance_p->*function_pointer) (args...) );
    } */
    
    // forwarded
    template <class... A>
    typename return_type<FunctionPointerType>::type // std::result_of<FunctionPointerType>::type // return_type<FunctionPointerType>::type 
            operator() (A && ... args) { // receive everything by reference
        if (!functor_instance_p) set_functor_instance();
        if (!function_pointer) {
            set_function_pointer();
            if (!function_pointer) {
                std::cout << std::flush;
                std::stringstream err_str;
                err_str << "no operator() in class '" << typeid(FunctorType).name() <<
                                    "' matching signature '" << typeid(FunctionPointerType).name() << "'";
                throw std::runtime_error(err_str.str().c_str());
            }
        }
        return ( (functor_instance_p->*function_pointer) (std::forward<A>(args)...) );
    }
    
    // ------------
    
    // Delete / destructor
    
    void clear (void) {
        if (new_created)
            delete functor_instance_p;
        functor_instance_p = NULL;
        new_created = false;
    }
    
    ~ExplicitFunctor () {
        clear();
    }
};




// ----------------------

template <class FunctorType, typename FunctionPointerType>
ExplicitFunctor<FunctorType,FunctionPointerType>
    CreateExplicitFunctor (FunctorType* functor_instance_p, FunctionPointerType function_pointer)
        { return (ExplicitFunctor<FunctorType,FunctionPointerType> (functor_instance_p, function_pointer)); }

// ====================================================

// To allocate memory when only pointer_type is known.
template <class T>
T* _new_p (T* dummy) {
    return (new T);
}
#define new_p(pointer_type)  _new_p(static_cast<pointer_type>(NULL))


#endif
