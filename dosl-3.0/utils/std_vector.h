#ifndef __DOSL_STD_VECTOR_H_
#define __DOSL_STD_VECTOR_H_

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
