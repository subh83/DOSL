#ifndef _PLANNER_BITS_HPP
#define _PLANNER_BITS_HPP

#define declare_alg_name(s)    static constexpr const char* AlgorithmName = s; \
                               std::string algorithm_name(void) { return (s); }


#endif
