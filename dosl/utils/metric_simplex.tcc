#ifndef __DOSL_METRIC_SIMPLEX_TCC__ 
#define __DOSL_METRIC_SIMPLEX_TCC__ 

// constants for comparing values
#define _MS_DOUBLE_EPS     1e-4    // make smaller if grid edge length is less than O(1)
#define _MS_DOUBLE_EPS_SQ  (_MS_DOUBLE_EPS*_MS_DOUBLE_EPS)

// includes
#include "macros_constants.tcc"
// --
#include "metric_simplex_bits/metric_simplex_base.tcc"
#include "metric_simplex_bits/metric_simplex_a.tcc"
#include "metric_simplex_bits/metric_simplex_collection.tcc"

#endif
