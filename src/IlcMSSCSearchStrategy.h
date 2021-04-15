/*
 * This is the branching strategy that guides the CP search. It is implemented as a goal to pass to the CP engine.
 * This goal places on the goals stack subsequent goals at each branching decision until a full solution has been instantiated.
 * Goals are generated in such a fashion as to produce a binary branching.
 * This search strategy has 3 operation modes: initial solution generation, subsequent search and tie-handling (when they occur).
 *
 * Main arguments: * vars, branching variables. In this context, representative integer variables of observations.
 *
 * Additional arguments: * data, refer to Data struct in Data.h for problem data nomenclature.
 *                       * searchParameters, see below for information on CustomCPSearchOptions
 *                       * solFound, bool from scope where engine IloCP is instantiated.
 *                             Is set to true once a first solution has been found using the engine's IloCP::next method.
 *
 * This search strategy uses elements from the work of:
 * Dao TBH., Duong KC., Vrain C. (2015) Constrained Minimum Sum of Squares Clustering by Constraint Programming.
 *     In: Pesant G. (eds) Principles and Practice of Constraint Programming. CP 2015.
 *     Lecture Notes in Computer Science, vol 9255. Springer, Cham
 *     doi:10.1007/978-3-319-23219-5_39
 *
 * This is part of my (Haouas, M.N.) MSc's research project under supervision of Pesant, G. & Aloise, D.
 *
 * Work pending publication in the CPAIOR 2020 conference proceeding as a full paper (submission 70).
 *     https://cpaior2020.dbai.tuwien.ac.at/papers/
 */


// Vector and vector operations
#include <algorithm>
#include <vector>

// Using (IBM ILOG CPLEX Optimization Studio) CP Optimizer Extensions
#include <ilcp/cpext.h>

// Problem data structure
#include "Data.h"

// Possible to use <limits> but eh...
#define __MAX_INT 2147483647


ILOSTLBEGIN


#ifndef __SEARCH_T
#define __SEARCH_T

namespace CustomCPSearchOptions {
    enum class InitialSolution {
        NONE, // Let main search generate an initial solution
        GREEDY_INIT, // Generate init sol through a greedy assignment that minimizes delta objective at each step
        MEMBERSHIPS_AS_INDICATED // Generate init sol as instructed in Data passed
    };

    enum class MainSearch {
        MAX_MIN_VAR // Branch on variable which induces the largest minimum delta objective
    };

    enum class TieHandling {
        NONE, // No tie handling
        UNBOUND_FARTHEST_TOTAL_SS, // Start empty cluster farthest to un-assigned points on the basis of total sum of squares
        FIXED_FARTHEST_DIST, // Start empty cluster farthest to fixed points on the basis of distance alone
        FIXED_MAX_MIN, // Start empty cluster at the point that has the maximum distance to its closest cluster
        FARTHEST_POINT_FROM_BIGGEST_CENTER, // Start empty cluster at the point that is farthest to the biggest cluster center
        MAX_MIN_POINT_FROM_ALL_CENTER // Start empty cluster at the point that has maximum minimum distance to all cluster centers
    };
}


struct SearchParameters {
    CustomCPSearchOptions::InitialSolution initialSolution;
    CustomCPSearchOptions::MainSearch mainSearch;
    CustomCPSearchOptions::TieHandling tieHandling;
};

#endif // !__SEARCH_T




int getDeltaObjective(IlcIntVarArray vars, IlcInt pt, IlcInt c, double const* const* const dissimilarities);
int getUnboundPointsTotalSS(IlcIntVarArray vars, IlcInt pt, double const* const* const dissimilarities);
int getIntDist(IlcInt i, IlcInt j, double const* const* const dissimilarities);