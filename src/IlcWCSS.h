/*
 * This constraint should be used in conjunction with an appropriate model to speedup exact Minimum Sum of Squares Clustering (MSSC).
 * The speedup is the result of simultaneous filtering of both the objective variable V (total Within Cluster Sum of Squares, WCSS) and the reprentative variables X of the observations.
 * The objective function is filtered through tightening of the lower bound.
 * The reprentative variables are filtered through cost-based filtering against an upper bound decided by the best solution to date.
 *
 * Main arguments: * X, array of integer representative variables that link observations to their cluster.
 *                 * V, WCSS of solution, must be constrained to take value of WCSS in Concert Technology model.
 *
 * Additional arguments: * data, refer to Data struct in Data.h for problem data nomenclature.
 *
 * Note: this constraint is also heavily dependent on the search strategy. Use IloMSSCSearchStrategy (IloGoal).
 *
 * Implementation of the work of:
 * Dao TBH., Duong KC., Vrain C. (2015) Constrained Minimum Sum of Squares Clustering by Constraint Programming.
 *     In: Pesant G. (eds) Principles and Practice of Constraint Programming. CP 2015.
 *     Lecture Notes in Computer Science, vol 9255. Springer, Cham
 *     doi:10.1007/978-3-319-23219-5_39
 *
 * We are grateful to Dao et al. for supplying us with the code of their own implementation using Gecode.
 *     This code is largely based on theirs.
 *
 * Note: Dao et al. have come up with a more efficient CP framework to solve general MSSC based on Repetitive Branch and Bound Alg (CP RBBA).
 *       This is not an implementation of that work.
 */


// Vector and vector operations
#include <algorithm>
#include <vector>

// Problem data structure
#include "Data.h"

// Using (IBM ILOG CPLEX Optimization Studio) CP Optimizer Extensions
#include <ilcp/cpext.h>


// Refer to cpp file for explanation of variables defined here


ILOSTLBEGIN


class IlcWCSSI : public IlcConstraintI {
protected:
    double const* const* const _dissimilarities;

    IlcInt _n, _k; // size of problem, nb of clusters
    IlcInt p, q; // size of sets P and U resp.

    IlcIntVarArray _X; // Point assignments
    IlcFloatVar _V; // total WCSS

    // In propagate
    std::vector<IlcInt> setU_unassigned;
    std::vector<IlcInt>* setP_assigned;

    IlcIntArray sizeCluster;

    IlcFloat** lb_schedule;

    IlcFloatArray S1;
    IlcFloat** s2;
    std::vector<IlcFloat>* s3;

    IlcFloat** lb_global;

    IlcFloatArray lb_except;
    IlcFloatArray lb_prime;

    double _epsc;

public:
    IlcWCSSI(IloCPEngine cp, IlcIntVarArray X, IlcFloatVar V, const Data& data);
    ~IlcWCSSI();
    virtual void propagate();
    virtual void post();
    IlcIntVarArray getEngineVars() { return _X; }
};


IlcConstraint IlcWCSS(IlcIntVarArray X, IlcFloatVar V, const Data& data);