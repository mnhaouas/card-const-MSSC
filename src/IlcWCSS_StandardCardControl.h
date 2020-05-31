/*
 * Author: Haouas, Mohammed Najib - mar 29th, 2019
 *                 > Last modified: may 30th, 2020
 *
 * This constraint should be used in conjunction with an appropriate model to speedup exact Minimum Sum of Squares Clustering with pre-set, strict cluster cardinalities.
 * Therefore, this constraint can be used for the special case of balanced Minimum Sum of Squares Clustering (bMSSC).
 * This constraint should provide a significant boost over using a general constraint (refer to IloWCSS) in conjunction with a GCC to constrain the cluster cardinalities.
 * This is because prior knowledge of target cardinalities enables the calculation of stronger bounds that can lead to further domain reductions.
 * The objective function is filtered through tightening of the lower bound.
 * The reprentative variables are filtered through cost-based filtering against an upper bound decided by the best incumbent solution.
 *
 * Main arguments: * X, array of integer representative variables that link observations to their cluster.
 *                 * V, WCSS of solution, must be constrained to take value of WCSS in Concert Technology model.
 *
 * Additional arguments: * data, refer to Data struct in Data.h for problem data nomenclature.
 *
 * Note: this constraint is also heavily dependent on the search strategy. Use IloMSSCSearchStrategy (IloGoal).
 *
 * This constraint is a modification of the work of:
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

// Problem data structure
#include "Data.h"

// Using (IBM ILOG CPLEX Optimization Studio) CP Optimizer Extensions
#include <ilcp/cpext.h>


// Refer to cpp file for explanation of variables defined here


ILOSTLBEGIN


class IlcWCSS_StandardCardControlI : public IlcConstraintI {
protected:
    double const* const* const _dissimilarities;
    IlcIntArray _targetCards;

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

    IlcFloat lb_global;

    IlcFloat lb_except;
    IlcFloat lb_prime;

    // Added
    IlcIntArray nb_points_to_add;
    IlcInt max_clust_completion;
    IlcFloat V_prime;

    double _epsc;

public:
    IlcWCSS_StandardCardControlI(IloCPEngine cp, IlcIntVarArray X, IlcFloatVar V, const Data& data);
    ~IlcWCSS_StandardCardControlI();
    virtual void propagate();
    virtual void post();
    IlcIntVarArray getEngineVars() { return _X; }
};

IlcConstraint IlcWCSS_StandardCardControl(IlcIntVarArray X, IlcFloatVar V, const Data& data);