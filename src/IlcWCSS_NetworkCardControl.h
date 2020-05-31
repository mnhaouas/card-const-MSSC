/*
 * Author: Haouas, Mohammed Najib - sept 5th, 2019
 *                 > Last modified: may 30th, 2020
 *
 * This constraint should be used in conjunction with an appropriate model to speedup exact Minimum Sum of Squares Clustering with pre-set, strict cluster cardinalities.
 * Therefore, this constraint can be used for the special case of balanced Minimum Sum of Squares Clustering (bMSSC).
 * This constraint should provide a significant boost over using IloWCSS_StandardCardControl because of more thorough bounds computed by means of MCF resolution.
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
 * This constraint uses elements from the work of:
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

// Using (IBM ILOG CPLEX Optimization Studio) CPLEX Optimizer through Concert Technology
#include <ilcplex/ilocplex.h>


// Refer to cpp file for explanation of variables defined here


ILOSTLBEGIN


class IlcWCSS_NetworkCardControlI : public IlcConstraintI {
protected:
    inline IlcFloat getDeltaObj(IlcInt origin_i, IlcInt origin_c, IlcInt target_c);
    IlcFloat deltaObj;

    IlcRevInt* destination;
    IlcRevBool* varWasFixed;

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

    IlcFloatArray S1;
    IlcFloat** s2;
    std::vector<IlcFloat>* s3;

    IlcInt cplex_var_map_counter;
    IlcInt** problem_to_cplex_var_map;
    IlcIntArray problem_to_cplex_cluster_var_map;

    IlcRevFloat* lb_global;

    IlcIntArray nb_points_to_add;
    IlcInt max_clust_completion;
    IlcInt cluster_not_filled_counter;
    IlcRevBool** hasFlow;

    double _epsc;

public:
    IlcWCSS_NetworkCardControlI(IloCPEngine cp, IlcIntVarArray X, IlcFloatVar V, const Data& data);
    ~IlcWCSS_NetworkCardControlI();
    virtual void propagate();
    virtual void post();
    IlcIntVarArray getEngineVars() { return _X; }
};

IlcConstraint IlcWCSS_NetworkCardControl(IlcIntVarArray X, IlcFloatVar V, const Data& data);