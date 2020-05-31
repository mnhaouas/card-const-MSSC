/*
 * Author: Haouas, Mohammed Najib - oct 29th, 2018
 *                 > Last modified: may 30th, 2020
 *
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


#include "IlcWCSS.h"


IloConstraint IloWCSS(IloEnv env, IloIntVarArray X, IloFloatVar V, const Data* data, const char* name = 0);