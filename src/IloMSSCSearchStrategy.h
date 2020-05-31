/*
 * Author: Haouas, Mohammed Najib - sept 10th, 2019
 *                 > Last modified: may 30th, 2020
 *
 * This is the branching strategy that guides the CP search. It is implemented as a goal to pass to the CP engine.
 * This goal places on the goals stack subsequent goals at each branching decision until a full solution has been instantiated.
 * Goals are generated in such a fashion as to produce a binary branching.
 * This search strategy has 3 operation modes: initial solution generation, subsequent search and tie-handling (when they occur).
 *
 * Main arguments: * vars, branching variables. In this context, representative integer variables of observations.
 *
 * Additional arguments: * data, refer to Data struct in Data.h for problem data nomenclature.
 *                       * searchParameters, refer to IlcMSSCSearchStrategy.h for information on CustomCPSearchOptions
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


#include "IlcMSSCSearchStrategy.h"


IloGoal IloMSSCSearchStrategy(IloEnv env, IloIntVarArray vars, const Data& data, const SearchParameters& searchParameters, const bool& solFound);