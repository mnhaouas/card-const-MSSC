/*
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


#include "IlcWCSS_NetworkCardControl.h"


IlcWCSS_NetworkCardControlI::IlcWCSS_NetworkCardControlI(IloCPEngine cp, IlcIntVarArray X, IlcFloatVar V, const Data& data) :
IlcConstraintI(cp), _X(X), _V(V), _dissimilarities(data.dissimilarities), _n(X.getSize()), _k(data.K) {
    // Populate _targetCards
    _targetCards = IlcIntArray(cp, _k);
    for (int c = 0; c < data.K; c++)
        _targetCards[c] = data.targetCardinalities[c];

    // Make sure everything's cute
    assert(_targetCards.getSize() == _k); // number of indicated cards is same as clusters
    IlcInt ctrl_nb_pts = 0;
    for (IlcInt c = 0; c < _k; c++)
        ctrl_nb_pts += _targetCards[c];
    assert(ctrl_nb_pts == _n); // number of indicated points is same as problem size

    // sets of points and their sizes
    setP_assigned = new (cp.getHeap()) std::vector<IlcInt>[_k]; // setP_assigned[c] = i means point i is assigned to cluster c

    // size of each cluster c for easy access
    sizeCluster = IlcIntArray(cp, _k);

    // Sum of dissimilarities of each cluster c
    S1 = IlcFloatArray(cp, _k); // S1[c] = sum of dissimilarities (squared) of cluster c 

    // Sum of dissimilarities (squared) between each unassigned point and each cluster
    s2 = new (cp.getHeap()) IlcFloat*[_n]; // Allocation on engine heap for efficient and unified memory control
    for (IlcInt i = 0; i < _n; i++)
        s2[i] = new (cp.getHeap()) IlcFloat[_k]; // s2[x][c] = sum of dissimilarities between x and all points in cluster c

    // Smallest 1/2 contribution of each unassigned point together with m other points
    s3 = new (cp.getHeap()) std::vector<IlcFloat>[_n]; // s3[x][m] = smallest contribution of the m other free points brought along in addition to x

    // Map between main problem variables and CPLEX model variables
    problem_to_cplex_var_map = new (cp.getHeap()) IlcInt*[_n];
    for (IlcInt i = 0; i < _n; i++)
        problem_to_cplex_var_map[i] = new (cp.getHeap()) IlcInt[_k]; // problem_to_cplex_var_map[i][c] = var number in CPLEX model representing arc between point i in U and cluster c

    // Map between main problem cluster numbers and CPLEX model variables
    problem_to_cplex_cluster_var_map = IlcIntArray(cp, _k); // problem_to_cplex_cluster_var_map[c] = var number in CPLEX model representing arc between cluster c and drain

    // Backup for which edges have flow after CPLEX is run
    hasFlow = new (cp.getHeap()) IlcRevBool*[_n];
    for (IlcInt i = 0; i < _n; i++)
        hasFlow[i] = new (cp.getHeap()) IlcRevBool[_k]; // indices are equivalent to problem_to_cplex_var_map

    // Global lower bound as computed by CPLEX
    lb_global = new (cp.getHeap()) IlcRevFloat(cp, 0.0);

    // Points to add in cluster c for completion
    nb_points_to_add = IlcIntArray(cp, _k);

    // This espsilon is subtracted from the computed lower bound to prevent false backtracking while comparing with upper bound due to rounding errors
    _epsc = 5e-3;

    // Keep a record of destination flow for each point in dataset prior to the current propagation (useful to know if there is a need to update MCF)
    destination = new (cp.getHeap()) IlcRevInt[_n];
    for (IlcInt i = 0; i < _n; i++)
        destination[i].setValue(cp, -1); // destination[i] = cluster to which the flow coming from i goes

    // Keep track of what variables were already fixed prior to the current propagation (useful to know if there is a need to update MCF)
    varWasFixed = new (cp.getHeap()) IlcRevBool[_n];
    for (IlcInt i = 0; i < _n; i++)
        varWasFixed[i].setValue(cp, IlcFalse);
}


IlcWCSS_NetworkCardControlI::~IlcWCSS_NetworkCardControlI() {
    // Any dynamically allocated elements are allocated in the engine heap which manages memory for us
}


void IlcWCSS_NetworkCardControlI::post() {
    for (IlcInt i = 0; i < _n; i++)
        _X[i].whenDomain(this);

    _V.whenRange(this);
}


void IlcWCSS_NetworkCardControlI::propagate() {
    /*
     * Preliminaries: propagation process relies on essential assumptions.
     *     This part, among other things, enforces those assumptions and handles special cases when they occur
     *     NOTE: The astute reader may notice that much of this part could be computed incrementally.
     *           In practice however, the overhead associated with tracking changes to operate incremental computations
     *           cancels any speedup they would produce.
     */
        
        // Reset sets & number of points assigned and unassigned
        setU_unassigned.clear(); q = 0;
        for (IlcInt c = 0; c < _k; c++)
            setP_assigned[c].clear();
        p = 0;

        // Populating sets and definition of partial problem
        for (IlcInt i = 0; i < _n; i++) {
            if (_X[i].isFixed()) {
                p++;
                setP_assigned[_X[i].getValue()].push_back(i);
            }
            else {
                q++;
                setU_unassigned.push_back(i);
            }
        }

        // Size of each cluster c and how many points to add
        max_clust_completion = 0; // Max number of points to add, all clusters considered
        for (IlcInt c = 0; c < _k; c++) {
            sizeCluster[c] = setP_assigned[c].size(); // sizeCluster[c] = m means c is size m
            nb_points_to_add[c] = _targetCards[c] - sizeCluster[c];

            if (nb_points_to_add[c] < 0) // Constraint is violated if a cluster is overfilled
                fail(); // Backtrack

            if (nb_points_to_add[c] > max_clust_completion)
                max_clust_completion = nb_points_to_add[c]; // Update max_clust_completion
        }

        // Preliminary filtering: remove possibility to assign points to cluster c if full
        // NOTE: This can be replaced with a GCC in the model only if there is a guarantee it would propagate before this constraint.
        bool prelimFilteringVarWasFixed;
        do {
            // Nothing has yet changed entering this loop
            prelimFilteringVarWasFixed = false;

            // Filter values if corresponding clusters are filled
            for (IlcInt c = 0; c < _k; c++) {
                if (nb_points_to_add[c] == 0) {
                    std::vector<IlcInt>::iterator setU_iter = setU_unassigned.begin();

                    while (setU_iter != setU_unassigned.end()) {
                        if (_X[*setU_iter].isInDomain(c)) {
                            _X[*setU_iter].removeValue(c); // Careful, a variable could get bound here, hence the whole reason this "do... while" exists in the first place

                            if (_X[*setU_iter].isFixed()) {
                                prelimFilteringVarWasFixed = true;
                                setP_assigned[_X[*setU_iter].getValue()].push_back(*setU_iter); p++;
                                setU_iter = setU_unassigned.erase(setU_iter); q--; 
                            } else {
                                ++setU_iter;
                            }
                        } else {
                            ++setU_iter;
                        }
                    }
                }
            }

            // If some variable was fixed, update sets and subproblem characteristics
            if (prelimFilteringVarWasFixed) {
                // Size of each cluster c and how many points to add
                max_clust_completion = 0; // max number of points to add all clusters considered
                for (IlcInt c = 0; c < _k; c++) {
                    sizeCluster[c] = setP_assigned[c].size(); // sizeCluster[c] = m means c is size m
                    nb_points_to_add[c] = _targetCards[c] - sizeCluster[c];

                    if (nb_points_to_add[c] < 0) // Constraint is violated if a cluster is overfilled
                        fail();

                    if (nb_points_to_add[c] > max_clust_completion)
                        max_clust_completion = nb_points_to_add[c]; // Update max_clust_completion
                }
            }
        } while (prelimFilteringVarWasFixed);

        // If no points are assigned, which can happen when posting this constraint, there is no work to be done
        if (q == _n) {
            _X[0].setValue(0); // The first point must be assigned to 0 in all cases when symmetry-breaking constraints are active
            return;
        }


    /*
     * Kitchen: prepare ingredients for further steps, ie relevant point contributions to each cluster
     */

        // Sum of dissimilarities of each cluster c
        for (IlcInt c = 0; c < _k; c++) {
            S1[c] = 0;
            for (IlcInt i = 0; i < (sizeCluster[c] - 1); i++)
                for (IlcInt j = i + 1; j < sizeCluster[c]; j++)
                    S1[c] += _dissimilarities[setP_assigned[c][i]][setP_assigned[c][j]];
        }

        // Variable mapping for clusters between CP realm and CPLEX realm
        cluster_not_filled_counter = 0;
        for (IlcInt c = 0; c < _k; c++) {
            if (nb_points_to_add[c] == 0) { // If cluster filled
                problem_to_cplex_cluster_var_map[c] = -1; // Exclude from MCF model
            }
            else {
                problem_to_cplex_cluster_var_map[c] = cluster_not_filled_counter;
                cluster_not_filled_counter++;
            }
        }

        // Initialize var counter for observation variables in CPLEX realm 
        cplex_var_map_counter = 0;

        // Sum of dissimilarities between each unassigned point and each cluster & variable mapping for observations between CP realm and CPLEX realm
        for (IlcInt i = 0; i < q; i++) { // for each unassigned point
            for (IlcInt c = 0; c < _k; c++) { // for each cluster
                if (_X[setU_unassigned[i]].isInDomain(c) && problem_to_cplex_cluster_var_map[c] != -1) {
                    // Second condition in this "if" statement is redundant
                    //     if nb_points_to_add[c] == 0, then necessarily _X[setU_unassigned[i]].isInDomain(c) is false
                    problem_to_cplex_var_map[i][c] = cplex_var_map_counter;
                    cplex_var_map_counter++;

                    s2[i][c] = 0;
                    for (IlcInt j = 0; j < sizeCluster[c]; j++) { // for each point j in cluster c
                        s2[i][c] += _dissimilarities[setU_unassigned[i]][setP_assigned[c][j]];
                    }
                }
                else { // else, unassigned point i can't be part of cluster c, set to infinity and exclude
                    problem_to_cplex_var_map[i][c] = -1;
                    s2[i][c] = IlcInfinity; // For completeness; it is not needed.
                }
            }
        }

        // Smallest 1/2 contributions of each unassigned point together with m other points
        for (IlcInt i = 0; i < q; i++) {
            s3[i].clear(); // Important to clear vector from previous run
            for (int j = 0; j < q; j++) {
                // At first, we put all points in that list
                s3[i].push_back(_dissimilarities[setU_unassigned[i]][setU_unassigned[j]] / 2);
            }

            std::sort(s3[i].begin(), s3[i].end()); // sort distances for each point, first element 0 because d(x,x) = 0

            for (int j = 1; j < max_clust_completion; j++)
                s3[i][j] += s3[i][j-1]; // Compute minimum 1/2 contributions, first element is 0
        }

        // Check whether meaningful change has occured that warrants fresh computations
        bool activeVarValHasChanged = false;

        for (IlcInt i = 0; i < _n; i++) {
            // This condition is true is no MCF has been solved to date
            //     In this case, it is necessary to make fresh computations anyways
            if (destination[i].getValue() == -1) {
                activeVarValHasChanged = true;
                break;
            }

            // If a point was assigned to a cluster other than one to which it was matched by the last MCF
            if (_X[i].isFixed() && _X[i].getValue() != destination[i].getValue()) {
                activeVarValHasChanged = true;
                break;
            }

            // If the cluster to which a point was matched by the last MCF no longer exists in the domain of the corresponding var
            if (!_X[i].isInDomain(destination[i].getValue())) {
                activeVarValHasChanged = true;
                break;
            }
        }

        // If a var was newly bound after the last propagation
        for (IlcInt i = 0; i < _n; i++) {
            if (_X[i].isFixed() && !varWasFixed[i].getValue()) {
                activeVarValHasChanged = true;
                varWasFixed[i].setValue(getCPEngine(), IlcTrue); // Update record for subsequent propagations
            }
        }

        // If there's an arc in the last MCF solution that can no longer be occupied
        for (IlcInt i = 0; i < q; i++) {
            if (destination[setU_unassigned[i]].getValue() == -1 || problem_to_cplex_var_map[i][destination[setU_unassigned[i]].getValue()] == -1 || problem_to_cplex_cluster_var_map[destination[setU_unassigned[i]].getValue()] == -1) {
                activeVarValHasChanged = true;
                break;
            }
        }


    /*
     * Core: Minimum-Cost Flow (MCF) formulation for lower-bound computations
     *     Resolution using CPLEX Optimizer through Concert Technology
     *     WARNING: Concert Technology environment and underlying objects are NOT automatically destroyed
     *              when they go out of scope, resulting in a memory leak. Effort must be made for appropriate
     *              cleanup beyond this point when leaving this scope (in particular when failures occur).
     */

        if (activeVarValHasChanged) { // A meaningful change has occured that warrants fresh computations 
            // CPLEX environment for minimum-cost flow
            IloEnv cpx_env;

            // Model
            IloModel cpx_model(cpx_env);

            // Variables
            IloInt _cn = q + cplex_var_map_counter + cluster_not_filled_counter; // Number of vars in CPLEX model
            IloNumVarArray _cX(cpx_env, _cn, 0.0, max_clust_completion); // CPLEX model variables, minimum arc capacity is 0. Most arcs should have capacity 1,
                                                                         // except arcs linking clusters to drain, capacity in that case should be nb_points_to_add[c]

            // NETWORK: primary source
            IloNumExpr _csource(cpx_env);

            for (IloInt i = 0; i < q; i++) {
                _csource += _cX[i];
                cpx_model.add(_cX[i] <= 1); // Max capacity is 1, see above
            }

            cpx_model.add(_csource == q);

            // NETWORK: transit unbound points
            for (int i = 0; i < q; i++) {
                IloNumExpr _cTempTransitPoint(cpx_env);
                _cTempTransitPoint -= _cX[i];

                for (IlcInt c = 0; c < _k; c++) {
                    if (problem_to_cplex_var_map[i][c] != -1) {
                        _cTempTransitPoint += _cX[q + problem_to_cplex_var_map[i][c]];
                        cpx_model.add(_cX[q + problem_to_cplex_var_map[i][c]] <= 1); // Max capacity is 1, see above
                    }
                }

                cpx_model.add(_cTempTransitPoint == 0);
            }

            // NETWORK: transit clusters
            IloNumExprArray incoming_cluster_flow_cost(cpx_env, _k); // Convenient way to get access to actual induced cost per cluster
            for (IloInt c = 0; c < _k; c++) {
                incoming_cluster_flow_cost[c] = IloNumExpr(cpx_env); // Empty handles are now null expressions
                incoming_cluster_flow_cost[c] += 0;
            }

            for (int c = 0; c < _k; c++) {
                if (problem_to_cplex_cluster_var_map[c] != -1) {
                    IloNumExpr _cTempTransitCluster(cpx_env);

                    bool ctrlAssignment = false;
                    for (int i = 0; i < q; i++) {
                        if (problem_to_cplex_var_map[i][c] != -1) {
                            ctrlAssignment = true;
                            _cTempTransitCluster -= _cX[q + problem_to_cplex_var_map[i][c]];
                            incoming_cluster_flow_cost[c] += ((s2[i][c] + s3[i][nb_points_to_add[c] - 1])*_cX[q + problem_to_cplex_var_map[i][c]]);
                        }
                    }

                    if (!ctrlAssignment) {
                        cpx_env.end(); // Cleanup
                        fail(); // If no incoming arcs to a cluster that must house observations, then this branch is unsuccessful
                    }

                    _cTempTransitCluster += _cX[q + cplex_var_map_counter + problem_to_cplex_cluster_var_map[c]];
                    cpx_model.add(_cTempTransitCluster == 0);

                    cpx_model.add(_cX[q + cplex_var_map_counter + problem_to_cplex_cluster_var_map[c]] <= nb_points_to_add[c]); // max capacity is nb_points_to_add[c]
                }
            }

            IloNumExpr lb_global_expr(cpx_env); // Convenient way to get access to global lower bound
            for (IlcInt c = 0; c < _k; c++) {
                lb_global_expr += ((S1[c] + incoming_cluster_flow_cost[c]) / _targetCards[c]);
            }

            // NETWORK: global drain
            IloNumExpr _cdrain(cpx_env);

            for (IlcInt c = 0; c < _k; c++)
                if (problem_to_cplex_cluster_var_map[c] != -1)
                    _cdrain -= _cX[q + cplex_var_map_counter + problem_to_cplex_cluster_var_map[c]];

            cpx_model.add(_cdrain == -q);

            // Objective
            cpx_model.add(IloMinimize(cpx_env, lb_global_expr));

            // Set up CPLEX engine
            IloCplex cplex(cpx_model);
            //cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Network); // ATTENTION: Uncomment to force usage of Network Simplex
                                                                                 //     In practice, it is much better to let CPLEX decide which alg to use
            cplex.setOut(cpx_env.getNullStream()); // Redirect IloAlgorithm cplex's output to the environment's null stream for quiet operation

            // Solve
            if (!cplex.solve()) {
                cpx_env.end();
                fail(); // If CPLEX can't solve model, it means this branch can't be successful because there is no valid assignment of free points
            }

            lb_global->setValue(getCPEngine(), cplex.getObjValue() - _epsc); // Protection against rounding errors, see constructor

            // Backup which arcs have flow, useful for further value filtering
            for (IlcInt i = 0; i < q; i++)
                for (IlcInt c = 0; c < _k; c++)
                    // Solution is integral, using inequality to shield against rounding errors
                    hasFlow[i][c].setValue(getCPEngine(), (cplex.getValue(_cX[q + problem_to_cplex_var_map[i][c]]) > 0.5));

            // No need for CPLEX beyond this point, release memory
            cpx_env.end();
        }


    /*
     * Filtering: cost-based filtering based on assignment assumptions
     *     Reusing most recent valid MCF solution for efficient computation
     */
        
        _V.setMin(lb_global->getValue()); // lb_global (ie lb_global_expr) and _V.getMax() have slightly different values (rounding errors). 
                                          // This means, sometimes, failure occurs when near a new, improving solution even though it shouldn't.
                                          // Removing a small epsilon solves the problem.
                                          // In an abundance of caution, apply as large an epsilon as possible. In this case, higher precision is superfluous.
                                          // Seriously, rounding errors are the devil.

        // Update records for next propagation if new MCF has been calculated
        if (activeVarValHasChanged)
            for (IlcInt i = 0; i < q; i++)
                for (IlcInt c = 0; c < _k; c++)
                    if (problem_to_cplex_var_map[i][c] != -1 && hasFlow[i][c].getValue())
                        destination[setU_unassigned[i]].setValue(getCPEngine(), c);

        if (activeVarValHasChanged)
            for (IlcInt i = 0; i < _n; i++)
                if (_X[i].isFixed())
                    destination[i].setValue(getCPEngine(), _X[i].getValue());

        // Variable fitering
        for (IlcInt c = 0; c < _k; c++) { // for each value c in domains of points, ie for each cluster
            for (IlcInt i = 0; i < q; i++) {
                if (problem_to_cplex_var_map[i][c] != -1 && !hasFlow[i][c].getValue()) {
                    // Get objective value increase if i-th point in U is assigned to cluster c
                    deltaObj = getDeltaObj(i, destination[setU_unassigned[i]].getValue(), c); // -1 means infeasible updated solution
                    
                    // If new objective exceeds incumbent cost
                    if (deltaObj < -0.1 || (lb_global->getValue() + deltaObj) > _V.getMax()) { // -0.1 is shield against rounding errors.
                        if (_X[setU_unassigned[i]].getSize() == 1) {
                            // For some reason, CP optimizer, in extremely rare cases, would NOT fail if the dom of a var is emptied
                            //     We force failure here
                            fail();
                        }

                        _X[setU_unassigned[i]].removeValue(c);
                    }
                }
            }
        }


    // Ended propagate
}


inline IlcFloat IlcWCSS_NetworkCardControlI::getDeltaObj(IlcInt origin_i, IlcInt origin_c, IlcInt targeted_c) {
    //                                                          ^ relocated point                 ^ headed to
    //                                                                           ^ whence it came

    // Bellman-Ford
    // Rarely, there could be negative-weight cycles in this graph (which is why the version of this alg with a queue fails; inf loop). That's not a problem.
    //     The path we're interested in (from targeted_c to origin_c) can never have negative cycles and so the true weight will always be returned.
    //     This is because if that were true, we can make its weight arbitrarily small, thus lowering the starting objective value (which can't happen; delta must be positive).
    //
    //     Here |V| = q + _k - 1 because we remove the vertex corresponding to origin_i and any edges that lead to/emanate from it.
    //     Graph here is bipartite, on the left are vertices representing points in U. On the right, partial clusters.
    //     Vertices are named 0..q-1 on the left-hand side for the points, q + c# on the right-hand side for the clusters.
    std::vector<IlcFloat> graphMinDist(q + _k, IlcInfinity);
    graphMinDist[q + targeted_c] = 0; // Origin is targeted_c, where we have excess flow to redirect

    bool hasChangedWeights;
    for (int pass = 1; pass <= (q + _k - 2); pass++) {
        hasChangedWeights = false;

        for (IlcInt i = 0; i < q; i++) {
            for (IlcInt c = 0; c < _k; c++) {
                if (i != origin_i && c != targeted_c && problem_to_cplex_var_map[i][c] != -1 && !hasFlow[i][c].getValue()) {
                    // Going right
                    // c != targeted_c: can never return to originating node. If there is a lower-weight path from targeted_c, it means negative-weight cycle.
                    if ((graphMinDist[i] + (s2[i][c] + s3[i][nb_points_to_add[c] - 1]) / _targetCards[c]) < graphMinDist[q + c]) {
                        graphMinDist[q + c] = graphMinDist[i] + (s2[i][c] + s3[i][nb_points_to_add[c] - 1]) / _targetCards[c];
                        hasChangedWeights = true;
                    }
                } else if (i != origin_i && c != origin_c && problem_to_cplex_var_map[i][c] != -1 && hasFlow[i][c].getValue()) {
                    // Going left
                    // c != origin_c: can never leave destination node. If there is a lower-weight path from origin_c, it means negative-weight cycle.
                    if ((graphMinDist[q + c] - (s2[i][c] + s3[i][nb_points_to_add[c] - 1]) / _targetCards[c]) < graphMinDist[i]) {
                        graphMinDist[i] = graphMinDist[q + c] - (s2[i][c] + s3[i][nb_points_to_add[c] - 1]) / _targetCards[c];
                        hasChangedWeights = true;
                    }
                }
            }
        }

        // If no change has occurred, no change will ever occur, stop
        if (!hasChangedWeights)
            break;
    }

    // Unreachable destination (ie, infeasible flow)
    if (!(graphMinDist[q + origin_c] < IlcInfinity)) // Destination is origin_c, where we have flow deficit
        return -1;

    // Delta obj is equal to delta contribution of origin_i (could be negative) + min-weight path weight of sending excess flow to destination (could be negative).
    //     However, delta obj can never be negative.
    return (s2[origin_i][targeted_c] + s3[origin_i][nb_points_to_add[targeted_c] - 1]) / _targetCards[targeted_c]
         - (s2[origin_i][origin_c] + s3[origin_i][nb_points_to_add[origin_c] - 1]) / _targetCards[origin_c]
         + graphMinDist[q + origin_c];
}


// Function which returns an engine handle (IlcConstraintI*) for the constraint implementation IlcWCSS_NetworkCardControlI
IlcConstraint IlcWCSS_NetworkCardControl(IlcIntVarArray X, IlcFloatVar V, const Data& data) {
    IlcCPEngine cp = X.getCPEngine(); // Get CP engine from variable array
    return new (cp.getHeap()) IlcWCSS_NetworkCardControlI(cp, X, V, data);
}


// Macro which wraps the engine constraint handle into a modeling layer (Concert Technology) extractable object
//     For some reason, ILOCPCONSTRAINTWRAPPERn doesn't work with reference types; have to use pointer to _datao
ILOCPCONSTRAINTWRAPPER3(IloWCSS_NetworkCardControl, cp, IloIntVarArray, _Xo, IloFloatVar, _Vo, const Data*, _datao) {
    use(cp, _Xo); // Force extraction of modeling layer extractables (ie, get engine level objects)
    use(cp, _Vo);
    return IlcWCSS_NetworkCardControl(cp.getIntVarArray(_Xo), cp.getFloatVar(_Vo), *_datao);
}