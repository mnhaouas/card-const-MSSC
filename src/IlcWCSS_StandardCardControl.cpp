/*
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


#include "IlcWCSS_StandardCardControl.h"


IlcWCSS_StandardCardControlI::IlcWCSS_StandardCardControlI(IloCPEngine cp, IlcIntVarArray X, IlcFloatVar V, const Data& data) : 
IlcConstraintI(cp), _X(X), _V(V), _dissimilarities(data.dissimilarities), _n(X.getSize()), _k(data.K) {
    // Populate _targetCards
    _targetCards = IlcIntArray(cp, _k);
    for (int c = 0; c < data.K; c++)
        _targetCards[c] = data.targetCardinalities[c];

    // Make sure everything's cute
    assert(_targetCards.getSize() == _k); // number of indicated cards is same as clusters

    IlcInt ctrl_nb_pts = 0;
    for (int c = 0; c < _k; c++)
        ctrl_nb_pts += _targetCards[c];
    assert(ctrl_nb_pts == _n); // number of indicated points is same as problem size

    // sets of points and their sizes
    setP_assigned = new (cp.getHeap()) std::vector<IlcInt>[_k]; // setP_assigned[c] = i means point i is assigned to cluster c

    // size of each cluster c
    sizeCluster = IlcIntArray(cp, _k);

    // lower bound of WCSS of each cluster c
    lb_schedule = new (cp.getHeap()) IlcFloat*[_k];
    for (int i = 0; i < _k; i++)
        lb_schedule[i] = new (cp.getHeap()) IlcFloat[2]; // lb_schedule[c][m] is lower bound on WCSS if c completed target (cardinality - m).

    // Sum of dissimilarities of each cluster c
    S1 = IlcFloatArray(cp, _k); // S1[c] = sum of dissimilarities (squared) of cluster c 

    // Sum of dissimilarities (squared) between each unassigned point and each cluster
    s2 = new (cp.getHeap()) IlcFloat*[_n]; // Allocation on engine heap for efficient and unified memory control
    for (int i = 0; i < _n; i++)
        s2[i] = new (cp.getHeap()) IlcFloat[_k]; // s2[x][c] = sum of dissimilarities between x and all points in cluster c

    // Smallest 1/2 contribution of each unassigned point together with m other points
    s3 = new (cp.getHeap()) std::vector<IlcFloat>[_n]; // s3[x][m] = smallest contribution of the m other free points brought along in addition to x

    // Global lower bound
    lb_global = 0;

    // Lower bound of all clusters except active cluster under assumption to house active point in filtering
    lb_except = 0;

    // Lower bound of active cluster under assumption to house active point in filtering
    lb_prime = 0;

    // Number of points to add in cluster c to completion
    nb_points_to_add = IlcIntArray(cp, _k);

    // This espsilon is subtracted from the computed lower bound to prevent false backtracking while comparing with upper bound due to rounding errors
    _epsc = 5e-5;
}


IlcWCSS_StandardCardControlI::~IlcWCSS_StandardCardControlI() {
    // Any dynamically allocated elements are allocated in the engine heap which manages memory for us
}


void IlcWCSS_StandardCardControlI::post() {
    for (int i = 0; i < _n; i++)
        _X[i].whenDomain(this);
    
    _V.whenRange(this);
}


void IlcWCSS_StandardCardControlI::propagate() {
    // Reset sets & number of points assigned and unassigned
    setU_unassigned.clear(); q = 0;
    for (int c = 0; c < _k; c++)
        setP_assigned[c].clear();
    p = 0;

    // Populating sets and definition of partial problem
    for (int i = 0; i < _n; i++) {
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

    // lower bound of WCSS of each cluster c if we add m points to it > INIT
    for (int c = 0; c < _k; c++)
        for (int m = 0; m < 2; m++)
            lb_schedule[c][m] = 0; // lb_schedule[c][m] is lower bound on WCSS if c completed target (cardinality - m).

    // Sum of dissimilarities of each cluster c
    for (int c = 0; c < _k; c++) {
        S1[c] = 0;
        for (int i = 0; i < (sizeCluster[c] - 1); i++)
            for (int j = i + 1; j < sizeCluster[c]; j++)
                S1[c] += _dissimilarities[setP_assigned[c][i]][setP_assigned[c][j]];
    }

    // Sum of dissimilarities between each unassigned point and each cluster
    for (int i = 0; i < q; i++) { // for each unassigned point
        for (int c = 0; c < _k; c++) { // for each cluster
            if (nb_points_to_add[c] > 0 && _X[setU_unassigned[i]].isInDomain(c)) { // if unassigned point i can be assigned to cluster c
                s2[i][c] = 0;
                for (int j = 0; j < sizeCluster[c]; j++) { // for each point j in cluster c
                    s2[i][c] += _dissimilarities[setU_unassigned[i]][setP_assigned[c][j]];
                }
            }
            else { // else, unassigned point i can't be part of cluster c, set to infinity to exclude
                s2[i][c] = IlcInfinity;
            }
        }
    }

    // Smallest 1/2 contribution of each unassigned point together with m other points
    //     We only need to study adding max_clust_completion
    for (int i = 0; i < q; i++) {
        s3[i].clear();
        for (int j = 0; j < q; j++) {
            // At first, we put all points in that list
            s3[i].push_back(_dissimilarities[setU_unassigned[i]][setU_unassigned[j]] / 2);
        }

        std::sort(s3[i].begin(), s3[i].end()); // sort distances for each point, first element 0 because d(x,x) = 0

        for (int j = 1; j < max_clust_completion; j++)
            s3[i][j] += s3[i][j - 1]; // Compute minimum 1/2 contributions, first element is 0
    }

    // Computing lower bound for each cluster
    //     We only need to study adding discreet amounts of points
    for (int c = 0; c < _k; c++) { // for each cluster
        for (int m = 0; m < 2; m++) { // if we add (nb_points_to_add[c] - m) points to c
            std::vector<IlcFloat> s; // s of every point x with current m possibilities
            for (int i = 0; i < q; i++) { // contribution for each unassigned point, must consider them all for sorting
                if ((nb_points_to_add[c] - m) > 0) { // if we have points to add
                    s.push_back(s2[i][c] + s3[i][(nb_points_to_add[c]) - 1]);
                }
                else {
                    s.push_back(0); // not adding anything, this is actually useless
                }
            }

            std::sort(s.begin(), s.end()); // sorting candidate points...

            IlcFloat S2 = 0;
            for (int i = 0; i < (nb_points_to_add[c] - m); i++) // ... of which we select the (nb_points_to_add[c] - m) points that induce the lowest cost
                S2 += s[i];

            lb_schedule[c][m] = (S1[c] + S2) / (nb_points_to_add[c] + sizeCluster[c] - m); // -m because at some point we take one fewer point (max m = 1).
        }
    }

    // Global lower bound
    //     No need for dynamic programming
    lb_global = 0;
    for (int c = 0; c < _k; c++)
        lb_global += lb_schedule[c][0];

    // Filter objective
    _V.setMin(lb_global - _epsc); // lb_global and _V.getMax() have slightly different values (rounding errors). 
                                  // This means, sometimes, failure occurs when near a new, improving solution even though it shouldn't.
                                  // Removing a small epsilon solves the problem.
                                  // In an abundance of caution, apply as large an epsilon as possible. In this case, higher precision is superfluous.
                                  // Seriously, rounding errors are the devil.

    for (int c = 0; c < _k; c++) { // for each value c in domains of points, ie for each cluster
        lb_except = lb_global - lb_schedule[c][0]; // Remove contribution of cluster c

        for (int i = 0; i < q; i++) {
            if (_X[setU_unassigned[i]].isInDomain(c)) { // for each i point in U such that c is in domain of var i
                lb_prime = ((nb_points_to_add[c] + sizeCluster[c] - 1)*lb_schedule[c][1] + s2[i][c] + s3[i][nb_points_to_add[c] - 1]) / (nb_points_to_add[c] + sizeCluster[c]);

                V_prime = lb_except + lb_prime; // Add updated contribution of cluster c

                // If new objective exceeds incumbent cost
                if (V_prime >= _V.getMax()) {
                    _X[setU_unassigned[i]].removeValue(c);
                }
            }
        }
    }
}


// Function which returns an engine handle (IlcConstraintI*) for the constraint implementation IlcWCSS_StandardCardControlI
IlcConstraint IlcWCSS_StandardCardControl(IlcIntVarArray X, IlcFloatVar V, const Data& data) {
    IlcCPEngine cp = X.getCPEngine();
    return new (cp.getHeap()) IlcWCSS_StandardCardControlI(cp, X, V, data);
}


// Macro which wraps the engine constraint handle into a modeling layer (Concert Technology) extractable object
ILOCPCONSTRAINTWRAPPER3(IloWCSS_StandardCardControl, cp, IloIntVarArray, _Xo, IloFloatVar, _Vo, const Data*, _datao) {
    use(cp, _Xo);
    use(cp, _Vo);
    return IlcWCSS_StandardCardControl(cp.getIntVarArray(_Xo), cp.getFloatVar(_Vo), *_datao);
}