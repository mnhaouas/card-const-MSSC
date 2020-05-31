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


IlcWCSSI::IlcWCSSI(IloCPEngine cp, IlcIntVarArray X, IlcFloatVar V, const Data& data) :
IlcConstraintI(cp), _X(X), _V(V), _dissimilarities(data.dissimilarities), _n(X.getSize()), _k(data.K) {
    // sets of points and their sizes
    setP_assigned = new (cp.getHeap()) std::vector<IlcInt>[_k]; // setP_assigned[c] = i means point i is assigned to cluster c

    // size of each cluster c
    sizeCluster = IlcIntArray(cp, _k);

    // lower bound of WCSS of each cluster c if we add m points to it
    lb_schedule = new (cp.getHeap()) IlcFloat*[_k];
    for (int i = 0; i < _k; i++)
        lb_schedule[i] = new (cp.getHeap()) IlcFloat[_n + 1]; // lb_schedule[c][m] is lower bound on WCSS if m points assigned to c

    // Sum of dissimilarities of each cluster c
    S1 = IlcFloatArray(cp, _k); // S1[c] = sum of dissimilarities (squared) of cluster c 

    // Sum of dissimilarities between each unassigned point and each cluster
    s2 = new (cp.getHeap()) IlcFloat*[_n]; // Allocation on engine heap for efficient and unified memory control
    for (int i = 0; i < _n; i++)
        s2[i] = new (cp.getHeap()) IlcFloat[_k]; // s2[x][c] = sum of dissimilarities between x and all points in cluster c

    // Smallest 1/2 contribution of each unassigned point together with m other points
    s3 = new (cp.getHeap()) std::vector<IlcFloat>[_n]; // s3[x][m] = smallest contribution of the m other free points brought along in addition to x

    // Array for dynamic prog for global lower bound computation, assign q points to clusters
    lb_global = new (cp.getHeap()) IlcFloat*[_k];
    for (int i = 0; i < _k; i++)
        lb_global[i] = new (cp.getHeap()) IlcFloat[_n + 1]; // lb_global[c][m] = lower bound on WCSS of c clusters if we assign m points to them

    // lb_except[m] = lb if we add m points to clusters 0 through _k-1 except c
    lb_except = IlcFloatArray(cp, _n);

    // lb_prime[m] = lb on wcss of cluster c when we add m points + point i to it
    lb_prime = IlcFloatArray(cp, _n);

    // This espsilon is subtracted from the computed lower bound to prevent false backtracking while comparing with upper bound due to rounding errors
    _epsc = 5e-5;
}

IlcWCSSI::~IlcWCSSI() {
    // Any dynamically allocated elements are allocated in the engine heap which manages memory for us
}


void IlcWCSSI::post() {
    for (int i = 0; i < _n; i++)
        _X[i].whenDomain(this);
    
    _V.whenRange(this);
}


void IlcWCSSI::propagate() {
    // Reset set & number of points assigned and unassigned
    setU_unassigned.clear(); q = 0;
    for (int i = 0; i < _k; i++)
        setP_assigned[i].clear();
    p = 0;

    // Populating sets
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

    // Size of each cluster c
    for (int c = 0; c < _k; c++)
        sizeCluster[c] = setP_assigned[c].size(); // sizeCluster[c] = m means c is size m

    // Lower bound of WCSS of each cluster c if we add m points to it > INIT
    for (int c = 0; c < _k; c++)
        for (int m = 0; m <= q; m++)
            lb_schedule[c][m] = 0; // lb_schedule[c][m] is lower bound on WCSS if m points assigned to c

    // Sum of dissimilarities of each cluster c
    for (int c = 0; c < _k; c++) {
        S1[c] = 0;
        for (int i = 0; i < sizeCluster[c] - 1; i++)
            for (int j = i + 1; j < sizeCluster[c]; j++)
                S1[c] += _dissimilarities[setP_assigned[c][i]][setP_assigned[c][j]];
    }

    // Sum of dissimilarities between each unassigned point and each cluster
    for (int i = 0; i < q; i++) { // for each unassigned point
        for (int c = 0; c < _k; c++) { // for each cluster
            if (_X[setU_unassigned[i]].isInDomain(c)) { // if unassigned point i can be assigned to cluster c
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
    for (int i = 0; i < q; i++) {
        s3[i].clear();
        for (int j = 0; j < q; j++) {
            // At first, we put all points in that list
            s3[i].push_back(_dissimilarities[setU_unassigned[i]][setU_unassigned[j]] / 2);
        }

        std::sort(s3[i].begin(), s3[i].end()); // sort distances for each point, first element 0 because d(x,x) = 0

        for (int j = 1; j < q; j++)
            s3[i][j] += s3[i][j - 1]; // Compute minimum 1/2 contributions, first element is 0
    }

    // Computing lower bound for each cluster
    for (int c = 0; c < _k; c++) { // for each cluster
        for (int m = 0; m <= q; m++) { // if we add m points to c, m = 0 we add nothing, m = 1 we add only x, m = q we add x and all of the unassigned points remaining
            std::vector<IlcFloat> s; // s of every point x with current m possibility << will be destroyed going out of scope
            for (int i = 0; i < q; i++) { // contribution for each unassigned point
                if (m > 0) {
                    assert(s3[i][0] == 0);
                    s.push_back(s2[i][c] + s3[i][m - 1]); // s3[i][m-1] = 0 for m = 1
                }
                else {
                    s.push_back(0); // not adding anything, this is actually useless
                }
            }

            std::sort(s.begin(), s.end()); // sorting candidate points...

            IlcFloat S2 = 0;
            for (int i = 0; i < m; i++) // ... of which we select the m points that induce the lowest cost
                S2 += s[i];

            if (sizeCluster[c] + m > 0)
                lb_schedule[c][m] = (S1[c] + S2) / (sizeCluster[c] + m);
            else
                lb_schedule[c][m] = 0;
        }
    }

    // Dynamic prog for global lower bound, assign q points to clusters
    for (int i = 0; i <= q; i++)
        lb_global[0][i] = lb_schedule[0][i]; // same because for both we assign points to one cluster, the cluster number 0

    for (int c = 1; c < _k; c++) {
        for (int m = 0; m <= q; m++) {
            lb_global[c][m] = IlcInfinity;

            for (int i = 0; i <= m; i++)
                if (lb_global[c - 1][i] + lb_schedule[c][m - i] < lb_global[c][m])
                    lb_global[c][m] = lb_global[c - 1][i] + lb_schedule[c][m - i];
        }
    }

    // Lower bound for all clusters
    _V.setMin(lb_global[_k - 1][q] - _epsc);

    for (int c = 0; c < _k; c++) { // for each value c in domains of points, ie for each cluster
        for (int m = 0; m < q; m++) {
            lb_except[m] = 0; // lb_except[m] = lb if we add m points to clusters 0 through _k-1 except c

            for (int j = m; j <= q; j++)
                if (lb_global[_k - 1][j] - lb_schedule[c][j - m] > lb_except[m])
                    lb_except[m] = lb_global[_k - 1][j] - lb_schedule[c][j - m];
        }

        for (int i = 0; i < q; i++) {
            if (_X[setU_unassigned[i]].isInDomain(c)) { // for each i point in U such that c is in domain of Gi
                for (int m = 0; m < q; m++)
                    lb_prime[m] = ((sizeCluster[c] + m)*lb_schedule[c][m] + s2[i][c] + s3[i][m]) / (sizeCluster[c] + m + 1);

                IlcFloat V_prime = IlcInfinity;

                // Find optimal updated WCSS
                for (int m = 0; m < q; m++) {
                    if (V_prime > (lb_except[q - 1 - m] + lb_prime[m]))
                        V_prime = (lb_except[q - 1 - m] + lb_prime[m]);
                }

                if (V_prime >= _V.getMax()) {
                    _X[setU_unassigned[i]].removeValue(c);
                }
            }
        }
    }
}


// Function which returns an engine handle (IlcConstraintI*) for the constraint implementation IlcWCSSI
IlcConstraint IlcWCSS(IlcIntVarArray X, IlcFloatVar V, const Data& data) {
    IlcCPEngine cp = X.getCPEngine();
    return new (cp.getHeap()) IlcWCSSI(cp, X, V, data);
}


// Macro which wraps the engine constraint handle into a modeling layer (Concert Technology) extractable object
ILOCPCONSTRAINTWRAPPER3(IloWCSS, cp, IloIntVarArray, _Xo, IloFloatVar, _Vo, const Data*, _datao) {
    use(cp, _Xo);
    use(cp, _Vo);
    return IlcWCSS(cp.getIntVarArray(_Xo), cp.getFloatVar(_Vo), *_datao);
}