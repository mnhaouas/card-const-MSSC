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


// Computes the delta objective when pt is assigned to cluster c
int getDeltaObjective(IlcIntVarArray vars, IlcInt pt, IlcInt c, double const* const* const dissimilarities) {
    double S1 = 0, S2 = 0;
    int card_cluster = 0;

    // Compute S1 WCSD for cluster c, S2 contribution of pt to c
    for (int i = 0; i < vars.getSize(); i++) {
        if (vars[i].isFixed() && (vars[i].getValue() == c)) {
            card_cluster++; // We found a point in c

            for (int j = i + 1; j < vars.getSize(); j++)
                if (vars[j].isFixed() && (vars[j].getValue() == c))
                    S1 += dissimilarities[i][j]; // we found a point in the same cluster so we add the distance between'em

            S2 += dissimilarities[i][pt]; // we add the distance between the point in c and candidate pt
        }
    }

    if (card_cluster == 0)
        return 0; // Assigning point to empty cluster has 0 cost

    // Multiply by 1000 for precision.
    return (int) (((S1 + S2) / (card_cluster + 1) - S1 / card_cluster) * 1000);
}


// Computes total SS between pt and all points in U.
int getUnboundPointsTotalSS(IlcIntVarArray vars, IlcInt pt, double const* const* const dissimilarities) {
    // Init total distance between pt and all points
    double total_dist = 0;

    for (int i = 0; i < vars.getSize(); i++)
        if (!vars[i].isFixed())
            total_dist += dissimilarities[i][pt];

    return (int) (total_dist * 100);
}


// Gives int score for (squared) distance between i and j
inline int getIntDist(IlcInt i, IlcInt j, double const* const* const dissimilarities) {
    return (int) (dissimilarities[i][j] * 1000);
}


// Strategy as a goal to be given to CP Optimizer engine
ILCGOAL4(IlcMSSCSearchStrategy, IlcIntVarArray, vars, const Data&, data, const SearchParameters&, searchParameters, const bool&, solFound) {
    /*
     * Initializations
     */

    IlcBool found = IlcFalse; // Check if a choice can be made. Otherwise, end of search for current solution reached

    IlcInt bestI; // Chosen variable
    IlcInt bestJ; // Chosen value for variable

    int min_contrib_loco; // min contribution of a variable
    int max_contrib_globo = 0; // max min_contrib_loco observed


    /*
     * Initial solution handling: generation of first solution against which further filtering will happen
     *     NOTE: part of this handling happens externally and is dependent on the order of the observations.
     */
    
    if (!solFound && searchParameters.initialSolution != CustomCPSearchOptions::InitialSolution::NONE) {
        switch (searchParameters.initialSolution) {
            case CustomCPSearchOptions::InitialSolution::GREEDY_INIT: {
                // find min size domain
                int minimum_domain_size = __MAX_INT;
                for (IlcInt i = 0; i < vars.getSize(); i++)
                    if (!vars[i].isFixed() && vars[i].getSize() < minimum_domain_size)
                        minimum_domain_size = vars[i].getSize();

                min_contrib_loco = __MAX_INT;

                for (IlcInt i = 0; i < vars.getSize(); i++) {
                    if (!vars[i].isFixed() && vars[i].getSize() == minimum_domain_size) {
                        found = IlcTrue;

                        for (IlcIntExpIterator iter(vars[i]); iter.ok(); ++iter) {
                            IlcInt j = *iter;

                            int aggr = getDeltaObjective(vars, i, j, data.dissimilarities);

                            if (aggr < min_contrib_loco) {
                                min_contrib_loco = aggr;
                                bestI = i;
                                bestJ = j;
                            }
                        }
                    }
                }
            }
            break;

            case CustomCPSearchOptions::InitialSolution::MEMBERSHIPS_AS_INDICATED: {
                // Choose unbound variable...
                for (IlcInt i = 0; i < vars.getSize(); i++) {
                    if (!vars[i].isFixed()) {
                        found = IlcTrue;
                        bestI = i;
                        break;
                    }
                }

                // ...bind it as instructed
                bestJ = data.memberships[bestI];
            }
            break;
        }
        
        // Return choice
        if (!found)
            return 0;

        // Check choice is consistent. Needed in case file loaded has bad memberships.
        assert(((bestI < vars.getSize()) && (bestI >= 0) && (bestJ < data.K) && (bestJ >= 0) && found) || (!found));

        // Binary branching
        return IlcOr(IlcAnd(vars[bestI] == bestJ, this),
                     IlcAnd(vars[bestI] != bestJ, this)
                     );
    }
    

    /*
     * Subsequent main search
     */

    switch (searchParameters.mainSearch) {
        case CustomCPSearchOptions::MainSearch::MAX_MIN_VAR: {
            IlcInt bestinterimJ;
            for (IlcInt i = 0; i < vars.getSize(); i++) {
                if (!vars[i].isFixed()) {
                    found = IlcTrue;
                    min_contrib_loco = __MAX_INT;

                    // Explore vars[i] domain
                    for (IlcIntExpIterator iter(vars[i]); iter.ok(); ++iter) {
                        IlcInt j = *iter;

                        int aggr = getDeltaObjective(vars, i, j, data.dissimilarities);

                        if (aggr < min_contrib_loco) {
                            min_contrib_loco = aggr;
                            bestinterimJ = j;
                        }
                    }

                    if (min_contrib_loco >= max_contrib_globo) {
                        max_contrib_globo = min_contrib_loco;
                        bestI = i;
                        bestJ = bestinterimJ;
                    }
                }
            }
        }
        break;
    }


    /*
     * Tie breaking: happens when, while backtracking, a cluster becomes empty.
     */

    if (max_contrib_globo == 0 && found) {
        /*
         * Before we continue, let's start by figuring out exactly which cluster(s) is/are empty
         */

        // Fun fact: it's not necessarily the highest assigned number + 1 like I thought because of sym breaking.
        //     It may not violated if clusters are skipped.
        int sk_cluster_jump = -1; // Cluster number where a jump occurs
        bool jump_happened = false;
        std::vector<int> occupied_clusters; occupied_clusters.reserve(data.K + 1); // List of non-empty clusters

        occupied_clusters.push_back(-1); // The back is always the highest occupied cluster index
        int sk_cluster_to_fill = -1; // Cluster selected to be filled for the tie-breaking

        for (int i = 0; i < vars.getSize(); i++) {
            if (vars[i].isFixed()) {
                if ((vars[i].getValue() - occupied_clusters.back()) >= 2) { // Skipped cluster, we're done, we know what the next empty cluster is
                    jump_happened = true;
                    sk_cluster_jump = occupied_clusters.back();
                    occupied_clusters.push_back(vars[i].getValue());
                    // Yet we don't break because occupied_clusters could be useful...
                }
                else if ((vars[i].getValue() - occupied_clusters.back()) == 1) {
                    occupied_clusters.push_back(vars[i].getValue());
                }
            }
        }

        if (jump_happened) {
            sk_cluster_to_fill = sk_cluster_jump + 1;
        }
        else if (occupied_clusters.back() <= (data.K - 2)) {
            sk_cluster_to_fill = occupied_clusters.back() + 1;
        }
        else { // There is no tie to break, piss off (odds of this happening are ultra slim)
            // Return choice
            if (!found)
                return 0;

            // Check choice is consistent. Needed in case file loaded has bad memberships.
            assert(((bestI < vars.getSize()) && (bestI >= 0) && (bestJ < data.K) && (bestJ >= 0) && found) || (!found));
            return IlcOr(IlcAnd(vars[bestI] == bestJ, this),
                         IlcAnd(vars[bestI] != bestJ, this)
                         ); // Binary branching
        }

        occupied_clusters.erase(occupied_clusters.begin());

        // Couple assertions to be sure...
        assert(occupied_clusters.size() < data.K);
        assert(sk_cluster_to_fill < data.K && sk_cluster_to_fill >= 0);


        /*
         * Now we can do our work...
         */

        switch (searchParameters.tieHandling) {
            // ***** Start empty cluster farthest to un-assigned points on the basis of total sum of squares ***** 
            case CustomCPSearchOptions::TieHandling::UNBOUND_FARTHEST_TOTAL_SS: {
                int comp_dist;
                int max_dist = 0;

                // Look for farthest point
                for (IlcInt i = 0; i < vars.getSize(); i++) {
                    if (!vars[i].isFixed() && vars[i].isInDomain(sk_cluster_to_fill)) {
                        comp_dist = getUnboundPointsTotalSS(vars, i, data.dissimilarities);

                        // ...break tie with farthest point.
                        if (comp_dist > max_dist) {
                            max_dist = comp_dist;
                            bestI = i;
                            bestJ = sk_cluster_to_fill;
                        }
                    }
                }
            }
            break;

            // ***** Start empty cluster farthest to fixed points on the basis of distance alone *****
            case CustomCPSearchOptions::TieHandling::FIXED_FARTHEST_DIST: {
                int max_dist = 0;

                // Look for farthest point
                for (int i = 0; i < vars.getSize(); i++) {
                    if (!vars[i].isFixed() && vars[i].isInDomain(sk_cluster_to_fill)) {
                        for (int j = 0; j < vars.getSize(); j++) {
                            if (vars[j].isFixed() && getIntDist(i, j, data.dissimilarities) > max_dist) {
                                max_dist = getIntDist(i, j, data.dissimilarities);
                                bestI = i;
                                bestJ = sk_cluster_to_fill;
                            }
                        }
                    }
                }
            }
            break;

            // ***** Start empty cluster at the point that has the maximum distance to its closest cluster *****
            case CustomCPSearchOptions::TieHandling::FIXED_MAX_MIN: {
                int max_dist_overall = 0;
                int max_dist_overall_indexof = bestI;
                for (int i = 0; i < vars.getSize(); i++) {
                    if (!vars[i].isFixed() && vars[i].isInDomain(sk_cluster_to_fill)) { // getting [PT]Min
                        int min_dist_to_all_clusters = __MAX_INT;
                        int min_dist_to_all_clusters_indexof;

                        for (std::vector<int>::iterator it = occupied_clusters.begin(); it != occupied_clusters.end(); it++) { // gettint [PT]nMin
                            int current_cluster = *it;
                            int min_dist_to_current_cluster = __MAX_INT;
                            int min_dist_to_current_cluster_indexof;

                            for (int j = 0; j < vars.getSize(); j++) {
                                if (vars[j].isFixed() && vars[j].getValue() == current_cluster && getIntDist(i, j, data.dissimilarities) < min_dist_to_current_cluster) {
                                    min_dist_to_current_cluster = getIntDist(i, j, data.dissimilarities);
                                    min_dist_to_current_cluster_indexof = i;
                                }
                            }

                            if (min_dist_to_current_cluster < min_dist_to_all_clusters) {
                                min_dist_to_all_clusters = min_dist_to_current_cluster;
                                min_dist_to_all_clusters_indexof = min_dist_to_current_cluster_indexof;
                            }
                        }

                        if (min_dist_to_all_clusters > max_dist_overall) {
                            max_dist_overall = min_dist_to_all_clusters;
                            max_dist_overall_indexof = min_dist_to_all_clusters_indexof;
                        }
                    }
                }

                bestI = max_dist_overall_indexof;
                bestJ = sk_cluster_to_fill;
            }
            break;

            // ***** Start empty cluster at the point that is farthest to the biggest cluster center *****
            case CustomCPSearchOptions::TieHandling::FARTHEST_POINT_FROM_BIGGEST_CENTER: {
                int* cards = new int[data.K];
                for (int i = 0; i < data.K; i++)
                    cards[i] = 0;

                // Determine cardinalities
                for (int i = 0; i < vars.getSize(); i++)
                    if (vars[i].isFixed())
                        cards[vars[i].getValue()] += 1;

                // Determine biggest cluster
                int biggest_cluster_index; int biggest_card = 0;
                for (int i = 0; i < data.K; i++) {
                    if (cards[i] > biggest_card) {
                        biggest_card = cards[i];
                        biggest_cluster_index = i;
                    }
                }

                if (biggest_card == 0) {
                    if (!found)
                        return 0;

                    assert(((bestI < vars.getSize()) && (bestI >= 0) && (bestJ < data.K) && (bestJ >= 0) && found) || (!found));
                    return IlcOr(IlcAnd(vars[bestI] == bestJ, this),
                                 IlcAnd(vars[bestI] != bestJ, this)
                                 ); // Binary branching
                }

                // Determine center of biggest cluster
                double* biggest_cluster_center = new double[data.S];
                for (int i = 0; i < data.S; i++)
                    biggest_cluster_center[i] = 0;

                for (int i = 0; i < vars.getSize(); i++) {
                    if (vars[i].isFixed() && vars[i].getValue() == biggest_cluster_index) {
                        for (int j = 0; j < data.S; j++) {
                            biggest_cluster_center[j] += data.coordinates[i][j];
                        }
                    }
                }

                for (int i = 0; i < data.S; i++)
                    biggest_cluster_center[i] /= biggest_card;

                // Find farthest point
                double biggest_distance = 0;
                for (int i = 0; i < vars.getSize(); i++) {
                    double temp_dist_pts = 0;
                    for (int j = 0; j < data.S; j++) {
                        temp_dist_pts += (biggest_cluster_center[j] - data.coordinates[i][j])*(biggest_cluster_center[j] - data.coordinates[i][j]);
                    }

                    if (temp_dist_pts > biggest_distance && vars[i].isInDomain(sk_cluster_to_fill)) {
                        biggest_distance = temp_dist_pts;
                        bestI = i;
                    }
                }
                bestJ = sk_cluster_to_fill;

                delete[] cards; delete[] biggest_cluster_center;
                cards = NULL; biggest_cluster_center = NULL;
            }
            break;

            // ***** Start empty cluster at the point that has maximum minimum distance to all cluster centers *****
            case CustomCPSearchOptions::TieHandling::MAX_MIN_POINT_FROM_ALL_CENTER: {
                int* cards = new int[data.K];
                for (int i = 0; i < data.K; i++)
                    cards[i] = 0;

                // Determine max occupied cluster and cardinalities
                for (int i = 0; i < vars.getSize(); i++)
                    if (vars[i].isFixed())
                        cards[vars[i].getValue()] += 1;

                // Determine center of clusters
                double** cluster_center = new double*[data.K];
                for (int i = 0; i < data.K; i++) {
                    cluster_center[i] = new double[data.S];
                    for (int j = 0; j < data.S; j++) {
                        cluster_center[i][j] = 0;
                    }
                }

                for (std::vector<int>::iterator iter = occupied_clusters.begin(); iter != occupied_clusters.end(); iter++)
                    assert(cards[*iter] > 0);

                for (std::vector<int>::iterator iter = occupied_clusters.begin(); iter != occupied_clusters.end(); iter++) {
                    int c = *iter;
                    for (int i = 0; i < vars.getSize(); i++) {
                        if (vars[i].isFixed() && vars[i].getValue() == c) {
                            for (int j = 0; j < data.S; j++) {
                                cluster_center[c][j] += data.coordinates[i][j];
                            }
                        }
                    }
                }

                for (std::vector<int>::iterator iter = occupied_clusters.begin(); iter != occupied_clusters.end(); iter++) {
                    int c = *iter;
                    for (int i = 0; i < data.S; i++)
                        cluster_center[c][i] /= cards[c];
                }

                // Find farthest point to centers
                double max_distance_global = 0;
                for (int i = 0; i < vars.getSize(); i++) {
                    if (!vars[i].isFixed() && vars[i].isInDomain(sk_cluster_to_fill)) {
                        double smallest_distance_local = std::numeric_limits<double>::infinity();
                        // for (int c = 0; c <= max_occupied_cluster; c++) {
                        for (std::vector<int>::iterator iter = occupied_clusters.begin(); iter != occupied_clusters.end(); iter++) {
                            int c = *iter;
                            double temp_dist_pts = 0;
                            for (int j = 0; j < data.S; j++) {
                                temp_dist_pts += (cluster_center[c][j] - data.coordinates[i][j])*(cluster_center[c][j] - data.coordinates[i][j]);
                            }

                            if (temp_dist_pts < smallest_distance_local) {
                                smallest_distance_local = temp_dist_pts;
                            }
                        }

                        if (smallest_distance_local > max_distance_global) {
                            max_distance_global = smallest_distance_local;
                            bestI = i;
                        }
                    }
                }
                bestJ = sk_cluster_to_fill;

                delete[] cards;
                for (int i = 0; i < data.K; i++)
                    delete[] cluster_center[i];

                delete[] cluster_center;
                cards = NULL; cluster_center = NULL;
            }
            break;

            case CustomCPSearchOptions::TieHandling::NONE: {
                // Nothing to do
            }
            break;
        }
    }


    /*
     * Return final choice
     */

    if (!found)
        return 0;

    assert(((bestI < vars.getSize()) && (bestI >= 0) && (bestJ < data.K) && (bestJ >= 0) && found) || (!found));
    return IlcOr(IlcAnd(vars[bestI] == bestJ, this),
                 IlcAnd(vars[bestI] != bestJ, this)
                 ); // Binary branching
}


// Macro which wraps the engine goal into a modeling layer (Concert Technology) object
ILOCPGOALWRAPPER4(IloMSSCSearchStrategy, cp, IloIntVarArray, varso, const Data&, datao, const SearchParameters&, searchParameterso, const bool&, solFoundo) {
    return IlcMSSCSearchStrategy(cp, cp.getIntVarArray(varso), datao, searchParameterso, solFoundo);
}