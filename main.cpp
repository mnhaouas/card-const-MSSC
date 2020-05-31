/*
 * Author: Haouas, Mohammed Najib - may 30th, 2020
 *                 > Last modified: may 30th, 2020
 *
 * Example Concert Technology model using card-const-MSSC framework to help you get started with card-const-MSSC.
 * 
 * Note: this file is an implementation example intended to give a general idea
 *       about how card-const-MSSC can be used. You are expected to modify it
 *       to suit your needs. Because data parsing is left to the user, this file will
 *       NOT compile as supplied here.
 * 
 * card-const-MSSC is my (Haouas, M.N.) MSc's research project under supervision of Pesant, G. & Aloise, D.
 *
 * Work pending publication in the CPAIOR 2020 conference proceeding as a full paper (submission 70).
 *     https://cpaior2020.dbai.tuwien.ac.at/papers/
 */


#include <ilcp/cp.h> // Using (IBM ILOG CPLEX Optimization Studio) CP Optimizer
#include "card-const-MSSC.h" // Using card-const-MSSC


int main(int argc, char** argv) {
    /* 
     * Start by building the problem data, eg. by reading a file. Refer to Data.h for information.
     */

    Data data;
    /* ... */


    /* 
     * Set up search preferences. Refer to IloMSSCSearchStrategy.h for information.
     */

    SearchParameters searchParameters;

    // For example...
    searchParameters.initialSolution = CustomCPSearchOptions::InitialSolution::MEMBERSHIPS_AS_INDICATED;
    searchParameters.mainSearch = CustomCPSearchOptions::MainSearch::MAX_MIN_VAR;
    searchParameters.tieHandling = CustomCPSearchOptions::TieHandling::UNBOUND_FARTHEST_TOTAL_SS;


    /*
     * CP model and optimizer.
     *     Resolution using CP Optimizer through Concert Technology.
     */

    IloEnv env; // Set up environment
    try {
        // Create model in env
        IloModel model(env);

        // VARIABLES: Representative and auxiliary variables
        IloIntVarArray x(env, data.N, 0, (data.K - 1)); // Observation representative variables, size N array, domains 0..K-1
        IloFloatVar V(env, 0, IloInfinity); // Objective variable, total within cluster sum of squares (WCSS). domain 0..inf
        IloIntVarArray cardinality(env, data.K, 1, data.N); // Clusters' cardinalities, size K array, domains 1..N
        model.add(V);
        model.add(x);

        // CONSTRAINT: Link cardinality to actual cardinalities through Global Cardinality Constraint (GCC)
        IloIntArray vals(env, data.K);
        for (int i = 0; i < data.K; i++)
            vals[i] = i;
        model.add(IloDistribute(env, cardinality, vals, x));

        // BRAIN: MSSC resolution constraints
            // // CONSTRAINT: Total WCSS lower bound
            // model.add(IloWCSS(env, x, V, &data));

            // // *or* CONSTRAINT: Total WCSS lower bound with external cardinality control
            // model.add(IloWCSS(env, x, V, &data));
            // for (int c = 0; c < data.K; c++)
            //     model.add(cardinality[c] == data.targetCardinalities[c]);

            // // *or* CONSTRAINT: Total WCSS lower bound with standard internal cardinality control
            // model.add(IloWCSS_StandardCardControl(env, x, V, &data));
            
            // *or* CONSTRAINT: Total WCSS lower bound with MCF-based internal cardinality control
            model.add(IloWCSS_NetworkCardControl(env, x, V, &data));

        // CONSTRAINT: Binding objective variable to WCSS using actual expression
        IloExprArray wcsd(env, data.K); // wcsd[c] = within cluster sum of dissimilarities for cluster c
        for (int i = 0; i < data.K; i++)
            wcsd[i] = IloNumExpr(env); // Empty handles are now null expressions

        for (int i = 0; i < data.K; i++) // for each cluster
            for (int pt1 = 0; pt1 < (data.N - 1); pt1++)
                for (int pt2 = (pt1 + 1); pt2 < data.N; pt2++) // for each pair of points
                    wcsd[i] += ((x[pt1] == i)*(x[pt2] == i)*data.dissimilarities[pt1][pt2]);

        IloNumExpr ub_V_exp(env); // V for current solution
        for (int i = 0; i < data.K; i++)
            ub_V_exp += (wcsd[i] / cardinality[i]);

        model.add(V == ub_V_exp);

        // SYM BREAKING: Pair-wise int value precedence for breaking value symmetry
        for (int i = 1; i < data.K; i++)
            model.add(IloIntPrecedeBinary(env, x, i-1, i));

        // OBJECTIVE: Minimize total WCSS
        model.add(IloMinimize(env, V));

        // SEARCH STRATEGY: Custom search heuristic
        bool solFound = false; // Witness for initial solution found
        IloGoal masterSearch = IloMSSCSearchStrategy(env, x, data, searchParameters, solFound); // Initial goal

        // ENGINE: Creating and configuring CP algorithm
        IloCP cp(model);
        // cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet); // Uncomment to suppress CP Optimizer output
        // cp.setParameter(IloCP::TimeLimit, INT_TIME_IN_SECONDS); // Uncomment to set time limit of INT_TIME_IN_SECONDS

        // RESOLUTION: Initialize solve process
        /* NOTE: per https://www.ibm.com/developerworks/community/forums/html/topic?id=02b1d19b-cc6b-4200-b474-277fdcb0b876,
         *       using the IloCP::solve autosearch is a bad idea to control bounds because there may be relaxations that
         *       lead to poorer bounds than the ones given by the incumbent solution.
         */
        cp.startNewSearch(masterSearch); // Start new search using custom search (by supplying initial goal)
                                         // NOTE: CP Optimizer moves towards initial fixed-point condition here
                                         //       All propagate member functions present are run.

        // RESOLUTION: Subsequent search
        while (cp.next()) {
            solFound = true; // At least one solution is found, so set to true

            // Display intermediate solutions, for example...
            cp.out() << std::endl << std::endl << "Status: " << cp.getStatus() << std::endl;
            cp.out() << "  V = " << cp.getObjValue() << std::endl;

            cp.out() << "  Corresponding memberships: " << std::endl << "  ";
            for (int i = 0; i < data.N; i++) {
                cp.out() << cp.getValue(x[i]) << " ";
                if ((i + 1) % 24 == 0)
                    cp.out() << "..." << std::endl << "  ";
            }
            cp.out() << std::endl;

            cp.out() << "  Cluster cardinalities: " << std::endl << "  ";
            for (int c = 0; c < data.K; c++)
                cp.out() << cp.getValue(cardinality[c]) << " ";
            cp.out() << std::endl;

            cp.out() << "  Cumulative solve duration: " << cp.getTime() << std::endl;
        }


        /*
         * Final print.
         */

        // Display final stats, for example...
        cp.out() << std::endl << std::endl << ">> Done. Status: " << cp.getStatus() << std::endl;
        cp.out() << "Number of branches  : " << cp.getInfo(IloCP::IntInfo::NumberOfBranches) << std::endl;
        cp.out() << "Number of fails     : " << cp.getInfo(IloCP::IntInfo::NumberOfFails) << std::endl;
        cp.out() << "Total solve duration: " << cp.getTime() << std::endl;
    }
    catch (IloException& ex) {
        env.out() << "Error: " << ex << std::endl;
    }


    /*
     * Disposition
     */

    env.end();
}