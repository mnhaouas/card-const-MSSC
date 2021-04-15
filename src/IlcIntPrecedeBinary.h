/*
 * This constraint ensures integer value precedence of value s over value t across
 *     integer variable array X and proceeds to the appropriate filtering.
 * This constraint is useful for breaking value symmetries.
 * This constraint maintains Generalized Arc Consistency (GAC).
 * A sequence of precedence can be maintained by posting this constraint for each pair or values.
 * 
 * Main arguments: * X, integer variables over which precedence must be maintained.
 *
 * Additional arguments: * s, antecedent value.
 *                       * t, subsequent value.
 *
 * Note: posting on all pairs of values ensures strictly stronger
 *     filtering but is inefficient in practice and does not produce
 *     better results. Therefore, it is recommended to post on pairs 
 *     of adjacent values only.
 *
 * Implementation of the work of:
 * Law Y.C., Lee J.H.M. (2004) Global Constraints for Integer and Set Value Precedence.
 *     In: Wallace M. (eds) Principles and Practice of Constraint Programming – CP 2004. CP 2004.
 *     Lecture Notes in Computer Science, vol 3258. Springer, Berlin, Heidelberg
 *     doi:10.1007/978-3-540-30201-8_28
 */


// Using (IBM ILOG CPLEX Optimization Studio) CP Optimizer Extensions
#include <ilcp/cpext.h>


// Variable names are kept consistent with Law et al.'s nomenclature.
//     Refer to cited paper above for any additional explanation of the filtering algorithm.


ILOSTLBEGIN


class IlcIntPrecedeBinaryI : public IlcConstraintI {
protected:
    void updateBeta();

    IlcIntVarArray _X;
    IlcInt _s, _t;

    IlcRevInt* alpha;
    IlcRevInt* beta;
    IlcRevInt* gamma;

    IlcInt _n; // Size of integer variable array over which this constraint is posted

    IlcCPEngine _cp;


public:
    IlcIntPrecedeBinaryI(IloCPEngine cp, IlcIntVarArray X, IlcInt s, IlcInt t);
    ~IlcIntPrecedeBinaryI();
    virtual void propagate();
    virtual void post();

    // Propagation depends on trigger. Therefore, specilized demons are called for filtering.
    //     inProcessIndex is the index of the variable in array X which triggered a propagation event.
    void mainDemon(IlcInt inProcessIndex);
    void gammaDemon(IlcInt inProcessIndex);
};


// Forward declaration of the function which returns the engine constraint handle
IlcConstraint IlcIntPrecedeBinary(IlcIntVarArray X, IlcInt s, IlcInt t);