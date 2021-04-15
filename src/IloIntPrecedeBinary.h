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


// Include this header file into your CP model source file to make this constraint available in Concert Technology.


#include "IlcIntPrecedeBinary.h"


// Forward declaration of the function which returns a modeling layer (Concert Technology) constraint object, which wraps an engine constraint handle.
IloConstraint IloIntPrecedeBinary(IloEnv env, IloIntVarArray X, IloInt s, IloInt t, const char* name = 0);