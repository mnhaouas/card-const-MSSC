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


#include "IlcIntPrecedeBinary.h"


IlcIntPrecedeBinaryI::IlcIntPrecedeBinaryI(IloCPEngine cp, IlcIntVarArray X, IlcInt s, IlcInt t) :
IlcConstraintI(cp), _cp(cp), _X(X), _n(_X.getSize()), _s(s), _t(t) {
    // Cannot use reversible objects as automatic objects
    alpha = new (_cp.getHeap()) IlcRevInt(_cp, 0);
    beta = new (_cp.getHeap()) IlcRevInt(_cp);
    gamma = new (_cp.getHeap()) IlcRevInt(_cp);

    // Initialize
    //     This could also be done in IlcIntPrecedeBinaryI::propagate since,
    //     in CP Optimizer, IlcConstraintI::propagate is called when the constraint
    //     is posted while solving towards initial fixed-point condition.
    while (alpha->getValue() < _n && !_X[alpha->getValue()].isInDomain(_s)) {
        _X[alpha->getValue()].removeValue(_t);
        alpha->setValue(_cp, alpha->getValue() + 1);
    }

    beta->setValue(_cp, alpha->getValue());
    gamma->setValue(_cp, alpha->getValue());

    if (alpha->getValue() < _n) {
        _X[alpha->getValue()].removeValue(_t);

        while (!(gamma->getValue() == _n || (_X[gamma->getValue()].isFixed() && _X[gamma->getValue()].getValue() == _t))) {
            gamma->setValue(_cp, gamma->getValue() + 1);
        }

        updateBeta();
    }
}


IlcIntPrecedeBinaryI::~IlcIntPrecedeBinaryI() {}


ILCCTDEMON1(IlcIntPrecedeBinaryI_mainDemon, IlcIntPrecedeBinaryI, mainDemon, IlcInt, inProcessIndex);
ILCCTDEMON1(IlcIntPrecedeBinaryI_gammaDemon, IlcIntPrecedeBinaryI, gammaDemon, IlcInt, inProcessIndex);
void IlcIntPrecedeBinaryI::post() {
    for (IlcInt i = 0; i < _n; i++) {
        _X[i].whenDomain(IlcIntPrecedeBinaryI_mainDemon(_cp, this, i));
        _X[i].whenValue(IlcIntPrecedeBinaryI_gammaDemon(_cp, this, i));
    }
}


void IlcIntPrecedeBinaryI::propagate() {
    // Propagate routines handed off to specialized constraint demons
}


// Run everytime a domain changes for a variable in X, through IlcIntPrecedeBinaryI_mainDemon
void IlcIntPrecedeBinaryI::mainDemon(IlcInt inProcessIndex) {
    if (beta->getValue() <= gamma->getValue()) {
        if (inProcessIndex == alpha->getValue() && !_X[inProcessIndex].isInDomain(_s)) {
            alpha->setValue(_cp, alpha->getValue() + 1);

            while (alpha->getValue() < beta->getValue()) {
                _X[alpha->getValue()].removeValue(_t);
                alpha->setValue(_cp, alpha->getValue() + 1);
            }

            while (alpha->getValue() < _n && !_X[alpha->getValue()].isInDomain(_s)) {
                _X[alpha->getValue()].removeValue(_t);
                alpha->setValue(_cp, alpha->getValue() + 1);
            }

            if (alpha->getValue() < _n)
                _X[alpha->getValue()].removeValue(_t);

            beta->setValue(_cp, alpha->getValue());

            if (alpha->getValue() < _n)
                updateBeta();
        }
        else if (inProcessIndex == beta->getValue() && !_X[inProcessIndex].isInDomain(_s)) {
            updateBeta();
        }
    }
}


// Run everytime a variable in X is bound, through IlcIntPrecedeBinaryI_gammaDemon
void IlcIntPrecedeBinaryI::gammaDemon(IlcInt inProcessIndex) {
    assert(_X[inProcessIndex].isFixed());

    if (beta->getValue() < gamma->getValue() && inProcessIndex < gamma->getValue() && _X[inProcessIndex].getValue() == _t) {
        gamma->setValue(_cp, inProcessIndex);

        if (beta->getValue() > inProcessIndex)
            _X[alpha->getValue()].setValue(_s);
    }
}


// Refer to cited paper above for information about this method
void IlcIntPrecedeBinaryI::updateBeta() {
    while ( !(beta->getValue() == _n || _X[beta->getValue()].isInDomain(_s)) )
        beta->setValue(_cp, beta->getValue() + 1);

    if (beta->getValue() > gamma->getValue())
        _X[alpha->getValue()].setValue(_s);
}


// Function which returns an engine handle (IlcConstraintI*) for the constraint implementation IlcIntPrecedeBinaryI
IlcConstraint IlcIntPrecedeBinary(IlcIntVarArray X, IlcInt s, IlcInt t) {
    IlcCPEngine cp = X.getCPEngine(); // Get CP engine from variable array
    return new (cp.getHeap()) IlcIntPrecedeBinaryI(cp, X, s, t); // Allocate implementation object on the engine heap for efficient memory management
}


// Macro which wraps the engine constraint handle into a modeling layer (Concert Technology) extractable object
ILOCPCONSTRAINTWRAPPER3(IloIntPrecedeBinary, cp, IloIntVarArray, _Xo, IloInt, _so, IloInt, _to) {
    use(cp, _Xo); // Force extraction of modeling layer extractables (ie, get engine level objects)
    return IlcIntPrecedeBinary(cp.getIntVarArray(_Xo), _so, _to);
}