/*
 * Author: Haouas, Mohammed Najib - may 30th, 2020
 *                 > Last modified: may 30th, 2020
 *
 * Include this header to your C++ sources to use the card-const-MSSC framework in your Concert Technology model.
 * Refer to git repository Readme file for usage instructions.
 *
 * This is part of my (Haouas, M.N.) MSc's research project under supervision of Pesant, G. & Aloise, D.
 *
 * Work pending publication in the CPAIOR 2020 conference proceeding as a full paper (submission 70).
 *     https://cpaior2020.dbai.tuwien.ac.at/papers/
 */


#ifndef __CARD_CONST_MSSC_H
#define __CARD_CONST_MSSC_H

// Problem data structure
#include "src/Data.h"

// Search strategy
#include "src/IloMSSCSearchStrategy.h"

// Constraints
#include "src/IloIntPrecedeBinary.h" // Symmetry breaking constraint, based on Integer Value Precedence.
#include "src/IloWCSS.h" // Constraint speeds up resolution of general MSSC through CP
#include "src/IloWCSS_StandardCardControl.h" // Constraint speeds up resolution of cardinality-constrained MSSC through CP, based on IloWCSS
#include "src/IloWCSS_NetworkCardControl.h" // Constraint speeds up resolution of cardinality-constrained MSSC through CP, based on MCF resolution

#endif // !__CARD_CONST_MSSC_H