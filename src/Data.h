/*
 * Author: Haouas, Mohammed Najib - sept 5th, 2019
 *                 > Last modified: may 30th, 2020
 *
 * Definition of the problem data structure which connects all parts of the resolution framework
 */


#ifndef __DATA_H
#define __DATA_H

#include <string>


struct Data {
	// Instance file ID (could be any uniquely identifying string, eg. file name)
	std::string fileID;

	// Number of observations, number of features, number of classes (resp.)
	int N, S, K;

	// Coordinates of observations in R^S (N-by-S)
	double** coordinates;
	// Squared euclidean distances between observations (N-by-N, symmetric)
	double** dissimilarities;

	// Target memberships (N-element array of integers between 0 and K-1 incl.)
	int* memberships;
	// Target class cardinalities (K-element array of non-zero integers, must agree with N and memberships above)
	int* targetCardinalities;
};

#endif // !__DATA_H