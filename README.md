
# card-const-MSSC

This is an extensible Constraint Programming (CP) framework for **exact resolution** of the Cardinality-Constrained Minimum Sum-of-Squares Clustering (MSSC) problem.

This work will appear in the *CPAIOR 2020* conference proceeding as a [full paper](https://cpaior2020.dbai.tuwien.ac.at/papers/) (submission 70).

## Usefulness

MSSC is a clustering problem where the task is to partition data points in &#8477;<sup>*s*</sup> such that the resulting clusters (also commonly called *classes*) have the lowest combined variance (ie, observations are closest to their respective clusters' centroids). MSSC (clustering in general) is an essential tool in data mining.

*card-const-MSSC* enables you to find the **globally optimal solution** to the constrained MSSC problem. In particular, it includes a powerful new Global Constraint for **efficient handling of strict cardinality constraints**.

These constraints can be encountered in various fields such as image segmentation, category management in business, document clustering, and workgroup composition. Cardinality constraints are also useful to reinforce the clustering process against the presence of outliers and eliminates the possibility of poor solutions with classes that are either too large or too small.

Traditionally, MSSC (an NP-Hard problem) is solved using heuristics and metaheuristics.

## Related work

*card-const-MSSC* is steeped in the following work:
> Dao TBH., Duong KC., Vrain C. (2017) Constrained clustering by constraint programming. *Artificial Intelligence*, *244*, 70–94. doi:[10.1016/j.artint.2015.05.006](https://doi.org/10.1016/j.artint.2015.05.006)

In "Constrained Clustering by Constraint Programming," Dao et al. define [a versatile CP framework](https://cp4clustering.github.io/) for the exact resolution of constrained clustering problems. Specifically, Dao et al. also define, in that same framework, a Global Constraint for solving the constrained Euclidean MSSC:
> Dao TBH., Duong KC., Vrain C. (2015) Constrained Minimum Sum of Squares Clustering by Constraint Programming. In: Pesant G. (eds) *Principles and Practice of Constraint Programming - CP 2015*. Lecture Notes in Computer Science, vol 9255. Springer, Cham. doi:[10.1007/978-3-319-23219-5_39](https://doi.org/10.1007/978-3-319-23219-5_39)

[A competing framework](https://github.com/Behrouz-Babaki/CCCG), based on Column Generation, also exists for the constrained MSSC:
> Babaki B., Guns T., Nijssen S. (2014) Constrained Clustering Using Column Generation. In: Simonis H. (eds) *Integration of AI and OR Techniques in Constraint Programming - CPAIOR 2014*. Lecture Notes in Computer Science, vol 8451. Springer, Cham. doi:[10.1007/978-3-319-07046-9_31](https://doi.org/10.1007/978-3-319-07046-9_31)

*card-const-MSSC*'s operation is fully compatible with that of the CP framework by Dao et al. It also includes, for convenience and easy comparison, a reimplementation of their 2015 MSSC Global Constraint above. *card-const-MSSC* is based on *IBM ILOG CPLEX Optimization Studio* whereas Dao et al. chose [Gecode](https://github.com/Gecode/gecode) for their work.

## Technical details

### Prerequisites

This framework has been implemented using *CP Optimizer Extensions* from *IBM ILOG CPLEX Optimization Studio*. It has been tested with *IBM ILOG CPLEX Optimization Studio* versions 12.8, 12.9. At this time, *CP Optimizer* engine extensions API is **only available in C++**.

If you are a faculty member or a student, you can get *IBM ILOG CPLEX Optimization Studio* for free from IBM.

### Complexity

*card-const-MSSC* offers 3 different constraints with propagation complexities as follow:
-  `IloWCSS` : time complexity of *O*(*kq*<sup>2</sup> log *q* + *qn*) and space complexity of *O*(*n*<sup>2</sup>);
-  `IloWCSS_StandardCardControl` : time complexity of *O*(*q*<sup>2</sup> log *q* + *qn*) and space complexity of *O*(*n*<sup>2</sup>);
-  `IloWCSS_NetworkCardControl` : time and space complexities dominated by CPLEX Optimizer's Network Simplex.

*card-const-MSSC* uses [IntegerValuePrecedence](https://github.com/mnhaouas/IntegerValuePrecedence) for value symmetry breaking.

## Usage

Include `card-const-MSSC.h` to your C++ sources to use the *card-const-MSSC* framework in your *Concert Technology* model.

To help you get started, `main.cpp` is supplied as an example of how to use *card-const-MSSC* in your application. The source code is commented throughout to explain the operation of its various components.

Different parts of *card-const-MSSC* communicate with each other via the `Data` struct which gathers all relevant information about the current problem instance:
- `Data::fileID` is an `std::string` which uniquely identifies the problem instance under resolution;
- `Data::N`, `Data::S`, `Data::K` are `int` values for the number of observations to partition, the number of features each observation has (dimension) and the number of classes (clusters) to produce (resp.);
- `Data::coordinates` is a `double**` `N`-by-`S` dynamic array where each line is an observation and each column is a feature (eg `coordinates[2][9]` is the 10<sup>th</sup> attribute of the 3<sup>rd</sup> observation);
- `Data::dissimilarities` is a `double**` `N`-by-`N` dynamic array where each element represents the squared Euclidean distance between the observations (eg `dissimilarities[2][9]` is the squared Euclidean distance between the 10<sup>th</sup> and the 3<sup>rd</sup> observation). This array is symmetric and its diagonal is made of zeros;
- `Data::memberships` is an `int*` `N`-element dynamic array which stores an initial solution (`memberships[i] = c` means observation `i` belongs to class `c` in the initial solution). `memberships` has elements between 0 and `K-1`;
- `Data::targetCardinalities` is an `int*` `K`-element dynamic array which stores the desired final class cardinalities (`targetCardinalities[c] = m` means class `c` should contain `m` observations). `targetCardinalities` should sum up to `N` and should agree with `memberships`.

The 3 constraints included in *card-const-MSSC* have the following prototypes:
```
IloConstraint                      IloWCSS(IloEnv env, IloIntVarArray X, IloFloatVar V, const  Data* data, const char* name = 0);
IloConstraint  IloWCSS_StandardCardControl(IloEnv env, IloIntVarArray X, IloFloatVar V, const  Data* data, const char* name = 0);
IloConstraint   IloWCSS_NetworkCardControl(IloEnv env, IloIntVarArray X, IloFloatVar V, const  Data* data, const char* name = 0);
```
where:
- `env` is the optimizer's `IloEnv` environment;
- `X` is the modeling layer handle `IloIntVarArray` for the `N`-variable array of integer representative variables that link observations to their cluster;
- `V` is the modeling layer handle `IloFloatVar` for the real variable representing the total Within Cluster Sum-of-Squares of the solution (which must be constrained in the *Concert Technology* model to take the value of the WCSS);
- `name` is an optional custom name given to the posted Constraint in the model.

`IloWCSS` is a reimplementation of the work of Dao et al. (2015). `IloWCSS_StandardCardControl` is an adaptation of `IloWCSS` where it is made more efficient for the case of Cardinality-Constrained MSSC. `IloWCSS_NetworkCardControl` leverages the resolution of Minimum Cost Flow (MCF) problems through CPLEX Optimizer (using Concert Technology) to more efficiently solve the Cardinality-Constrained MSSC.

The search strategy is passed to the engine via `IloCP::startNewSearch` as a [goal](https://www.ibm.com/support/knowledgecenter/SSSA5P_12.10.0/ilog.odms.cpo.help/CP_Optimizer/Advanced_user_manual/topics/goals_understand_overview.html) with the following prototype:
```
IloGoal  IloMSSCSearchStrategy(IloEnv env, IloIntVarArray vars, const Data& data, const SearchParameters& searchParameters, const bool& solFound);
```
where:
- `env` is the optimizer's `IloEnv` environment;
- `vars` is the modeling layer handle `IloIntVarArray` for the branching variables. In this context, they're the representative integer variables of observations `X`;
- `data` is the `Data` struct which gathers all information about the current problem instance (see above);
- `searchParameters` is the `SearchParameters` struct which contains search heuristic preferences (see `IlcMSSCSearchStrategy.h` for information);
- `solFound` is a `bool` which takes the value `true` once a first solution has been found using the engine's `IloCP::next` method (it exists in the scope where CP Optimizer engine `IloCP` is instantiated).

## Acknowledgement

I am grateful to my brilliant supervisors, [Pesant G.](https://www.polymtl.ca/expertises/en/pesant-gilles) and [Aloise D.](https://www.gerad.ca/en/people/daniel-aloise), for their support throughout my graduate studies. Thank you to [Babaki B.](https://behrouz-babaki.github.io/) as well as to [Olivier P.](https://github.com/PhilippeOlivier) who have been available to answer my questions.

I am also grateful to the *Natural Sciences and Engineering Research Council* of Canada for their financial support in making this project possible.