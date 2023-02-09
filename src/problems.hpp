#ifndef GSEMO_PROBLEM
#define GSEMO_PROBLEM 

#define FMT_HEADER_ONLY
#include "ioh.hpp"
extern std::vector<std::array<double, 2> > pre_pareto_y;
extern std::array<double, 2> pre_hv_reference;

double onemax_constraint(const std::vector<int> &x);

double COCZ_constraint(const std::vector<int> &x);

double leadingones_constraint(const std::vector<int> &x);

double onejump(const std::vector<int> &x);

double jump_constraint(const std::vector<int> &x);

std::shared_ptr<ioh::problem::IntegerSingleObjective> get_problem(const int id, const int dimension);

#endif