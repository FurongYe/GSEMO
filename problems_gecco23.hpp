#define FMT_HEADER_ONLY
#include "ioh.hpp"

double onemax_constraint(const std::vector<int> &x);

double leadingones_constraint(const std::vector<int> &x);

double onejump(const std::vector<int> &x);

double jump_constraint(const std::vector<int> &x);

std::shared_ptr<ioh::problem::IntegerSingleObjective> get_problem(const int id, const int dimension);
