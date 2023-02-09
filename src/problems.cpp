#include "problems.hpp"


std::vector<std::array<double, 2> > pre_pareto_y;
std::array<double, 2> pre_hv_reference = {-1,-1};

double onemax_constraint(const std::vector<int> &x)
{
    return (double)x.size() - std::accumulate(x.begin(), x.end(), 0.0);
}

double leadingones_constraint(const std::vector<int> &x)
{
    int result = 0;
    for (int i = x.size() - 1; i >= 0; --i) {
        if (x[i] == 0) {
            result++;
        }
        else {
            break;
        }
    }
    return static_cast<double>(result);
}


double COCZ(const std::vector<int> &x)
{
    return std::accumulate(x.begin(), x.end(), 0.0);
}

double COCZ_constraint(const std::vector<int> &x) {
    int result = 0;
    for (int i = 0; i != x.size(); ++i) {
        if (i < x.size()/2) {
            result += (x[i]);
        }
        else {
            result += (1-x[i]);
        }
    }
    return (double)(result);
}

double onejump(const std::vector<int> &x)
{
    int k = 2;
    int xones = std::accumulate(x.begin(), x.end(), 0.0);
    double result;

    if (xones == x.size() || xones <= x.size() - k)
    {
        result = xones + k;
    }
    else
    {
        result = x.size() - xones;
    }
    return static_cast<double>(result);
}

double jump_constraint(const std::vector<int> &x)
{
    int k = 2;
    int xzeros = x.size() - std::accumulate(x.begin(), x.end(), 0.0);
    double result;

    if (xzeros == x.size() || xzeros <= x.size() - k)
    {
        result = xzeros + k;
    }
    else
    {
        result = x.size() - xzeros;
    }
    return static_cast<double>(result);
}

std::shared_ptr<ioh::problem::FunctionalConstraint<int>> make_constraint(std::function<double(std::vector<int>)> fn){
    return std::make_shared<ioh::problem::FunctionalConstraint<int>>(
        fn, 1.0, 1, ioh::problem::constraint::Enforced::HIDDEN);
}

std::shared_ptr<ioh::problem::IntegerSingleObjective> get_problem(const int id, const int dimension)
{
    std::shared_ptr<ioh::problem::IntegerSingleObjective> problem;
    switch (id)
    {
    case 1:
    {
        auto &oneminmaxfactory = ioh::problem::ProblemRegistry<ioh::problem::PBO>::instance();
        problem = oneminmaxfactory.create(id, 1, dimension);
        problem->add_constraint(make_constraint(onemax_constraint));
        if (pre_pareto_y.size() != 0) {pre_pareto_y.clear();}
        for (size_t i = 0; i <= dimension; ++i) {
            std::array<double, 2> y = {static_cast<double>(i), static_cast<double>(dimension- i)};
            pre_pareto_y.push_back(y);

        }
        break;
    }
    case 2:
    {
        auto &leadingonesfactory = ioh::problem::ProblemRegistry<ioh::problem::PBO>::instance();
        problem = leadingonesfactory.create(id, 1, dimension);
        problem->add_constraint(make_constraint(leadingones_constraint));
        
        if (pre_pareto_y.size() != 0) {pre_pareto_y.clear();}
        std::array<double, 2> y = {0.0, static_cast<double>(dimension)};
        pre_pareto_y.push_back(y);
        for (size_t i = 1; i <= dimension; ++i) {
            y[0] = static_cast<double>(i);
            y[1] = static_cast<double>(dimension -i);
            pre_pareto_y.push_back(y);
        }
        // y[0] = static_cast<double>(dimension);
        // y[1] = 0.0;
        // pre_pareto_y.push_back(y);
        break;
    }
    case 3:
    {
        ioh::problem::wrap_function<int, double>(&onejump,                           // the new function
                                                 "onejumpzerojump",                  // name of the new function
                                                 ioh::common::OptimizationType::MAX, // optimization type
                                                 0,                                  // lowerbound
                                                 1                                   // upperbound
        );
        auto &jumpfactory = ioh::problem::ProblemRegistry<ioh::problem::IntegerSingleObjective>::instance();
        problem = jumpfactory.create("onejumpzerojump", 1, dimension);
        problem->add_constraint(make_constraint(jump_constraint));

        int k = 2;
        if (pre_pareto_y.size() != 0) {pre_pareto_y.clear();}
        std::array<double, 2> y = {0.0 + k, static_cast<double>(dimension) + k};
        pre_pareto_y.push_back(y);
        for (size_t i = 2 * k; i <= dimension; ++i) {
            y[0] = static_cast<double>(i);
            y[1] = static_cast<double>(2* k + dimension - i);
            pre_pareto_y.push_back(y);
        }
        y[0] = static_cast<double>(dimension) + k;
        y[1] = 0.0 + k;
        pre_pareto_y.push_back(y);
        break;
    }
    case 4:
    {
        ioh::problem::wrap_function<int, double>(&COCZ,                           // the new function
                                                 "COCZ",                  // name of the new function
                                                 ioh::common::OptimizationType::MAX, // optimization type
                                                 0,                                  // lowerbound
                                                 1                                   // upperbound
        );
        auto &coczfactory = ioh::problem::ProblemRegistry<ioh::problem::IntegerSingleObjective>::instance();
        problem = coczfactory.create("COCZ", 1, dimension);
        problem->add_constraint(make_constraint(COCZ_constraint));

        if (pre_pareto_y.size() != 0) {pre_pareto_y.clear();}
        std::array<double, 2> y = {0,0};
        for (size_t i = (dimension /2) ; i <= dimension; ++i) {
            y[0] = static_cast<double>(i);
            y[1] = dimension - i  + (dimension / 2);
            pre_pareto_y.push_back(y);
        }
        break;
    }
    default:
        std::cerr << "unknown problem id" << std::endl;
    }
    return problem;
}
