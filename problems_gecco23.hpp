#include "ioh.hpp"

extern int dimension;

double onemax_constraint(const std::vector<int>& x) {
    return (double)x.size() - std::accumulate(x.begin(), x.end(), 0.0);
}

double leadingones_constraint(const std::vector<int>& x) {
    auto result = 0;
    for (auto i = x.size() -1; i > 0; --i)   
        if (x[i] == 0)
            result++;
        else
            break;

    return static_cast<double>(result);
}

double onejump(const std::vector<int> &x)
{ 
    int k = 2;
    int xones = std::accumulate(x.begin(), x.end(), 0.0);
    double result;

    if(xones != x.size() && xones <= x.size() - k) {
        result = xones + k;
    }
    else {
        result = x.size() - xones;
    }
    return static_cast<double>(result);
}

double jump_constraint(const std::vector<int> &x)
{ 
    int k = 2;
    int xzeros = x.size() - std::accumulate(x.begin(), x.end(), 0.0);
    double result;

    if(xzeros != x.size() && xzeros <= x.size() - k) {
        result = xzeros + k;
    }
    else {
        result = x.size() - xzeros;
    }
    return static_cast<double>(result);
}

auto minmax = std::make_shared<ioh::problem::FunctionalConstraint<int>>(
         onemax_constraint, 1.0, 1, ioh::problem::constraint::Enforced::HIDDEN);

auto trailing_zeros = std::make_shared<ioh::problem::FunctionalConstraint<int>>(
        leadingones_constraint, 1.0, 1, ioh::problem::constraint::Enforced::HIDDEN);
auto zerojump = std::make_shared<ioh::problem::FunctionalConstraint<int>>(
        jump_constraint, 1.0, 1, ioh::problem::constraint::Enforced::HIDDEN);
        
    
std::shared_ptr<ioh::problem::Integer> get_problem(const int id, const int dimension) {
    std::shared_ptr<ioh::problem::Integer> problem;
    switch(id)
    {   
        case 1: {
            auto &oneminmaxfactory = ioh::problem::ProblemRegistry<ioh::problem::PBO>::instance();
            problem = oneminmaxfactory.create(id, 1, dimension);
            problem->add_constraint(minmax);
            break; }
        case 2: {
            auto &leadingonesfactory = ioh::problem::ProblemRegistry<ioh::problem::PBO>::instance();   
            problem = leadingonesfactory.create(id, 1, dimension);
            problem->add_constraint(trailing_zeros);
            break;}
        case 3:{
            ioh::problem::wrap_function<int>(&onejump,  // the new function
                                      "onejumpzerojump", // name of the new function
                                      ioh::common::OptimizationType::MAX, // optimization type
                                      0,  // lowerbound  
                                      1  // upperbound
                                    );
            auto &jumpfactory = ioh::problem::ProblemRegistry<ioh::problem::Integer>::instance();
            problem = jumpfactory.create("onejumpzerojump",1,dimension);
            problem->add_constraint(zerojump);
            break;}
        default:
            std::cerr << "unknown problem id" << std::endl;
    }
    return problem;
}

